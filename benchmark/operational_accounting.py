"""Descriptive failure and timing accounting for every declared raw opportunity."""
from collections import Counter
import math

from benchmark.operational_records import execution_state


def _number(value):
    if (not isinstance(value, (int, float)) or isinstance(value, bool)
            or not math.isfinite(value) or value < 0):
        raise ValueError('Invalid observed timing')
    return float(value)


def _index(tasks, outcomes):
    slots = {task['slot']: task for task in tasks}
    keys = {task['source_input_key']: task for task in tasks}
    if len(slots) != len(tasks) or len(keys) != len(tasks):
        raise ValueError('Duplicate declared slot or primitive identity')
    if set(outcomes) - set(slots):
        raise ValueError('Undeclared raw outcome')
    for slot, result in outcomes.items():
        if result['slot'] != slot or result['request'] != slots[slot]:
            raise ValueError('Outcome does not belong to its declared task')
    return slots, keys


def _stage(result, state):
    reason = state['reason']
    if reason == 'workflow_no_estimate':
        execution = result['controller_execution']
        if execution.get('timed_out') or execution.get('exit_code') != 0:
            return 'execution_failure_or_full_item_timeout'
        thermal = result.get('diagnostics', {}).get('thermal', {})
        if thermal.get('status') == 'incomplete':
            return 'minimum_or_thermal_QC_refusal'
        return 'workflow_refusal_unclassified'
    if reason in {'external_resource_interruption', 'common_deadline_censored'}:
        return 'resource_or_common_budget_censoring'
    if not state['terminal']:
        return 'pending_or_unstarted'
    if reason == 'molecular_record_complete':
        return 'completed_record_requires_physical_authentication'
    return 'execution_failure_or_full_item_timeout'


def summarize_raw(tasks, outcomes):
    """Keep censored work distinct from molecular failures and elapsed-time gaps.

    This describes recorded opportunities. It neither authenticates molecular
    energies nor estimates coverage, GPU utilization or reservation occupancy.
    """
    slots, _ = _index(tasks, outcomes)
    rows = {}
    for slot, task in slots.items():
        result = outcomes.get(slot)
        state = execution_state(result)
        row = dict(slot=slot, source_input_key=task['source_input_key'],
                   opportunity=state, stage=_stage(result, state),
                   observed_controller_seconds=None, backend_construction_seconds=None)
        if result is not None and state['terminal']:
            # The resource-interruption branch writes elapsed_seconds before
            # the ordinary controller_execution stamp. Retain this real cost.
            if 'controller_elapsed_seconds' in result:
                row['observed_controller_seconds'] = _number(result['controller_elapsed_seconds'])
                row['timing_field'] = 'controller_elapsed_seconds'
            elif (state['reason'] == 'external_resource_interruption'
                  and 'elapsed_seconds' in result):
                row['observed_controller_seconds'] = _number(result['elapsed_seconds'])
                row['timing_field'] = 'elapsed_seconds_resource_branch'
            if 'model_initialization_seconds' in result:
                row['backend_construction_seconds'] = _number(result['model_initialization_seconds'])
        rows[slot] = row
    recorded = [r['observed_controller_seconds'] for r in rows.values()
                if r['observed_controller_seconds'] is not None]
    construction = [r['backend_construction_seconds'] for r in rows.values()
                    if r['backend_construction_seconds'] is not None]
    return dict(rows=rows, stage_counts=dict(Counter(r['stage'] for r in rows.values())),
                all_declared_opportunities_complete=all(r['opportunity']['opportunity_complete']
                                                        for r in rows.values()),
                observed_controller_seconds=math.fsum(recorded),
                outcomes_with_observed_controller_time=len(recorded),
                terminal_outcomes_without_observed_controller_time=sum(
                    r['opportunity']['terminal'] and r['observed_controller_seconds'] is None
                    for r in rows.values()),
                nonterminal_opportunities=sum(not r['opportunity']['terminal'] for r in rows.values()),
                observed_backend_construction_seconds=math.fsum(construction),
                backend_construction_times_available=len(construction),
                scope='Recorded raw-job elapsed sums, including observed resource interruptions. '
                      'Missing times remain missing. Construction covers imports/loader/backend '
                      'creation; first forward may initialize further. No population coverage, '
                      'physical validity, GPU occupancy or standalone reaction-runtime inference.')


def dependency_report(prerequisites, tasks, outcomes):
    """Retain all failure/pending causes; no dominant-cause zero filling."""
    _, keys = _index(tasks, outcomes)
    rows = []
    for key in sorted(set(prerequisites)):
        if key not in keys:
            raise ValueError('Undeclared prerequisite')
        task = keys[key]
        result = outcomes.get(task['slot'])
        state = execution_state(result)
        rows.append(dict(source_input_key=key, slot=task['slot'],
                         opportunity=state, stage=_stage(result, state)))
    # A complete failed opportunity fixes unavailability under the declared
    # one-opportunity rule, even if other dependencies are still pending.
    failed = [r for r in rows if r['opportunity']['opportunity_complete']
              and r['opportunity']['reason'] != 'molecular_record_complete']
    return dict(prerequisites=rows,
                has_completed_failed_prerequisite=bool(failed),
                has_resource_or_common_budget_censoring=any(
                    r['stage'] == 'resource_or_common_budget_censoring' for r in rows),
                has_pending_prerequisites=any(not r['opportunity']['terminal'] for r in rows),
                molecular_records_still_require_authentication=True)
