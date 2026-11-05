import copy

import pytest

from benchmark.operational_accounting import dependency_report, summarize_raw


def task(slot):
    return dict(slot=slot, source_input_key='key:' + slot)


def outcome(t, status, **fields):
    return dict(slot=t['slot'], request=copy.deepcopy(t), status=status, **fields)


def test_resource_interruption_cost_is_retained_and_censoring_is_not_qc_failure():
    a, b, c = map(task, ['a', 'b', 'c'])
    results = dict(a=outcome(a, 'external_resource_interruption',
                            execution={'external_GPU_contention': [{'pid': 7}]}, elapsed_seconds=32.3),
                   b=outcome(b, 'execution_timeout', controller_execution={
                       'exit_code': None, 'timed_out': True, 'common_deadline_censored': True},
                             controller_elapsed_seconds=29.4))
    report = summarize_raw([a, b, c], results)
    assert report['observed_controller_seconds'] == pytest.approx(61.7)
    assert report['stage_counts'] == dict(resource_or_common_budget_censoring=2, pending_or_unstarted=1)
    assert report['outcomes_with_observed_controller_time'] == 2
    assert not report['all_declared_opportunities_complete']
    assert report['rows']['c']['observed_controller_seconds'] is None


def test_completed_refusal_and_pending_dependencies_are_reported_together():
    a, b, c = map(task, ['a', 'b', 'c'])
    results = dict(a=outcome(a, 'workflow_no_estimate', controller_execution={'exit_code': 0},
                            diagnostics={'thermal': {'status': 'incomplete'}},
                            controller_elapsed_seconds=80., model_initialization_seconds=11.),
                   b=outcome(b, 'molecular_record_complete', elapsed_seconds=500.))
    report = dependency_report(['key:a', 'key:b', 'key:c'], [a, b, c], results)
    assert report['has_completed_failed_prerequisite']
    assert report['has_pending_prerequisites']
    assert not report['has_resource_or_common_budget_censoring']
    raw = summarize_raw([a, b, c], results)
    assert raw['observed_controller_seconds'] == 80.
    assert raw['observed_backend_construction_seconds'] == 11.
    assert raw['rows']['b']['observed_controller_seconds'] is None
    assert raw['rows']['b']['opportunity']['reason'] == 'controller_has_not_authenticated_completion'
    results['a']['controller_execution']['exit_code'] = 1
    failed_execution = summarize_raw([a, b, c], results)
    assert failed_execution['rows']['a']['stage'] == 'execution_failure_or_full_item_timeout'


def test_unknown_input_and_nonfinite_elapsed_time_fail_closed():
    a = task('a')
    wrong = outcome(task('other'), 'external_resource_interruption', elapsed_seconds=1.)
    with pytest.raises(ValueError, match='declared task'):
        summarize_raw([a], {'a': wrong})
    with pytest.raises(ValueError, match='Undeclared prerequisite'):
        dependency_report(['unowned'], [a], {})
    bad = outcome(a, 'external_resource_interruption', elapsed_seconds=float('nan'))
    with pytest.raises(ValueError, match='observed timing'):
        summarize_raw([a], {'a': bad})
