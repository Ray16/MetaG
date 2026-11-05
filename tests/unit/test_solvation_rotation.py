"""Keep native failures in the panel and distinguish precision from Gibbs accuracy."""
import importlib.util
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[2]
spec=importlib.util.spec_from_file_location('solvation_rotation_audit',ROOT/'benchmark/solvation_rotation.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
PACKET=ROOT/'benchmark/solvation_rotation_validation'


def test_every_requested_pair_retains_success_or_native_failure():
    report=module.audit(PACKET)
    assert (report['requested_pairs'],report['successful_pairs'],report['failed_pairs'],report['unattempted_pairs'])==(432,426,6,0)
    assert report['actual_native_executions']==858
    assert report['unattempted_downstream_water_calls']==6
    assert not report['full_frozen_complete_panel_endpoint_available']
    assert not report['primary_analysis_gate_weakened']


def test_dense_grid_improves_complete_subcases_without_filling_failed_values():
    report=module.audit(PACKET);rows=report['complete_subcase_diagnostics']
    assert rows['default|default|full_red']['orientation_range_kj']==pytest.approx(.4861097038576645)
    assert rows['extreme|default|full_red']['orientation_range_kj']<.001
    assert rows['extreme|default|charged_ox']['orientation_range_kj'] is None
    assert report['cycle_diagnostics']['charge|extreme|default']['status']=='withheld_due_failed_native_pairs'
    assert not report['production_numerical_policy_adopted']
    assert report['new_experimental_references_admitted']==0


def test_original_complete_panel_gate_still_refuses_without_native_execution():
    report=module.strict_endpoint_replay(PACKET)
    assert report['original_gate_unchanged'] and not report['quantum_execution']
    assert report['status']=='original_complete_panel_endpoint_withheld'
