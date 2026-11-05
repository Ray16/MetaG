"""Reconstruct the complete recorded bank and keep numerical and accuracy claims separate."""
import importlib.util
import json
from pathlib import Path
import zipfile

import pytest

ROOT=Path(__file__).resolve().parents[2]
spec=importlib.util.spec_from_file_location('global_amino_audit',ROOT/'benchmark/global_amino_reference.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
PACKET=ROOT/'benchmark/global_amino_reference_validation'


def test_full_raw_bank_and_equal_target_adverse_result():
    report=module.audit(PACKET)
    assert report['authenticated_raw_records']==596 and report['common_labels']==112
    assert report['maximum_raw_partition_difference_kj']<1e-7
    assert report['independent_common_MAE_kj']['baseline']==pytest.approx(16.678183744158485)
    assert report['independent_common_MAE_kj']['candidate']==pytest.approx(17.107155011067544)
    assert not report['production_reference_policy_adopted']


def test_cold_reassembly_reproduces_all_requests_without_quantum_imports():
    report=module.replay(PACKET)
    assert report['requests']==364 and report['reproduced_estimates']==148
    assert report['quantum_imports']==0


def test_independent_decimals_and_projected_hessians():
    report=module.replay(PACKET,physical=True)
    assert report['estimated_reactions_verified']==148 and report['n_source_states']==6
    for row in report['independent_Hessian_controls'].values():
        assert row['n_strong_imag']==0 and row['force_max_eV_A']<=.001
        assert abs(row['difference_kj'])<.0001


def test_original_refusal_is_preserved_and_recovery_is_a_complete_record():
    with zipfile.ZipFile(PACKET/'numerical_controls.zip') as z:
        old=json.loads(z.read('quantum_extension/results/alanine_S_parent.json'))
        new=json.loads(z.read('recovery_calculation/results/alanine_S_parent.json'))
    assert old['status']=='workflow_no_estimate' and old['energy_kj'] is None
    assert new['status']=='molecular_record_complete' and new['energy_kj'] is not None
    catalog=json.loads((PACKET/'catalog596.json').read_text())
    row=next(r for r in catalog['records'] if r['raw_profile']=='curvature_recovery')
    assert row['sha256']==new['cache_sha256']
    with zipfile.ZipFile(PACKET/'raw_records.zip') as z:raw=json.loads(z.read(row['name']))
    e=raw['meta']['thermal']['ensemble']
    assert e['n_candidates_validated']==len(e['validated_minima'])==1
    assert e['n_reference_thermal_approximation']==0
    trial=e['optimizer_trials'][-1]
    assert trial['strict_thermal_accepted'] and trial['max_force_eV_A']<=.001
