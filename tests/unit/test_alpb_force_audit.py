"""Verify retained derivative evidence and reject corrupted native records."""
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/'benchmark/audit_alpb_forces.py'
spec=importlib.util.spec_from_file_location('alpb_force_auditor',SOURCE)
audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)
PACKET=SOURCE.parent/'aqueous_force_validation'

def fixture():
    saved=json.loads((PACKET/'result.json').read_text());case=saved['cases'][0]
    with gzip.open(PACKET/'ethanol.raw.jsonl.gz','rt') as stream:record=json.loads(next(stream))
    return record,case,json.loads((PACKET/'protocol.json').read_text())

def test_all_refined_grid_derivatives_reconstruct():
    report=audit.audit(PACKET)
    assert report['all_cases_passed'] and report['n_cases']==4
    assert report['n_native_executions']==800

def test_default_grid_failure_is_preserved():
    report=audit.audit(PACKET/'default_grid_trial')
    assert report['n_cases']==1 and report['n_native_executions']==218
    assert not report['all_cases_passed']
    assert report['cases'][0]['net_torque_eV']>.005

def test_packet_and_decompressed_native_bytes_are_authenticated():
    manifest=json.loads((PACKET/'manifest.json').read_text())
    for name,digest in manifest['files'].items():
        assert hashlib.sha256((PACKET/name).read_bytes()).hexdigest()==digest
    for name,digest in manifest['raw_uncompressed_sha256'].items():
        assert hashlib.sha256(gzip.decompress((PACKET/(name+'.gz')).read_bytes())).hexdigest()==digest

@pytest.mark.parametrize('mutation',[
    lambda r:r.update(charge=1),
    lambda r:r.update(solvent=True),
    lambda r:r.update(returncode=1),
    lambda r:r.update(energy_Eh=r['energy_Eh']+.1),
    lambda r:r.update(gradient_requested=False),
    lambda r:r.update(input_xyz=r['input_xyz'].replace('C ','N ',1)),
    lambda r:r['argv'].__setitem__(r['argv'].index('--uhf')+1,'2'),
])
def test_invalid_record_is_rejected(mutation):
    record,case,protocol=fixture();mutation(record)
    with pytest.raises(ValueError):
        audit.check_record(record,case['symbols'],audit.np.array(case['positions_A']),case['case']['charge'],False,True,protocol)

def test_gradient_coordinate_mismatch_is_rejected():
    record,case,protocol=fixture()
    with pytest.raises(ValueError):
        audit.extract_gradient(record['gradient_file'],case['symbols'],audit.np.array(case['positions_A'])+.01,protocol)
