"""Controller completion and molecular arithmetic are separate from file presence."""
import copy
import hashlib
import json
from pathlib import Path
import pytest
from benchmark.operational_records import authenticate,execution_state

MODEL={k:'0'*64 for k in ['checkpoint_sha256','references_sha256','inference_sha256']}
PRODUCER=dict(model=MODEL,cacheable=True)


def test_child_only_success_cannot_be_promoted_and_resource_interruptions_are_censored():
    assert not execution_state(dict(status='molecular_record_complete'))['terminal']
    base=dict(status='molecular_record_complete',controller_execution=dict(exit_code=0,timed_out=False))
    assert execution_state(base)['opportunity_complete']
    for key,value in [('common_deadline_censored',True),('external_GPU_contention',[{'pid':1}])]:
        row=copy.deepcopy(base);row['controller_execution'][key]=value
        assert not execution_state(row)['opportunity_complete']
    quality=dict(status='workflow_no_estimate',controller_execution=dict(exit_code=0,timed_out=False))
    assert execution_state(quality)['opportunity_complete']
    with pytest.raises(ValueError,match='Unknown terminal'):execution_state(dict(status='unrecognized',controller_execution={}))


def fixture(tmp_path):
    slot='fixture';path=tmp_path/'cache'/slot/'actual_record.json';path.parent.mkdir(parents=True)
    geometry=dict(schema='synthetic-test-fixture',symbols=['H','H'],positions_angstrom=[[0.,0.,0.],[.74,0.,0.]],charge=0,spin_multiplicity=1)
    geometry['sha256']=hashlib.sha256(json.dumps(geometry,sort_keys=True,allow_nan=False).encode()).hexdigest()
    minimum=dict(electronic_kj=-100.,solvation_kj=-2.,thermal_kj=-3.,total_kj=-105.,geometry=geometry)
    raw=dict(smi='[H][H]',q=0,method='implicit',G=-105.,sigma=.1,settings=dict(thermal_validation_max=96,producer=PRODUCER),meta=dict(thermal=dict(ensemble=dict(validated_minima=[minimum],n_reference_thermal_approximation=0))))
    path.write_text(json.dumps(raw))
    task=dict(slot=slot,kind='species',source_input_key='0|[H][H]',charge=0,smiles='[H][H]',spin=1,n_atoms=2)
    result=dict(status='molecular_record_complete',cache_path=str(path),cache_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),energy_kj=-105.,sampling_sigma_kj=.1,producer=PRODUCER,controller_execution=dict(exit_code=0,timed_out=False))
    return path,task,result,raw


def test_reconstructs_the_complete_partition_and_rejects_a_changed_value_with_updated_hashes(tmp_path):
    path,task,result,raw=fixture(tmp_path)
    _,audit=authenticate(task,result,tmp_path,dict(thermal_validation_max=96),MODEL)
    assert audit['n_validated_minima']==1 and audit['maximum_decimal_partition_error_kj']<1e-7
    raw['G']+=1.;result['energy_kj']+=1.;path.write_text(json.dumps(raw));result['cache_sha256']=hashlib.sha256(path.read_bytes()).hexdigest()
    with pytest.raises(ValueError,match='partition mismatch'):authenticate(task,result,tmp_path,dict(thermal_validation_max=96),MODEL)


def test_wrong_spin_is_not_an_interchangeable_cache_record(tmp_path):
    path,task,result,raw=fixture(tmp_path);raw['settings']['spin']=3;path.write_text(json.dumps(raw));result['cache_sha256']=hashlib.sha256(path.read_bytes()).hexdigest()
    with pytest.raises(ValueError,match='spin mismatch'):authenticate(task,result,tmp_path,dict(thermal_validation_max=96),MODEL)


def test_original_absolute_cache_path_cannot_supply_an_outside_fallback(tmp_path):
    path,task,result,raw=fixture(tmp_path)
    other=tmp_path/'separate_phase';other.mkdir()
    with pytest.raises(ValueError,match='path/hash'):authenticate(task,result,other,dict(thermal_validation_max=96),MODEL)
