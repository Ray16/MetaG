"""Numerical audit remains separate from experimental validation and coverage claims."""
from pathlib import Path
import hashlib
import json
import shutil
import zipfile

import pytest
from benchmark.boundary_reference import audit,replay

PACKET=Path(__file__).resolve().parents[2]/'benchmark/boundary_reference_validation'


def test_every_raw_partition_reaction_and_shared_uncertainty_reconstructs():
    result=audit(PACKET)
    assert result['raw_molecular_records']==25
    assert result['availability']['primary']=={'full':5,'anchor14':5,'anchor12and14':8}
    assert result['availability']['augmented']=={'full':11,'anchor14':11,'anchor12and14':12}
    assert result['augmented_invariance_checks']==136
    assert set(result['augmented_cycles_per_arm'].values())=={1}
    assert result['curvature_controls']['n_fixed_geometries']==4
    assert result['curvature_controls']['maximum_retained_RRHO_difference_kj']<.0001
    assert result['independent_experimental_labels']==0 and not result['production_policy_adopted']


def test_frozen_both_profile_replay_requires_no_quantum_runtime_or_network():
    result=replay(PACKET)
    assert result['requested_per_arm']==12
    assert result['availability']['augmented']['anchor12and14']==12
    assert result['quantum_imports_forbidden'] and result['network_forbidden']


def test_altered_prediction_fails_even_when_container_checksums_are_updated(tmp_path):
    packet=tmp_path/'packet';shutil.copytree(PACKET,packet,ignore=shutil.ignore_patterns('__pycache__'))
    path=packet/'assemblies.zip'
    with zipfile.ZipFile(path) as archive:entries={n:archive.read(n) for n in archive.namelist()}
    name='augmented/anchor12and14.json';doc=json.loads(entries[name])
    doc['results']['rxn47022']['result']['dG_unrounded']+=1.
    entries[name]=(json.dumps(doc,indent=2)+'\n').encode()
    with zipfile.ZipFile(path,'w',compression=zipfile.ZIP_DEFLATED) as archive:
        for name,data in entries.items():archive.writestr(name,data)
    m=json.loads((packet/'manifest.json').read_text())
    m['files']['assemblies.zip']=dict(bytes=path.stat().st_size,sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    (packet/'manifest.json').write_text(json.dumps(m,indent=2)+'\n')
    with pytest.raises(AssertionError):audit(packet)
