"""The release-source checker must fail on missing or changed evidence."""
import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

SCRIPT=Path(__file__).resolve().parents[2]/'benchmark/check_manuscript_source_manifest.py'
spec=importlib.util.spec_from_file_location('check_manuscript_source_manifest',SCRIPT)
checker=importlib.util.module_from_spec(spec);spec.loader.exec_module(checker)


def make_manifest(tmp_path, rows):
    manuscript=tmp_path/'manuscript';(manuscript/'evidence').mkdir(parents=True)
    for name in checker.MANIFESTS:
        (manuscript/'evidence'/name).write_text(json.dumps({'sources':rows}))
    for name in checker.MAPPED_MANIFESTS:
        (manuscript/'evidence'/name).write_text(json.dumps({'files':{}}))
    (manuscript/'evidence'/checker.ROUTE_MANIFEST).write_text(json.dumps({'files':{}}))
    (manuscript/'evidence'/checker.GEOMETRY_MANIFEST).write_text(json.dumps({'sources':[],'records':{'component_replay':{'sources_sha256':{}},'exchange_replay':{'sources_sha256':{}}}}))
    return manuscript


def test_presence_and_hash_are_required(tmp_path):
    source=tmp_path/'analysis'/'result.json';source.parent.mkdir();source.write_text('original')
    digest=hashlib.sha256(source.read_bytes()).hexdigest()
    manuscript=make_manifest(tmp_path,[{'path':'analysis/result.json','sha256':digest}])
    assert checker.check_sources(manuscript)['verified']==1
    source.write_text('changed')
    assert checker.check_sources(manuscript)['changed']==['analysis/result.json']
    source.unlink()
    assert checker.check_sources(manuscript)['missing']==['analysis/result.json']


def test_path_traversal_is_rejected(tmp_path):
    manuscript=make_manifest(tmp_path,[{'path':'../other','sha256':'0'*64}])
    with pytest.raises(ValueError,match='Unsafe'):
        checker.check_sources(manuscript)


def test_legacy_physical_path_resolves_inside_staged_root(tmp_path):
    source=tmp_path/'analysis'/'physical.json';source.parent.mkdir();source.write_text('physical')
    digest=hashlib.sha256(source.read_bytes()).hexdigest();manuscript=make_manifest(tmp_path,[])
    (manuscript/'evidence'/checker.PHYSICAL_MANIFEST).write_text(json.dumps({'files':{'physical':{'path':'/old/root/MetaG_JCIM/analysis/physical.json','sha256':digest}}}))
    assert checker.check_sources(manuscript)['verified']==1
