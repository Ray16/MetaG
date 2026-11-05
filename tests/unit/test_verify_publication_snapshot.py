"""Portable replay must detect original-workspace reads even when a script catches the error."""
import importlib.util
import json
from pathlib import Path
import sys
import zipfile
import pytest

ROOT=Path(__file__).resolve().parents[2]
def load(name):
    spec=importlib.util.spec_from_file_location(name,ROOT/'tools'/f'{name}.py')
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod
builder=load('publication_snapshot');verifier=load('verify_publication_snapshot')


@pytest.mark.parametrize('attempt_fallback',[False,True])
def test_replay_requires_owned_inputs_and_reports_caught_fallbacks(tmp_path,attempt_fallback):
    original=tmp_path/'original';original.mkdir();(original/'outside.txt').write_text('outside fixture')
    snapshot=tmp_path/'snapshot';snapshot.mkdir()
    script="from pathlib import Path\nassert (Path(__file__).parents[1]/'owned.txt').read_text()=='owned input'\n"
    if attempt_fallback:
        script+=f"try:\n    Path({str(original/'outside.txt')!r}).read_text()\nexcept PermissionError:\n    pass\n"
    w=builder.Writer(snapshot/verifier.ARCHIVES[0])
    w.add('MetaG_submission/MetaG/manuscript/owned.txt',b'owned input','fixture')
    w.add('MetaG_submission/MetaG/manuscript/tools/check_owned.py',script.encode(),'fixture')
    w.close()
    for i,name in enumerate(verifier.ARCHIVES[1:]):
        with zipfile.ZipFile(snapshot/name,'w') as z:z.writestr(f'MetaG_submission/analysis/study{i}/raw.txt','fixture')
    manifest=dict(source_map=dict(code_commit='fixture',paper_commit='fixture'),archive_members=w.members,
        files={name:dict(bytes=(snapshot/name).stat().st_size,sha256=builder.digest(snapshot/name)) for name in verifier.ARCHIVES})
    (snapshot/'SNAPSHOT_MANIFEST.json').write_text(json.dumps(manifest))
    result=verifier.verify(snapshot,tmp_path/'extracted',tmp_path/'logs',[original],sys.executable,30)
    assert result['n_checks']==result['n_passed']==1
    assert result['source_members_unchanged'] and not result['missing_direct_guard_receipts']
    assert result['scientific_accuracy_validated'] is False
    assert result['status']==('failed' if attempt_fallback else 'passed')
    assert bool(result['blocked_fallback_or_network_attempts'])==attempt_fallback


def test_extractor_rejects_a_path_escape_without_writing_it(tmp_path):
    path=tmp_path/'input.zip'
    with zipfile.ZipFile(path,'w') as z:z.writestr('../outside.txt','no')
    with pytest.raises(ValueError,match='Unsafe archive member'):
        verifier.extract(path,tmp_path/'destination')
    assert not (tmp_path/'outside.txt').exists()
