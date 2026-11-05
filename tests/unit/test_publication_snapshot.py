"""Distribution snapshots must not mix changing checkouts with historical study code."""
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import warnings
import zipfile
import pytest

spec=importlib.util.spec_from_file_location('publication_snapshot',Path(__file__).resolve().parents[2]/'tools/publication_snapshot.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def git(path,*args):
    return subprocess.check_output(['git','-C',str(path),'-c','user.name=Snapshot test',
        '-c','user.email=snapshot@example.invalid','-c','commit.gpgsign=false','-c','core.hooksPath=/dev/null',*args],text=True).strip()


def repository(path,name,text):
    path.mkdir();git(path,'init','--initial-branch=main')
    file=path/name;file.parent.mkdir(parents=True,exist_ok=True);file.write_text(text)
    git(path,'add','.');git(path,'commit','-m','fixture')
    return git(path,'rev-parse','HEAD')


def fixture(tmp_path,bad=None):
    code,paper,legacy=tmp_path/'code',tmp_path/'paper',tmp_path/'legacy'
    c=repository(code,'src/metag/core.py','pinned code\n')
    p=repository(paper,'sections/abstract.tex','pinned manuscript\n')
    legacy.mkdir();codezip=legacy/'MetaG_code_and_evidence.zip'
    entries=[('MetaG_submission/MetaG/src/metag/core.py',b'historical runtime\n'),
             ('MetaG_submission/MetaG/manuscript/sections/abstract.tex',b'old paper\n'),
             ('MetaG_submission/analysis/old-study/raw.json',b'{"energy":123}'),
             ('MetaG_submission/verification/old-review.json',b'{"old":true}')]
    if bad=='unsafe':entries.append(('MetaG_submission/../../escape',b'no'))
    if bad=='duplicate':entries.append(entries[-1])
    with warnings.catch_warnings():
        warnings.simplefilter('ignore',UserWarning)
        with zipfile.ZipFile(codezip,'w') as z:
            for name,data in entries:z.writestr(name,data)
    for name in module.COMPANIONS:
        with zipfile.ZipFile(legacy/name,'w') as z:z.writestr('raw.json',b'{"retained":true}')
    manifest=dict(public_source_commit='old-code',paper_commit='old-paper',runtime_source_sha256='0'*64,
        files={f.name:dict(bytes=f.stat().st_size,sha256=module.digest(f)) for f in legacy.glob('*.zip')})
    (legacy/'UPLOAD_MANIFEST.json').write_text(json.dumps(manifest))
    return code,paper,c,p,legacy


def test_snapshot_uses_named_commits_and_namespaces_old_runtime(tmp_path):
    code,paper,c,p,legacy=fixture(tmp_path)
    (code/'src/metag/core.py').write_text('newer committed code\n')
    git(code,'add','.');git(code,'commit','-m','newer head')
    (paper/'sections/abstract.tex').write_text('uncommitted collaborator work\n')
    (code/'untracked.txt').write_text('not part of the pinned source\n')
    output=tmp_path/'snapshot';r=module.build(code,paper,c,p,legacy,output)
    with zipfile.ZipFile(output/'MetaG_code_and_evidence.zip') as z:
        assert z.read('MetaG_submission/MetaG/src/metag/core.py')==b'pinned code\n'
        assert z.read('MetaG_submission/MetaG/manuscript/sections/abstract.tex')==b'pinned manuscript\n'
        assert z.read('MetaG_submission/historical_snapshots/review_20261104/src/metag/core.py')==b'historical runtime\n'
        assert z.read('MetaG_submission/analysis/old-study/raw.json')==b'{"energy":123}'
        assert not any('untracked.txt' in name for name in z.namelist())
    for name in module.COMPANIONS:assert (output/name).read_bytes()==(legacy/name).read_bytes()
    assert r['source_map']['code_commit']==c and r['source_map']['paper_commit']==p
    assert r['archive_member_readback_verified'] and not r['submission_ready']
    assert (paper/'sections/abstract.tex').read_text()=='uncommitted collaborator work\n'
    with pytest.raises(FileExistsError):module.build(code,paper,c,p,legacy,output)


def test_modified_legacy_archive_is_rejected_before_creating_output(tmp_path):
    args=fixture(tmp_path);legacy=args[-1]
    with (legacy/module.COMPANIONS[0]).open('ab') as stream:stream.write(b'changed')
    output=tmp_path/'snapshot'
    with pytest.raises(ValueError,match='differs from its declared manifest'):
        module.build(*args,output)
    assert not output.exists()


@pytest.mark.parametrize('kind',['unsafe','duplicate'])
def test_unsafe_or_duplicate_legacy_members_cannot_become_a_valid_snapshot(tmp_path,kind):
    args=fixture(tmp_path,kind);output=tmp_path/'snapshot'
    with pytest.raises(ValueError,match='Unsafe archive name|Duplicate legacy member'):
        module.build(*args,output)
    assert not (output/'SNAPSHOT_MANIFEST.json').exists()
    assert not (tmp_path.parent/'escape').exists()
