"""Build a commit-pinned publication review package without relabelling study sources."""
from __future__ import annotations
import argparse
import hashlib
import io
import json
import shutil
import subprocess
import tarfile
from pathlib import Path, PurePosixPath
import zipfile

PREFIX='MetaG_submission/'
SOURCE_PREFIX=PREFIX+'MetaG/'
PAPER_PREFIX=SOURCE_PREFIX+'manuscript/'
COMPANIONS=('MetaG_current_study.zip','MetaG_historical_selection.zip','MetaG_dgbyg_comparison.zip')
COMPRESSED={'.zip','.gz','.zst','.xz','.bz2','.png','.jpg','.jpeg','.npz','.pdf','.whl'}


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for data in iter(lambda:stream.read(1024*1024),b''):h.update(data)
    return h.hexdigest()


def safe_name(name):
    p=PurePosixPath(name)
    if not name or p.is_absolute() or '..' in p.parts or '\\' in name or any(':' in part for part in p.parts):
        raise ValueError(f'Unsafe archive name: {name}')
    if p.as_posix()!=name.rstrip('/'):
        raise ValueError(f'Noncanonical archive name: {name}')
    return p


def commit(repo, revision):
    return subprocess.check_output(['git','-C',str(repo),'rev-parse','--verify',revision+'^{commit}'],text=True).strip()


def git_members(repo, revision):
    # Capture only the named committed tree, never a working-tree glob.
    blob=subprocess.check_output(['git','-C',str(repo),'archive','--format=tar',revision])
    with tarfile.open(fileobj=io.BytesIO(blob),mode='r:') as archive:
        for entry in archive:
            safe_name(entry.name)
            if entry.isdir():continue
            if not entry.isfile():raise ValueError('Nonregular committed member: '+entry.name)
            yield entry.name,archive.extractfile(entry).read(),entry.mode


class Writer:
    def __init__(self,path):
        self.archive=zipfile.ZipFile(path,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=1)
        self.members={}
    def add(self,name,data,origin,mode=0o644):
        safe_name(name)
        if name in self.members:raise ValueError('Duplicate output member: '+name)
        info=zipfile.ZipInfo(name,date_time=(1980,1,1,0,0,0))
        info.external_attr=(0o100000 | (mode & 0o777))<<16
        info.compress_type=zipfile.ZIP_STORED if PurePosixPath(name).suffix.lower() in COMPRESSED else zipfile.ZIP_DEFLATED
        self.archive.writestr(info,data,compresslevel=1)
        self.members[name]=dict(bytes=len(data),sha256=hashlib.sha256(data).hexdigest(),origin=origin,mode=mode & 0o777)
    def close(self):self.archive.close()


def build(code_repo,paper_repo,code_revision,paper_revision,legacy_dir,output):
    code_repo,paper_repo,legacy_dir,output=map(Path,[code_repo,paper_repo,legacy_dir,output])
    code_sha,paper_sha=commit(code_repo,code_revision),commit(paper_repo,paper_revision)
    old=json.loads((legacy_dir/'UPLOAD_MANIFEST.json').read_text())
    for name in ['MetaG_code_and_evidence.zip',*COMPANIONS]:
        expected=old['files'][name];path=legacy_dir/name
        if path.stat().st_size!=expected['bytes'] or digest(path)!=expected['sha256']:
            raise ValueError('Legacy archive differs from its declared manifest: '+name)
    if output.exists():raise FileExistsError('Preserve existing snapshots; choose a new output directory')
    output.mkdir(parents=True)
    archive_path=output/'MetaG_code_and_evidence.zip';writer=Writer(archive_path)
    code_count=paper_count=legacy_count=0
    try:
        for name,data,mode in git_members(code_repo,code_sha):
            if name.startswith('manuscript/'):
                raise ValueError('Manuscript must have its own pinned repository')
            if PurePosixPath(name).suffix.lower() in {'.pt','.ckpt','.safetensors'}:
                raise ValueError('Model weights are outside this source-package scope')
            writer.add(SOURCE_PREFIX+name,data,'current_code:'+code_sha,mode);code_count+=1
        for name,data,mode in git_members(paper_repo,paper_sha):
            writer.add(PAPER_PREFIX+name,data,'manuscript:'+paper_sha,mode);paper_count+=1
        with zipfile.ZipFile(legacy_dir/'MetaG_code_and_evidence.zip') as source:
            names=source.namelist()
            if len(names)!=len(set(names)):raise ValueError('Duplicate legacy member')
            for entry in source.infolist():
                name=entry.filename;safe_name(name)
                if entry.is_dir():continue
                if (entry.external_attr>>16)&0o170000==0o120000:
                    raise ValueError('Symlink in legacy evidence')
                if not name.startswith(PREFIX):raise ValueError('Unknown legacy archive layout')
                if name.startswith(SOURCE_PREFIX+'src/'):
                    target=PREFIX+'historical_snapshots/review_20261104/'+name[len(SOURCE_PREFIX):]
                elif name.startswith(SOURCE_PREFIX):
                    continue  # Current committed code/paper replace the old checkout completely.
                elif name.startswith(PREFIX+'distributions/'):
                    target=PREFIX+'historical_snapshots/review_20261104/'+name[len(PREFIX):]
                elif name.startswith(PREFIX+'verification/'):
                    target=PREFIX+'verification/history/review_20261104/'+name[len(PREFIX+'verification/'):]
                elif len(PurePosixPath(name).parts)==2:
                    target=PREFIX+'verification/history/review_20261104/'+PurePosixPath(name).name
                else:
                    target=name  # Historical raw evidence keeps its original study-relative layout.
                writer.add(target,source.read(entry),'retained_legacy_evidence',(entry.external_attr>>16)&0o777);legacy_count+=1
        mapping=dict(schema='metag-publication-source-map-v1',code_commit=code_sha,paper_commit=paper_sha,
            historical_parent_code_commit=old['public_source_commit'],historical_parent_paper_commit=old['paper_commit'],
            historical_parent_runtime_source_sha256=old['runtime_source_sha256'],
            current_code_scope='Current development source; not the source of every historical molecular calculation.',
            study_source_scope='Each retained study keeps its own source/settings/producer identities. No scientific result is relabelled to the current commit.',
            independent_experimental_benchmark_complete=False,production_reference_policy_finalized=False,
            author_metadata_complete=False,submission_ready=False)
        writer.add(PREFIX+'SOURCE_MAP.json',(json.dumps(mapping,indent=2)+'\n').encode(),'snapshot_metadata')
        writer.add(PREFIX+'README.md',(
            '# MetaG publication review snapshot\n\n'
            f'Code commit: `{code_sha}`. Manuscript commit: `{paper_sha}`.\n\n'
            'This is a review snapshot, not a submission-ready release. The independent experimental benchmark, '
            'final supported production reference policy, representative numerical ModelSEED validation and '
            'author-final metadata remain incomplete. See SOURCE_MAP.json and the manuscript limitations.\n\n'
            'Extract the three unchanged companion archives into the same parent. Their legacy study layouts '
            'remain intact. Current code is under MetaG/; the earlier review runtime is separately retained '
            'under historical_snapshots/. Every newer control retains its own frozen source in its evidence packet.\n\n'
            'Use the declared CPU audit environment to run MetaG/manuscript/tools/check_*.py. '
            'Reassembly of retained records is separate from new quantum calculations. Model weights are not included.\n'
        ).encode(),'snapshot_metadata')
    finally:
        writer.close()
    for name in COMPANIONS:
        shutil.copyfile(legacy_dir/name,output/name)
        assert digest(output/name)==old['files'][name]['sha256']
    manuscript=Writer(output/'MetaG_manuscript_source.zip')
    try:
        for name,data,mode in git_members(paper_repo,paper_sha):manuscript.add(name,data,'manuscript:'+paper_sha,mode)
    finally:manuscript.close()
    # Read back every new archive member, including bytes and executable mode.
    for path,expected in [(archive_path,writer.members),(output/'MetaG_manuscript_source.zip',manuscript.members)]:
        with zipfile.ZipFile(path) as archive:
            assert set(archive.namelist())==set(expected)
            for name,row in expected.items():
                data=archive.read(name)
                assert len(data)==row['bytes'] and hashlib.sha256(data).hexdigest()==row['sha256']
                assert ((archive.getinfo(name).external_attr>>16)&0o777)==row['mode']
    result=dict(schema='metag-publication-review-snapshot-v1',source_map=mapping,
        inventory=dict(current_code_files=code_count,current_manuscript_files=paper_count,retained_legacy_files=legacy_count),
        archive_members=writer.members,manuscript_archive_members=manuscript.members,
        files={p.name:dict(bytes=p.stat().st_size,sha256=digest(p)) for p in sorted(output.glob('*.zip'))},
        legacy_companions_byte_identical=True,archive_member_readback_verified=True,
        manuscript_build_and_checks='Required separately on this exact manuscript commit; no PDF is accepted from an unverified working tree.',
        scientific_validation_complete=False,submission_ready=False)
    (output/'SNAPSHOT_MANIFEST.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for arg in ['code-repo','paper-repo','code-commit','paper-commit','legacy-dir','output']:
        parser.add_argument('--'+arg,required=True)
    args=parser.parse_args()
    result=build(args.code_repo,args.paper_repo,args.code_commit,args.paper_commit,args.legacy_dir,args.output)
    print(json.dumps(dict(inventory=result['inventory'],source_map=result['source_map'],files=result['files']),indent=2))


if __name__=='__main__':main()
