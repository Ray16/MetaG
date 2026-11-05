"""Verify and replay an extracted publication snapshot with original-project reads blocked."""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import signal
import subprocess
import sys
import time
import zipfile

ARCHIVES=('MetaG_code_and_evidence.zip','MetaG_current_study.zip','MetaG_historical_selection.zip','MetaG_dgbyg_comparison.zip')


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()


def extract(path,destination,expected=None):
    with zipfile.ZipFile(path) as z:
        entries=z.infolist()
        if len({e.filename for e in entries})!=len(entries):raise ValueError('Duplicate archive member')
        if sum(e.file_size for e in entries)>5*1024**3:raise ValueError('Archive expansion exceeds review limit')
        if expected is not None and set(z.namelist())!=set(expected):raise ValueError('Member inventory mismatch')
        for e in entries:
            p=PurePosixPath(e.filename)
            if (p.is_absolute() or '..' in p.parts or '\\' in e.filename or any(':' in part for part in p.parts)
                    or p.as_posix()!=e.filename.rstrip('/')):raise ValueError('Unsafe archive member')
            if e.is_dir():continue
            kind=(e.external_attr>>16)&0o170000
            if kind not in {0,0o100000}:raise ValueError('Nonregular archive member')
            target=destination/Path(*p.parts)
            data=z.read(e)
            if expected is not None:
                row=expected[e.filename]
                if len(data)!=row['bytes'] or hashlib.sha256(data).hexdigest()!=row['sha256']:
                    raise ValueError('Member hash mismatch: '+e.filename)
                if ((e.external_attr>>16)&0o777)!=row['mode']:raise ValueError('Member mode mismatch')
            if target.exists():
                if target.is_symlink() or digest(target)!=hashlib.sha256(data).hexdigest():
                    raise ValueError('Conflicting extracted member')
            else:
                target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(data)
                target.chmod((e.external_attr>>16)&0o777 or 0o644)
    return len(entries)


def write_guard(folder,allowed,forbidden,receipts):
    folder.mkdir(parents=True);receipts.mkdir(parents=True)
    body='''import atexit,json,os,sys
from pathlib import Path
allowed=ALLOWED
forbidden=FORBIDDEN
receipt_root=RECEIPTS
blocked=[]
def inside(value,root):return value==root or value.startswith(root+os.sep)
def audit(event,args):
    if event=='open' and not isinstance(args[0],int):
        path=os.path.realpath(os.fsdecode(args[0]))
        if any(inside(path,r) for r in forbidden) and not any(inside(path,r) for r in allowed):
            blocked.append(dict(event=event,path=path));raise PermissionError('Original project read blocked: '+path)
    if event in {'socket.connect','socket.getaddrinfo'}:
        blocked.append(dict(event=event));raise PermissionError('Network access blocked in portable replay')
def finish():
    p=Path(receipt_root)/(str(os.getpid())+'.json')
    p.write_text(json.dumps(dict(pid=os.getpid(),guard_active=True,blocked_attempts=blocked)))
sys.dont_write_bytecode=True
sys.addaudithook(audit)
atexit.register(finish)
'''
    body=body.replace('ALLOWED',repr([str(Path(p).resolve()) for p in [*allowed,folder,receipts]]))
    body=body.replace('FORBIDDEN',repr([str(Path(p).resolve()) for p in forbidden]))
    body=body.replace('RECEIPTS',repr(str(receipts.resolve())))
    path=folder/'sitecustomize.py';path.write_text(body);return path


def verify(snapshot,destination,logs,forbidden_roots,python=sys.executable,timeout=180):
    snapshot,destination,logs=map(lambda p:Path(p).resolve(),[snapshot,destination,logs])
    if destination.exists() or logs.exists():raise FileExistsError('Use new extraction and log directories')
    if not forbidden_roots:raise ValueError('Declare original project roots for the no-fallback check')
    manifest=json.loads((snapshot/'SNAPSHOT_MANIFEST.json').read_text())
    for name in ARCHIVES:
        row=manifest['files'][name];p=snapshot/name
        if p.stat().st_size!=row['bytes'] or digest(p)!=row['sha256']:raise ValueError('Archive hash mismatch: '+name)
    destination.mkdir(parents=True);logs.mkdir(parents=True)
    inventory={name:extract(snapshot/name,destination,manifest['archive_members'] if name==ARCHIVES[0] else None) for name in ARCHIVES}
    code=destination/'MetaG_submission/MetaG';paper=code/'manuscript'
    scripts=sorted((paper/'tools').glob('check_*.py'))
    if not scripts:raise ValueError('No manuscript checks in the snapshot')
    guard=write_guard(logs/'guard',[destination],forbidden_roots,logs/'guard_receipts')
    environment=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',METAG_EVIDENCE_PORTABLE_ONLY='1',
        OMP_NUM_THREADS='1',OMP_THREAD_LIMIT='1',MKL_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',
        PYTHONPATH=os.pathsep.join(map(str,[guard.parent,code/'src',code,paper/'tools'])))
    checks=[]
    for script in scripts:
        started=time.monotonic()
        process=subprocess.Popen([str(python),str(script)],cwd=paper,env=environment,
            stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,start_new_session=True)
        try:
            output,_=process.communicate(timeout=timeout);exit_code=process.returncode
        except subprocess.TimeoutExpired:
            os.killpg(process.pid,signal.SIGTERM)
            try:output,_=process.communicate(timeout=3)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid,signal.SIGKILL);output,_=process.communicate()
            exit_code=124
        (logs/(script.stem+'.log')).write_text(output)
        check=dict(script=script.name,script_sha256=digest(script),exit_code=exit_code,
            elapsed_seconds=time.monotonic()-started,process_id=process.pid)
        checks.append(check)
        (logs/'checks.json').write_text(json.dumps(checks,indent=2)+'\n')
        print(('PASS ' if exit_code==0 else 'FAIL ')+script.name,flush=True)
    # Check the captured source/evidence again after every script has run.
    unchanged=True
    for name,row in manifest['archive_members'].items():
        path=destination/name
        if not path.is_file() or path.is_symlink() or digest(path)!=row['sha256']:
            unchanged=False;break
    receipts=[json.loads(p.read_text()) for p in (logs/'guard_receipts').glob('*.json')]
    guarded={r['pid'] for r in receipts if r['guard_active']}
    missing=[r['process_id'] for r in checks if r['process_id'] not in guarded]
    blocked=[dict(pid=r['pid'],attempt=a) for r in receipts for a in r['blocked_attempts']]
    result=dict(schema='metag-portable-publication-replay-v1',source_map=manifest['source_map'],
        snapshot_manifest_sha256=digest(snapshot/'SNAPSHOT_MANIFEST.json'),archive_entry_counts=inventory,
        n_checks=len(checks),n_passed=sum(c['exit_code']==0 for c in checks),source_members_unchanged=unchanged,
        guard_source_sha256=digest(guard),guarded_python_processes=len(receipts),missing_direct_guard_receipts=missing,
        blocked_fallback_or_network_attempts=blocked,scientific_accuracy_validated=False)
    result['status']='passed' if result['n_checks']==result['n_passed'] and unchanged and not missing and not blocked else 'failed'
    (logs/'portable_replay_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ['snapshot','destination','logs']:parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--forbid-root',type=Path,action='append',required=True)
    parser.add_argument('--python',default=sys.executable);parser.add_argument('--timeout',type=int,default=180)
    args=parser.parse_args();result=verify(args.snapshot,args.destination,args.logs,args.forbid_root,args.python,args.timeout)
    print(json.dumps(result,indent=2))
    if result['status']!='passed':raise SystemExit(1)


if __name__=='__main__':main()
