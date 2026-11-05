"""Bounded two-worker CPU queue; retain every gas/ALPB native log and refusal."""
from pathlib import Path
import argparse
import datetime
import hashlib
import json
import math
import os
import re
import shutil
import signal
import subprocess
import tempfile
import time
import traceback

HERE=Path(__file__).resolve().parent
sha=lambda path:hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(name):return json.loads((HERE/name).read_text())
def write(path,value):
    path=Path(path);temporary=path.with_name(path.name+f'.tmp-{os.getpid()}')
    temporary.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');temporary.replace(path)
def owner():
    return dict(host=os.uname().nodename,pid=os.getpid(),start_ticks=Path(f'/proc/{os.getpid()}/stat').read_text().rsplit(')',1)[1].split()[19],
        cpu_affinity=sorted(os.sched_getaffinity(0)),utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
parser=argparse.ArgumentParser();parser.add_argument('worker',type=int,choices=[0,1]);args=parser.parse_args()
protocol=read('protocol.json');launch=read('launch_manifest.json');tasks=read('tasks.json')
for name,digest in read('freeze.json')['sha256'].items():assert sha(HERE/name)==digest,name
assert sha(protocol['xtb'])==protocol['xtb_sha256']
assert sorted(os.sched_getaffinity(0))==[protocol['worker_CPUs'][args.worker]]
assert all(os.environ.get(name)=='1' for name in ['OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS'])
assert max(os.getloadavg())+3<=os.cpu_count()
assert int(Path('/proc/meminfo').read_text().split('MemAvailable:')[1].split()[0])*1024>=8*1024**3
worker_dir=HERE/'workers'/str(args.worker);worker_dir.mkdir(exist_ok=False)
write(worker_dir/'owner.json',owner())
deadline=datetime.datetime.fromisoformat(launch['deadline_utc']);child=None
def remaining():return (deadline-datetime.datetime.now(datetime.timezone.utc)).total_seconds()
def join_child():
    global child
    if child is not None:
        if child.poll() is None:
            os.killpg(child.pid,signal.SIGTERM)
            try:child.wait(timeout=10)
            except subprocess.TimeoutExpired:os.killpg(child.pid,signal.SIGKILL);child.wait()
        child=None
def stop(signum,frame):
    join_child();raise SystemExit(128+signum)
signal.signal(signal.SIGTERM,stop);signal.signal(signal.SIGINT,stop)
number=r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?'
status=dict(state='running',completed=[],failures=[],active=None)
failures=0
try:
    for task in tasks:
        key=task['key']
        if int(hashlib.sha256(key.encode()).hexdigest(),16)%2!=args.worker:continue
        if (HERE/'results'/(key+'.json')).exists():continue
        if remaining()<=0:raise TimeoutError('Original study deadline reached')
        assert shutil.disk_usage('/tmp').free>5*1024**3
        claim=HERE/'claims'/key
        try:claim.mkdir()
        except FileExistsError:continue
        write(claim/'owner.json',owner());status['active']=key
        attempt=HERE/'attempts'/key
        try:
            attempt.mkdir(exist_ok=False);write(attempt/'owner.json',owner())
            print('START '+key,flush=True);started=time.monotonic()
            xyz=str(len(task['symbols']))+'\n\n'+''.join(s+' '+' '.join(f'{float(v):.6f}' for v in point)+'\n'
                for s,point in zip(task['symbols'],task['positions_angstrom']))
            native=[]
            with tempfile.TemporaryDirectory(prefix='metag-solvation-rotation-') as temp:
                work=Path(temp);(work/'input.xyz').write_text(xyz)
                for phase in ['gas','water']:
                    argv=[protocol['xtb'],str(work/'input.xyz'),'--gfn','2','--chrg',str(task['charge']),'--uhf','0','--sp']
                    if task['accuracy']!='default':argv+=['--acc',task['accuracy']]
                    if phase=='water':
                        argv+=['--alpb','water']
                        if task['grid']=='extreme':
                            (work/'solvent.in').write_text('$gbsa\n grid=extreme\n$end\n')
                            argv+=['--input',str(work/'solvent.in')]
                    restart=work/'xtbrestart'
                    preexisting_restart_sha256=sha(restart) if restart.exists() else None
                    child=subprocess.Popen(argv,cwd=work,stdin=subprocess.DEVNULL,stdout=subprocess.PIPE,stderr=subprocess.PIPE,
                        text=True,start_new_session=True)
                    identity=dict(host=os.uname().nodename,pid=child.pid,
                        start_ticks=Path(f'/proc/{child.pid}/stat').read_text().rsplit(')',1)[1].split()[19],
                        cpu_affinity=sorted(os.sched_getaffinity(child.pid)))
                    write(attempt/(phase+'.owner.json'),identity)
                    try:stdout,stderr=child.communicate(timeout=max(.1,min(protocol['per_native_timeout_seconds'],remaining())))
                    except subprocess.TimeoutExpired:
                        timed_out_child=child;join_child()
                        stdout,stderr=timed_out_child.communicate(timeout=5)
                        write(attempt/(phase+'.timeout.json'),dict(phase=phase,argv=argv,owner=identity,input_xyz=xyz,
                            stdout=stdout,stderr=stderr,returncode=timed_out_child.returncode,timed_out=True))
                        raise TimeoutError('Native point/deadline budget exhausted')
                    code=child.returncode;child=None
                    match=re.search(r'TOTAL ENERGY\s+('+number+r')\s+Eh',stdout)
                    value=float(match[1].replace('D','E')) if match else None
                    record=dict(phase=phase,argv=argv,owner=identity,input_xyz=xyz,returncode=code,energy_Eh=value,
                        preexisting_restart_sha256=preexisting_restart_sha256,stdout=stdout,stderr=stderr)
                    write(attempt/(phase+'.json'),record);native.append(record)
                    assert code==0 and value is not None and math.isfinite(value)
                    assert 'convergence criteria satisfied' in stdout
                    if phase=='water':assert re.search(r'Grid points\s+'+('5810' if task['grid']=='extreme' else '230')+r'\s+per atom',stdout)
            delta=(native[1]['energy_Eh']-native[0]['energy_Eh'])*protocol['Hartree_to_kj_per_mol']
            result=dict(key=key,task=task,gas_Eh=native[0]['energy_Eh'],water_Eh=native[1]['energy_Eh'],
                solvation_kj=delta,elapsed_seconds=time.monotonic()-started,
                native_record_sha256={phase:sha(attempt/(phase+'.json')) for phase in ['gas','water']},
                protocol_sha256=sha(HERE/'protocol.json'),worker_sha256=sha(__file__),passed=True)
            write(HERE/'results'/(key+'.json'),result);status['completed'].append(key);failures=0
            print('COMPLETE '+key+f' dGsolv={delta:.9g}',flush=True)
        except Exception as error:
            failure=dict(key=key,error=str(error),traceback=traceback.format_exc())
            write(attempt/'failure.json',failure);status['failures'].append(failure);failures+=1
            if failures>=protocol['consecutive_failure_limit']:raise
            if remaining()<=0:raise
        finally:
            join_child();shutil.rmtree(claim);status['active']=None;write(worker_dir/'status.json',status)
    status['state']='complete' if not status['failures'] else 'terminal_with_failures'
except BaseException as error:
    status.update(state='stopped',error=str(error),traceback=traceback.format_exc());raise
finally:
    join_child();write(worker_dir/'status.json',status)
