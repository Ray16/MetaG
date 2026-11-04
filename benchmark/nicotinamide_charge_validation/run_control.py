"""Bounded CPU-only controller; one CUDA child and one exact claim at a time."""
from pathlib import Path
import datetime
import os
import shutil
import signal
import subprocess
import time
import traceback
from control_common import HERE,authenticate,read,owner,sha,write

inputs,protocol,numeric=authenticate();manifest=read(HERE/'launch_manifest.json')
assert sorted(os.sched_getaffinity(0))==manifest['cpus']
assert os.environ.get('CUDA_VISIBLE_DEVICES') is not None
write(HERE/'controller_owner.json',owner())
remaining=(datetime.datetime.fromisoformat(manifest['deadline_utc'])-datetime.datetime.now(datetime.timezone.utc)).total_seconds()
deadline=time.monotonic()+min(1800,remaining)
child=None
def terminate_child():
    global child
    if child is not None and child.poll() is None:
        os.killpg(child.pid,signal.SIGTERM)
        try:child.wait(timeout=15)
        except subprocess.TimeoutExpired:
            os.killpg(child.pid,signal.SIGKILL);child.wait()
    child=None
def stop(signum,frame):
    terminate_child();raise SystemExit(128+signum)
signal.signal(signal.SIGTERM,stop);signal.signal(signal.SIGINT,stop)
def run(slot):
    global child
    remaining=deadline-time.monotonic()
    if remaining<=0:raise TimeoutError('Overall thirty-minute budget exhausted')
    argv=[manifest['python'],str(HERE/'aqueous_point.py'),slot]
    with (HERE/(slot+'.console.log')).open('x') as stream:
        child=subprocess.Popen(argv,stdin=subprocess.DEVNULL,stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
        write(HERE/(slot+'.predictor_owner.json'),dict(pid=child.pid,host=os.uname().nodename,
            start_ticks=Path(f'/proc/{child.pid}/stat').read_text().rsplit(')',1)[1].split()[19],
            cpu_affinity=sorted(os.sched_getaffinity(child.pid)),argv=argv))
        try:code=child.wait(timeout=remaining)
        except subprocess.TimeoutExpired:
            terminate_child();raise TimeoutError(f'Overall budget exhausted during{slot}')
        finally:
            if child is not None and child.poll() is not None:child=None
    if code:raise RuntimeError(f'{slot} returned{code}; retain console log and partials')
    value=read(HERE/'points'/(slot+'.json'));assert value['passed']
    assert value['source_sha256']==sha(HERE/'aqueous_point.py')
    if slot!='smoke':assert value['geometry']==inputs[slot]['geometry']
status=dict(expected_slots=['oxidized','reduced'],completed=[],failures=[],status='running')
write(HERE/'control_status.json',status)
try:
    run('smoke')
    assert read(HERE/'points/smoke.json')['passed']
    failures=0
    for slot in ['oxidized','reduced']:
        assert int(Path('/proc/meminfo').read_text().split('MemAvailable:')[1].split()[0])*1024>=256*1024**3
        assert shutil.disk_usage('/tmp').free>=50*1024**3
        claims=HERE/'claims';claims.mkdir(exist_ok=True);claim=claims/slot;claim.mkdir()
        write(claim/'owner.json',owner())
        (HERE/'attempts').mkdir(exist_ok=True);attempt=HERE/'attempts'/slot;attempt.mkdir()
        write(attempt/'owner.json',owner())
        try:
            print('POINT START '+slot,flush=True);run(slot)
            status['completed'].append(slot);failures=0
            print('POINT COMPLETE '+slot,flush=True)
        except Exception as error:
            failures+=1;record=dict(slot=slot,error=str(error),traceback=traceback.format_exc())
            status['failures'].append(record);write(attempt/'failure.json',record)
            if failures>=3 or isinstance(error,TimeoutError):raise
        finally:
            # The child is joined before ownership is released, including termination.
            terminate_child();shutil.rmtree(claim);write(HERE/'control_status.json',status)
    status['status']='complete' if len(status['completed'])==2 else 'terminal_with_failures'
except BaseException as error:
    status.update(status='stopped',error=str(error),traceback=traceback.format_exc());raise
finally:
    terminate_child();write(HERE/'control_status.json',status)
