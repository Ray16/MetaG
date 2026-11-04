"""CPU-only controller; one CUDA child at a time under a single gate lease."""
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time
import traceback
from control_common import HERE, authenticate, owner, sha, write

inputs, protocol = authenticate()
manifest = json.loads((HERE/'launch_manifest.json').read_text())
assert sorted(os.sched_getaffinity(0)) == manifest['cpus']
assert os.environ.get('CUDA_VISIBLE_DEVICES') is not None, 'reservation gate required'
assert max(os.getloadavg())+len(manifest['cpus']) <= os.cpu_count()
write(HERE/'controller_owner.json', owner())
deadline = time.monotonic()+protocol['runtime_budget_seconds']
child = None

def stop(signum, frame):
    global child
    if child is not None and child.poll() is None:
        os.killpg(child.pid, signal.SIGTERM)
        try: child.wait(timeout=15)
        except subprocess.TimeoutExpired:
            os.killpg(child.pid, signal.SIGKILL)
            child.wait()
    raise SystemExit(128+signum)

signal.signal(signal.SIGTERM, stop)
signal.signal(signal.SIGINT, stop)

def run(argv, name, env):
    global child
    remaining = deadline-time.monotonic()
    if remaining <= 0: raise TimeoutError('Overall control budget exhausted')
    with (HERE/(name+'.console.log')).open('x') as stream:
        child = subprocess.Popen(argv,env=env,stdin=subprocess.DEVNULL,stdout=stream,
                                 stderr=subprocess.STDOUT,start_new_session=True)
        write(HERE/(name+'.predictor_owner.json'), dict(pid=child.pid,
              start_ticks=Path(f'/proc/{child.pid}/stat').read_text().rsplit(')',1)[1].split()[19],
              argv=argv,host=os.uname().nodename,cpu_affinity=sorted(os.sched_getaffinity(child.pid))))
        try: code = child.wait(timeout=remaining)
        except subprocess.TimeoutExpired:
            os.killpg(child.pid,signal.SIGTERM)
            try:child.wait(timeout=15)
            except subprocess.TimeoutExpired:
                os.killpg(child.pid,signal.SIGKILL);child.wait()
            raise TimeoutError(f'Overall budget exhausted in {name}')
        finally:
            if child.poll() is not None:child=None
    if code:raise RuntimeError(f'{name} returned {code}; inspect durable console log')

environment = dict(os.environ)
environment.update(PYTHONDONTWRITEBYTECODE='1',HF_HUB_OFFLINE='1',
    PYTHONPATH=os.pathsep.join([str(HERE/'uma_source/src'),str(HERE/'fairchem_sdk/src')]))
summary = dict(expected_slots=sorted(inputs['inputs']), completed=[], failures=[], status='running')
write(HERE/'control_status.json',summary)
try:
    run([manifest['uma_python'],str(HERE/'replay_uma.py')],'uma_replay',environment)
    assert json.loads((HERE/'uma_replay.json').read_text())['passed']
    environment.pop('PYTHONPATH')
    environment['LD_LIBRARY_PATH']=manifest['redox_library_path']
    run([manifest['dft_python'],str(HERE/'reference_point.py'),'smoke'],'gpu_smoke',environment)
    assert json.loads((HERE/'points/smoke.json').read_text())['passed']
    failures = 0
    order = sorted(inputs['inputs'],key=lambda key:(key.startswith('full'),key))
    for slot in order:
        if (HERE/'attempts'/slot).exists():
            raise RuntimeError(f'Attempt already exists for {slot}; no automatic rerun')
        assert int(Path('/proc/meminfo').read_text().split('MemAvailable:')[1].split()[0])*1024 >= 256*1024**3
        assert shutil.disk_usage('/tmp').free >= 50*1024**3
        claims = HERE/'claims';claims.mkdir(exist_ok=True)
        claim = claims/slot
        claim.mkdir()
        write(claim/'owner.json',owner())
        (HERE/'attempts').mkdir(exist_ok=True)
        (HERE/'attempts'/slot).mkdir()
        write(HERE/'attempts'/slot/'owner.json',owner())
        try:
            print(f'POINT START {slot}',flush=True)
            run([manifest['dft_python'],str(HERE/'reference_point.py'),slot],slot,environment)
            value=json.loads((HERE/'points'/(slot+'.json')).read_text())
            assert value['passed'] and value['geometry']['sha256']==inputs['inputs'][slot]['geometry']['sha256']
            assert value['source_sha256']==sha(HERE/'reference_point.py')
            summary['completed'].append(slot);failures=0
            print(f'POINT COMPLETE {slot}',flush=True)
        except Exception as error:
            failures+=1
            record=dict(slot=slot,error=str(error),traceback=traceback.format_exc())
            summary['failures'].append(record)
            write(HERE/'attempts'/slot/'failure.json',record)
            if failures>=protocol['consecutive_actual_failure_limit']:raise
            if isinstance(error,TimeoutError):raise
        finally:
            # run() has joined its exact child before claim release.
            shutil.rmtree(claim)
            write(HERE/'control_status.json',summary)
    summary['status']='complete' if len(summary['completed'])==len(inputs['inputs']) else 'terminal_with_failures'
except BaseException as error:
    summary.update(status='stopped',controller_error=str(error),traceback=traceback.format_exc())
    raise
finally:
    write(HERE/'control_status.json',summary)
