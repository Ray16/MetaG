"""Retain lease until the exact numerical child has terminated and been joined."""
from pathlib import Path
import datetime
import json
import os
import signal
import subprocess

HERE = Path(__file__).resolve().parent
m = json.loads((HERE/'launch_manifest.json').read_text())
assert sorted(os.sched_getaffinity(0)) == m['cpus']
assert os.environ.get('CUDA_VISIBLE_DEVICES') is not None
def record(pid):
    return dict(host=os.uname().nodename, pid=pid,
        start_ticks=Path(f'/proc/{pid}/stat').read_text().rsplit(')',1)[1].split()[19])
(HERE/'controller_owner.json').write_text(json.dumps(record(os.getpid()),indent=2)+'\n')
child = None
def terminate():
    if child is not None and child.poll() is None:
        os.killpg(child.pid,signal.SIGTERM)
        try: child.wait(timeout=15)
        except subprocess.TimeoutExpired:
            os.killpg(child.pid,signal.SIGKILL)
            child.wait()
def stop(signum,frame):
    terminate()
    raise SystemExit(128+signum)
signal.signal(signal.SIGTERM,stop)
signal.signal(signal.SIGINT,stop)
remaining = (datetime.datetime.fromisoformat(m['deadline_utc'])-datetime.datetime.now(datetime.timezone.utc)).total_seconds()
assert 0 < remaining <= 1200
with (HERE/'predictor.console.log').open('x') as log:
    child = subprocess.Popen([m['python'],str(HERE/'predict_component.py')],stdin=subprocess.DEVNULL,
        stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
    (HERE/'predictor_process_owner.json').write_text(json.dumps(record(child.pid),indent=2)+'\n')
    try: code = child.wait(timeout=remaining)
    except subprocess.TimeoutExpired:
        terminate()
        (HERE/'budget_exhausted.json').write_text(json.dumps(dict(child_joined=True))+'\n')
        raise
    finally: terminate()
raise SystemExit(code)
