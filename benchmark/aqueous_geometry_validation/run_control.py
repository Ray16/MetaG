"""One reserved CUDA worker; serial fixed endpoints with exact process ownership."""
from pathlib import Path
from datetime import datetime,timezone
import os,signal,subprocess,time,traceback,shutil
from control_common import HERE,authenticate,read,write,owner,sha
inputs,numeric=authenticate();manifest=read(HERE/'launch_manifest.json');protocol=read(HERE/'protocol.json')
assert sorted(os.sched_getaffinity(0))==manifest['cpus']
assert os.environ.get('CUDA_VISIBLE_DEVICES') is not None
write(HERE/'controller_owner.json',owner())
deadline=time.monotonic()+max(0,(datetime.fromisoformat(manifest['deadline_utc'])-datetime.now(timezone.utc)).total_seconds())
child=None

def terminate_child():
 global child
 if child is not None and child.poll() is None:
  os.killpg(child.pid,signal.SIGTERM)
  try:child.wait(timeout=15)
  except subprocess.TimeoutExpired:os.killpg(child.pid,signal.SIGKILL);child.wait()
 child=None

def stop(signum,frame):
 terminate_child();raise SystemExit(128+signum)
signal.signal(signal.SIGTERM,stop);signal.signal(signal.SIGINT,stop)

def execute(slot):
 global child
 remaining=deadline-time.monotonic()
 if remaining<=0:raise TimeoutError('Declared overall budget exhausted')
 argv=[manifest['worker_python'],str(HERE/'reference_point.py'),slot]
 with (HERE/(slot+'.console.log')).open('x') as stream:
  child=subprocess.Popen(argv,stdin=subprocess.DEVNULL,stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
  record=dict(host=os.uname().nodename,pid=child.pid,start_ticks=Path(f'/proc/{child.pid}/stat').read_text().rsplit(')',1)[1].split()[19],cpu_affinity=sorted(os.sched_getaffinity(child.pid)),argv=argv)
  write(HERE/(slot+'.predictor_owner.json'),record)
  assert record['cpu_affinity']==manifest['cpus']
  try:code=child.wait(timeout=remaining)
  except subprocess.TimeoutExpired:terminate_child();raise TimeoutError('Overall control budget exhausted during '+slot)
  finally:
   if child is not None and child.poll() is not None:child=None
 if code:raise RuntimeError(f'{slot} exited{code}; complete/partial outputs and logs preserved')
 result=read(HERE/'points'/(slot+'.json'));assert result['passed']
 assert result['source_sha256']==sha(HERE/'reference_point.py')
 assert result['protocol_sha256']==sha(HERE/'fixed_geometry_protocol.json')
 assert result['input_sha256']==sha(HERE/'geometry_inputs.json')
 if slot!='smoke':
  assert result['geometry']==inputs['inputs'][slot]['geometry']
  for phase in ['gas','water']:
   assert result[phase]['n_molecular_orbitals']==result['n_AO']
   assert result[phase]['converged'] and result[phase]['orbital_gradient_norm']<=numeric['stationary_gradient_maximum']
 for name,expected in read(HERE/'dependency_identity.json')['versions'].items():assert result['versions'][name]==expected,(name,result['versions'][name],expected)

status=dict(status='running',expected_slots=protocol['expected_slots'],completed=[],failures=[],smoke_passed=False)
write(HERE/'control_status.json',status)
try:
 execute('smoke');status['smoke_passed']=True;write(HERE/'control_status.json',status)
 consecutive=0
 for slot in protocol['expected_slots']:
  if max(os.getloadavg())+len(manifest['cpus'])>os.cpu_count():raise RuntimeError('CPU load exceeds declared capacity')
  if int(Path('/proc/meminfo').read_text().split('MemAvailable:')[1].split()[0])*1024<256*1024**3:raise RuntimeError('RAM precondition no longer met')
  if shutil.disk_usage('/tmp').free<50*1024**3:raise RuntimeError('Scratch precondition no longer met')
  claim=HERE/'claims'/slot;claim.parent.mkdir(exist_ok=True);claim.mkdir();write(claim/'owner.json',owner())
  attempt=HERE/'attempts'/slot;attempt.parent.mkdir(exist_ok=True);attempt.mkdir();write(attempt/'owner.json',owner())
  try:
   print('POINT START '+slot,flush=True);execute(slot)
   status['completed'].append(slot);consecutive=0;print('POINT COMPLETE '+slot,flush=True)
  except Exception as error:
   consecutive+=1;record=dict(slot=slot,error=str(error),traceback=traceback.format_exc())
   status['failures'].append(record);write(attempt/'failure.json',record)
   if consecutive>=3 or isinstance(error,TimeoutError):raise
  finally:
   terminate_child();shutil.rmtree(claim);write(HERE/'control_status.json',status)
 status['status']='complete' if len(status['completed'])==len(protocol['expected_slots']) else 'terminal_with_failures'
except BaseException as error:
 status.update(status='stopped',error=str(error),traceback=traceback.format_exc());raise
finally:
 terminate_child();write(HERE/'control_status.json',status)
