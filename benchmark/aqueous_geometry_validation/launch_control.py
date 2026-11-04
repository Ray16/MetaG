"""Preflight, pin producer/input ownership, then reserve exactly one CUDA worker."""
from pathlib import Path
from datetime import datetime,timezone,timedelta
import ast,json,os,shutil,subprocess,sys
from control_common import HERE,sha,read,write
assert not (HERE/'launch_manifest.json').exists(),'Inspect an existing run before any restart'
protocol=read(HERE/'protocol.json');resources=protocol['resource_ceiling'];cpus=resources['cpus']
assert os.uname().nodename==resources['host']
worker_python='/homes/rzhu/miniforge3/envs/redox/bin/python'
controller_python='/homes/rzhu/miniforge3/envs/uma/bin/python'
helper='/nfs/lambda_stor_01/homes/rzhu/.codex-personal/skills/lambda-fleet/scripts/preflight_node.py'
gate='/nfs/lambda_stor_01/homes/rzhu/bin/gpu_reserve'
environment=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',PYTHONNOUSERSITE='1',LD_LIBRARY_PATH='/homes/rzhu/miniforge3/envs/redox/lib')
for name in ['OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS']:environment[name]='1'
argv=[controller_python,helper,'local','--python',worker_python,'--cpu-set',','.join(map(str,cpus)),
 '--min-available-ram-mib','262144','--module','numpy','--module','pyscf','--module','ase','--bin',gate,
 '--bin','/usr/bin/nvidia-smi','--path',str(HERE/'geometry_inputs.json'),'--library-path',environment['LD_LIBRARY_PATH']]
r=subprocess.run(argv,env=environment,capture_output=True,text=True,timeout=60)
(HERE/'dependency_preflight.json').write_text(r.stdout)
if r.returncode:print(r.stdout);print(r.stderr);r.check_returncode()
assert read(HERE/'dependency_preflight.json')['ok']
deps=read(HERE/'dependency_identity.json');assert deps['passed']
assert deps['geometry_inputs_sha256']==sha(HERE/'geometry_inputs.json')
assert deps['numeric_settings_sha256']==sha(HERE/'fixed_geometry_protocol.json')
assert max(os.getloadavg())+len(cpus)<=os.cpu_count()
assert shutil.disk_usage('/tmp').free>=50*1024**3
for path in HERE.glob('*.py'):ast.parse(path.read_text(),filename=str(path))
gpu_index=resources['gpu_index']
query=subprocess.run(['/usr/bin/nvidia-smi',f'--id={gpu_index}','--query-gpu=index,uuid,name,memory.total,memory.used,utilization.gpu','--format=csv,noheader,nounits'],capture_output=True,text=True,check=True).stdout.strip().split(', ')
assert int(query[0])==gpu_index and int(query[3])>=32000 and int(query[4])<=500 and int(query[5])<=5
uuid=query[1]
resident=subprocess.run(['/usr/bin/nvidia-smi','--query-compute-apps=gpu_uuid,pid,process_name','--format=csv,noheader'],capture_output=True,text=True,check=True).stdout
assert uuid not in resident,'Resource contention; no geometry attempt consumed'
write(HERE/'launch_resources.json',dict(host=os.uname().nodename,cpus=cpus,load=list(os.getloadavg()),available_ram_mib=int(Path('/proc/meminfo').read_text().split('MemAvailable:')[1].split()[0])/1024,scratch_free_bytes=shutil.disk_usage('/tmp').free,gpu=query,resident_compute_snapshot=resident))
files=['protocol.json','fixed_geometry_protocol.json','geometry_inputs.json','control_common.py','reference_point.py','run_control.py','launch_control.py','prepare_inputs.py','preflight_basis.py','cpu_smoke.py','cpu_water_smoke.json','dependency_identity.json','dependency_preflight.json','resources.json','launch_resources.json']
files += [str(path.relative_to(HERE)) for path in sorted((HERE/'parents').rglob('*')) if path.is_file()]
now=datetime.now(timezone.utc)
argv=[gate,'run',str(gpu_index),'--label','metag_aqueous_geometry_parent','--','/usr/bin/taskset','-c',','.join(map(str,cpus)),controller_python,str(HERE/'run_control.py')]
manifest=dict(utc=now.isoformat(),deadline_utc=(now+timedelta(seconds=resources['overall_seconds'])).isoformat(),host=os.uname().nodename,gpu_index=gpu_index,gpu_uuid=uuid,cpus=cpus,worker_python=worker_python,controller_python=controller_python,argv=argv,expected_slots=protocol['expected_slots'],files_sha256={name:sha(HERE/name) for name in files},scope=protocol['scope'])
write(HERE/'launch_manifest.json',manifest)
with (HERE/'launcher.log').open('x') as log:
 process=subprocess.Popen(argv,env=environment,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
identity=dict(host=os.uname().nodename,pid=process.pid,start_ticks=Path(f'/proc/{process.pid}/stat').read_text().rsplit(')',1)[1].split()[19])
write(HERE/'launcher_owner.json',identity);print(json.dumps(identity))
