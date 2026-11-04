"""Preflight, freeze ownership, and atomically reserve the one declared GPU."""
from pathlib import Path
import ast
import datetime
import hashlib
import json
import os
import shutil
import subprocess
import sys

HERE=Path(__file__).resolve().parent
assert not (HERE/'launch_manifest.json').exists(), 'Existing run: inspect before resuming'
cpu_set=[2,3,4,5]
uma_python='/nfs/lambda_stor_01/homes/rzhu/.cache/metag-environments/uma-nitro-20261023-645af7d/venv/bin/python'
dft_python='/homes/rzhu/miniforge3/envs/redox/bin/python'
redox_libraries='/homes/rzhu/miniforge3/envs/redox/lib'
helper='/nfs/lambda_stor_01/homes/rzhu/.codex-personal/skills/lambda-fleet/scripts/preflight_node.py'
gate='/nfs/lambda_stor_01/homes/rzhu/bin/gpu_reserve'
environment=dict(os.environ)
for name in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    environment[name]='1'
environment.update(PYTHONDONTWRITEBYTECODE='1',HF_HUB_OFFLINE='1',
    PYTHONPATH=os.pathsep.join([str(HERE/'uma_source/src'),str(HERE/'fairchem_sdk/src')]))
for name,python,modules in [('uma',uma_python,['torch','ase','fairchem.core','metag.energetics.uma']),
                            ('dft',dft_python,['numpy','pyscf','ase'])]:
    argv=[sys.executable,helper,'local','--python',python,'--cpu-set',','.join(map(str,cpu_set)),
          '--min-available-ram-mib','262144','--bin','/usr/bin/nvidia-smi','--bin',gate,
          '--path',str(HERE/'geometry_inputs.json'),'--library-path',redox_libraries]
    for module in modules:argv.extend(['--module',module])
    result=subprocess.run(argv,capture_output=True,text=True,env=environment,timeout=60)
    (HERE/(name+'_final_preflight.json')).write_text(result.stdout)
    if result.returncode:
        print(result.stdout);print(result.stderr);result.check_returncode()
for path in HERE.glob('*.py'):ast.parse(path.read_text(),filename=str(path))
assert max(os.getloadavg())+len(cpu_set)<=os.cpu_count()
assert shutil.disk_usage('/tmp').free>=50*1024**3
processes=subprocess.run(['/usr/bin/nvidia-smi','--query-compute-apps=gpu_uuid,pid,process_name',
                          '--format=csv,noheader'],capture_output=True,text=True,check=True).stdout
uuid='GPU-fab99e67-e6e2-56bb-e07d-ea240da3e84a'
assert uuid not in processes,'Chosen GPU acquired a resident process; no input attempt consumed'
argv=[gate,'run','7','--label','metag_nicotinamide_reference_control','--',
      '/usr/bin/taskset','-c',','.join(map(str,cpu_set)),sys.executable,str(HERE/'run_control.py')]
files=['geometry_inputs.json','fixed_geometry_protocol.json','dependency_sources.json','cpu_smoke.py',
       'cpu_water_smoke.json','control_common.py','replay_uma.py','reference_point.py','run_control.py',
       'launch_control.py','uma_final_preflight.json','dft_final_preflight.json']
if (HERE/'uma_replay.json').exists():
    files.extend(['uma_replay.json','infrastructure_attempt1/launch_manifest.json','infrastructure_revision.json'])
if (HERE/'infrastructure_attempt2/launch_manifest.json').exists():
    files.append('infrastructure_attempt2/launch_manifest.json')
original_start=datetime.datetime.now(datetime.timezone.utc)
if (HERE/'infrastructure_attempt1/launch_manifest.json').exists():
    original_start=datetime.datetime.fromisoformat(json.loads((HERE/'infrastructure_attempt1/launch_manifest.json').read_text())['utc'])
deadline=original_start+datetime.timedelta(seconds=5400)
assert datetime.datetime.now(datetime.timezone.utc)<deadline
manifest=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),host=os.uname().nodename,
    gpu_index=7,gpu_uuid=uuid,cpus=cpu_set,argv=argv,uma_python=uma_python,dft_python=dft_python,
    redox_library_path=redox_libraries,overall_runtime_budget_seconds=5400,
    deadline_utc=deadline.isoformat(),
    expected_geometry_slots=sorted(json.loads((HERE/'geometry_inputs.json').read_text())['inputs']),
    owned_file_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in files},
    dependency_preflight='uma_final_preflight.json; dft_final_preflight.json',
    scope='One bounded fixed-geometry reference control; no production change or empirical fit.')
(HERE/'launch_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
with (HERE/'launcher.log').open('x') as log:
    process=subprocess.Popen(argv,env=environment,stdin=subprocess.DEVNULL,stdout=log,
                             stderr=subprocess.STDOUT,start_new_session=True)
identity=dict(host=os.uname().nodename,pid=process.pid,
    start_ticks=Path(f'/proc/{process.pid}/stat').read_text().rsplit(')',1)[1].split()[19])
(HERE/'launcher_owner.json').write_text(json.dumps(identity,indent=2)+'\n')
print(json.dumps(identity))
