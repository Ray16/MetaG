"""Read-only preflight, frozen manifest and atomic single-GPU submission."""
from pathlib import Path
import ast
import datetime
import json
import os
import shutil
import subprocess
import sys
from control_common import HERE,read,sha,write

assert not (HERE/'launch_manifest.json').exists(),'Inspect an existing run before any restart'
cpus=[2,3,4,5];python='/homes/rzhu/miniforge3/envs/redox/bin/python'
gate='/nfs/lambda_stor_01/homes/rzhu/bin/gpu_reserve'
helper='/nfs/lambda_stor_01/homes/rzhu/.codex-personal/skills/lambda-fleet/scripts/preflight_node.py'
environment=dict(os.environ,LD_LIBRARY_PATH='/homes/rzhu/miniforge3/envs/redox/lib',PYTHONDONTWRITEBYTECODE='1',
    OMP_NUM_THREADS='1',OMP_THREAD_LIMIT='1',MKL_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1')
argv=[sys.executable,helper,'local','--python',python,'--cpu-set',','.join(map(str,cpus)),
    '--min-available-ram-mib','262144','--module','numpy','--module','pyscf','--module','ase',
    '--bin',gate,'--bin','/usr/bin/nvidia-smi','--path',str(HERE/'geometry_inputs.json'),
    '--library-path',environment['LD_LIBRARY_PATH']]
r=subprocess.run(argv,env=environment,capture_output=True,text=True,timeout=60)
(HERE/'dependency_preflight.json').write_text(r.stdout)
if r.returncode:print(r.stdout);print(r.stderr);r.check_returncode()
assert max(os.getloadavg())+len(cpus)<=os.cpu_count()
assert shutil.disk_usage('/tmp').free>=50*1024**3
for path in HERE.glob('*.py'):ast.parse(path.read_text(),filename=str(path))
uuid='GPU-fab99e67-e6e2-56bb-e07d-ea240da3e84a'
resident=subprocess.run(['/usr/bin/nvidia-smi','--query-compute-apps=gpu_uuid,pid','--format=csv,noheader'],capture_output=True,text=True,check=True).stdout
assert uuid not in resident,'Contention; no geometry attempt consumed'
now=datetime.datetime.now(datetime.timezone.utc)
argv=[gate,'run','7','--label','metag_nicotinamide_phosphate_charge','--','/usr/bin/taskset','-c',','.join(map(str,cpus)),sys.executable,str(HERE/'run_control.py')]
files=['protocol.json','geometry_inputs.json','basis_preflight.json','reference_geometry_inputs.json',
    'numerical_protocol.json','parent_reference_point.py','cpu_water_smoke.json','cpu_smoke.py',
    'control_common.py','aqueous_point.py','run_control.py','launch_control.py','dependency_preflight.json']
files.extend(str(p.relative_to(HERE)) for p in (HERE/'parent_points').glob('*.json'))
manifest=dict(utc=now.isoformat(),deadline_utc=(now+datetime.timedelta(seconds=1800)).isoformat(),
    host=os.uname().nodename,gpu_index=7,gpu_uuid=uuid,cpus=cpus,python=python,argv=argv,
    expected_slots=['oxidized','reduced'],attempts_per_slot=1,consecutive_failure_limit=3,
    files_sha256={p:sha(HERE/p) for p in files},scope='Two fixed-geometry aqueous electronic points only; no biochemical free-energy or production-policy change.')
write(HERE/'launch_manifest.json',manifest)
with (HERE/'launcher.log').open('x') as log:
    process=subprocess.Popen(argv,env=environment,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
record=dict(host=os.uname().nodename,pid=process.pid,start_ticks=Path(f'/proc/{process.pid}/stat').read_text().rsplit(')',1)[1].split()[19])
write(HERE/'launcher_owner.json',record);print(json.dumps(record))
