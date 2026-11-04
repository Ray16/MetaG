"""Bounded geometry-parent calculation ownership and immutable input checks."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,os,socket
HERE=Path(__file__).resolve().parent

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):return json.loads(Path(path).read_text())
def write(path,value):
 path=Path(path);temporary=path.with_name(path.name+f'.tmp-{os.getpid()}')
 temporary.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');temporary.replace(path)
def owner():
 return dict(host=socket.gethostname(),pid=os.getpid(),start_ticks=Path(f'/proc/{os.getpid()}/stat').read_text().rsplit(')',1)[1].split()[19],utc=datetime.now(timezone.utc).isoformat(),cpu_affinity=sorted(os.sched_getaffinity(0)),managed_cuda_visible_devices=os.environ.get('CUDA_VISIBLE_DEVICES'),native_thread_caps={key:os.environ.get(key) for key in ['OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS']})
def authenticate():
 manifest=read(HERE/'launch_manifest.json')
 for name,digest in manifest['files_sha256'].items():assert sha(HERE/name)==digest,name
 inputs=read(HERE/'geometry_inputs.json');protocol=read(HERE/'protocol.json')
 assert sha(HERE/'fixed_geometry_protocol.json')==protocol['numerical_protocol_sha256']
 assert sha(HERE/'reference_point.py')==protocol['parent_worker_source_sha256']
 for relative,record in inputs['source_records'].items():assert sha(HERE/relative)==record['sha256'],relative
 for slot,item in inputs['inputs'].items():
  candidate=read(HERE/item['source_candidate']);geometry=item['geometry']
  assert geometry['symbols']==candidate['symbols'] and geometry['positions_angstrom']==candidate['positions_A']
  assert geometry['charge']==candidate['charge'] and geometry['spin_multiplicity']==1
  assert sha(HERE/item['source_candidate'])==item['original_candidate_sha256']
  value={k:v for k,v in geometry.items() if k!='sha256'}
  assert hashlib.sha256(json.dumps(value,sort_keys=True,allow_nan=False).encode()).hexdigest()==geometry['sha256']
 deps=read(HERE/'dependency_identity.json')
 for path,digest in deps['installed_file_sha256'].items():assert sha(path)==digest,path
 return inputs,read(HERE/'fixed_geometry_protocol.json')
