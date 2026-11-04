"""Owned records for the separately bounded phosphate charge-state control."""
from pathlib import Path
from collections import Counter
import datetime
import hashlib
import json
import os

HERE=Path(__file__).resolve().parent
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):return json.loads(Path(path).read_text())
def write(path,value):
    path=Path(path);temporary=path.with_name(path.name+f'.tmp-{os.getpid()}')
    temporary.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');temporary.replace(path)
def owner():
    return dict(host=os.uname().nodename,pid=os.getpid(),
        start_ticks=Path(f'/proc/{os.getpid()}/stat').read_text().rsplit(')',1)[1].split()[19],
        utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        cpu_affinity=sorted(os.sched_getaffinity(0)),managed_cuda_visible_devices=os.environ.get('CUDA_VISIBLE_DEVICES'),
        native_thread_caps={k:os.environ.get(k) for k in ['OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS']})
def authenticate():
    manifest=read(HERE/'launch_manifest.json')
    for name,digest in manifest['files_sha256'].items():assert sha(HERE/name)==digest,name
    protocol=read(HERE/'protocol.json');inputs=read(HERE/'geometry_inputs.json')['states']
    assert sha(HERE/'numerical_protocol.json')==protocol['parent_protocol_sha256']
    assert sha(HERE/'reference_geometry_inputs.json')==protocol['parent_input_sha256']
    parent=read(HERE/'reference_geometry_inputs.json')['inputs']
    for family,item in inputs.items():
        source=parent[item['parent_slot']]['geometry'];geometry=item['geometry']
        assert source['sha256']==item['parent_geometry_sha256']
        keep=[i for i in range(len(source['symbols'])) if i not in item['removed_H_indices']]
        assert len(item['removed_H_indices'])==2
        assert all(source['symbols'][i]=='H' for i in item['removed_H_indices'])
        assert geometry['symbols']==[source['symbols'][i] for i in keep]
        assert geometry['positions_angstrom']==[source['positions_angstrom'][i] for i in keep]
        assert geometry['charge']==source['charge']-2
        value={k:v for k,v in geometry.items() if k!='sha256'}
        assert hashlib.sha256(json.dumps(value,sort_keys=True,allow_nan=False).encode()).hexdigest()==geometry['sha256']
        point=HERE/'parent_points'/(item['parent_slot']+'.json')
        assert sha(point)==protocol['parent_reference_points'][item['parent_slot']]
        assert read(point)['source_sha256']==sha(HERE/'parent_reference_point.py')
    return inputs,protocol,read(HERE/'numerical_protocol.json')
