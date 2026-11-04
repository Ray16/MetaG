"""Durable records shared by this bounded, fixed-geometry control."""
from pathlib import Path
import datetime
import hashlib
import json
import os
import socket

HERE = Path(__file__).resolve().parent

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def write(path, value):
    path = Path(path)
    temporary = path.with_name(path.name + f'.tmp-{os.getpid()}')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)

def owner():
    return dict(host=socket.gethostname(), pid=os.getpid(),
                start_ticks=Path(f'/proc/{os.getpid()}/stat').read_text().rsplit(')', 1)[1].split()[19],
                utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                cpu_affinity=sorted(os.sched_getaffinity(0)),
                managed_cuda_visible_devices=os.environ.get('CUDA_VISIBLE_DEVICES'),
                native_thread_caps={key: os.environ.get(key) for key in
                    ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS',
                     'OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS')})

def authenticate():
    manifest = json.loads((HERE/'launch_manifest.json').read_text())
    for relative, expected in manifest['owned_file_sha256'].items():
        assert sha(HERE/relative) == expected, relative
    inputs = json.loads((HERE/'geometry_inputs.json').read_text())
    for source in inputs['source_records'].values():
        assert sha(source['path']) == source['sha256'], source['path']
    for key, item in inputs['inputs'].items():
        source = inputs['source_records'][item['family']]
        raw = json.loads(Path(source['path']).read_text())
        minimum = raw['meta']['thermal']['ensemble']['validated_minima'][item['minimum_index']]
        assert minimum['geometry'] == item['geometry'], key
        assert minimum['reference_smiles'] == item['reference_smiles'], key
        for component, value in item['original_components_kj'].items():
            assert minimum[component] == value, (key, component)
    deps = json.loads((HERE/'dependency_sources.json').read_text())
    for relative, expected in deps['files_sha256'].items():
        assert sha(HERE/relative) == expected, relative
    return inputs, json.loads((HERE/'fixed_geometry_protocol.json').read_text())
