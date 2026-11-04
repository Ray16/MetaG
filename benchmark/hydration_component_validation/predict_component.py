"""One pinned, reserved GPU; retain every source-geometry single point."""
from pathlib import Path
import hashlib
import json
import os
import shutil
import signal
import time
import traceback

HERE = Path(__file__).resolve().parent
read = lambda p: json.loads(Path(p).read_text())
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p, value):
    p = Path(p)
    temp = p.with_suffix(p.suffix + '.tmp')
    temp.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temp.replace(p)
def owner():
    return dict(host=os.uname().nodename, pid=os.getpid(),
        start_ticks=Path(f'/proc/{os.getpid()}/stat').read_text().rsplit(')', 1)[1].split()[19],
        affinity=sorted(os.sched_getaffinity(0)), CUDA_VISIBLE_DEVICES=os.environ.get('CUDA_VISIBLE_DEVICES'))
def stop(signum, frame):
    raise SystemExit(128 + signum)
signal.signal(signal.SIGTERM, stop)
signal.signal(signal.SIGINT, stop)
protocol = read(HERE / 'protocol.json')
launch = read(HERE / 'launch_manifest.json')
assert sorted(os.sched_getaffinity(0)) == protocol['resources']['cpus']
assert os.environ.get('CUDA_VISIBLE_DEVICES') is not None
for name, digest in launch['files_sha256'].items():
    assert sha(HERE / name) == digest, name
for name, digest in launch['engine_files_sha256'].items():
    assert sha(Path(protocol['source_engine']) / name) == digest, name
for name, digest in protocol['source_sha256'].items():
    assert sha(HERE / name) == digest, name
rows = {r['key']: r for r in read(HERE / 'inputs.json')}
rows = [rows[key] for key in protocol['expected_keys']]
for row in rows:
    assert sha(HERE / row['coordinate_file']) == row['coordinate_sha256']
write(HERE / 'predictor_owner.json', owner())
import numpy as np
import torch
from ase import Atoms
torch.set_num_threads(1)
torch.set_num_interop_threads(1)
from metag.energetics.uma import load_uma, _predict
calculator = load_uma(protocol['inference']['model'], profile=protocol['inference']['profile'])
write(HERE / 'model_provenance.json', dict(producer=calculator.metag_provenance,
    runtime=calculator.metag_runtime_profile, torch=torch.__version__, numpy=np.__version__))
def atoms(row):
    a = Atoms(row['symbols'], positions=row['positions_A'])
    a.info.update(charge=0, spin=1)
    return a

status = dict(status='running', expected=len(rows), completed=0, failures=[], completed_batches=[])
write(HERE / 'gpu_status.json', status)
started = time.monotonic()
consecutive = 0
(HERE / 'predictions').mkdir(exist_ok=True)
(HERE / 'attempts').mkdir(exist_ok=True)
(HERE / 'claims').mkdir(exist_ok=True)
checks = {rows[i]['key'] for i in [0, len(rows)//2, len(rows)-1]}
try:
    size = protocol['inference']['batch_size']
    for index in range(0, len(rows), size):
        assert time.monotonic() - started < protocol['resources']['seconds'], 'budget exhausted'
        batch = rows[index:index+size]
        slot = f'{index//size:04d}'
        claim = HERE / 'claims' / slot
        claim.mkdir()
        write(claim / 'owner.json', owner())
        attempt = HERE / 'attempts' / slot
        attempt.mkdir()
        write(attempt / 'owner.json', owner())
        try:
            assert not (HERE / 'predictions' / (slot + '.json')).exists()
            print('PREDICT START', slot, batch[0]['path'], flush=True)
            energies, force, bi = _predict(calculator, [atoms(row) for row in batch])
            energies = energies.cpu().numpy()
            force = force.cpu().numpy().astype(float)
            records = []
            offset = 0
            for row, energy in zip(batch, energies):
                n = len(row['symbols'])
                record = dict(key=row['key'], path=row['path'], coordinate_sha256=row['coordinate_sha256'],
                    energy_eV=float(energy), forces_eV_A=force[offset:offset+n].tolist())
                assert np.isfinite(energy) and np.isfinite(force[offset:offset+n]).all()
                if row['key'] in checks:
                    individual, _, _ = _predict(calculator, [atoms(row)])
                    diff = (float(individual.cpu().numpy()[0]) - float(energy)) * protocol['eV_to_kj']
                    record['individual_minus_batch_kj'] = diff
                    assert abs(diff) <= protocol['inference']['max_individual_batch_difference_kj'], diff
                records.append(record)
                offset += n
            assert offset == len(force)
            write(HERE / 'predictions' / (slot + '.json'), dict(slot=slot,
                producer_sha256=sha(__file__), protocol_sha256=sha(HERE/'protocol.json'),
                records=records, elapsed_seconds=time.monotonic()-started))
            status['completed'] += len(batch)
            status['completed_batches'].append(slot)
            consecutive = 0
            print('PREDICT COMPLETE', slot, status['completed'], '/', len(rows), flush=True)
        except Exception as error:
            consecutive += 1
            failure = dict(slot=slot, keys=[r['key'] for r in batch], error=str(error), traceback=traceback.format_exc())
            write(attempt/'failure.json', failure)
            status['failures'].append(failure)
            if consecutive >= protocol['resources']['consecutive_failure_limit']:
                raise
        finally:
            shutil.rmtree(claim)
            write(HERE/'gpu_status.json', status)
    status['status'] = 'complete' if status['completed'] == len(rows) else 'terminal_with_failures'
except BaseException as error:
    status.update(status='stopped', error=str(error), traceback=traceback.format_exc())
    raise
finally:
    status['elapsed_seconds'] = time.monotonic() - started
    write(HERE/'gpu_status.json', status)
