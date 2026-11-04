"""Rebind the original electronic energies to exact stored geometries and producer."""
import json
import math
import time
from pathlib import Path
from control_common import HERE, authenticate, owner, write

start = time.monotonic()
inputs, protocol = authenticate()
import torch
torch.set_num_threads(1)
torch.set_num_interop_threads(1)
from ase import Atoms
from metag.energetics.uma import load_uma, batched_energies

write(HERE/'uma_owner.json', owner())
calculator = load_uma('uma-s-1p2p1', profile='nitro')
for source in inputs['source_records'].values():
    raw = json.loads(Path(source['path']).read_text())
    assert calculator.metag_provenance == raw['settings']['producer']['model']
rows = {}
for key, item in inputs['inputs'].items():
    geometry = item['geometry']
    atoms = Atoms(geometry['symbols'], positions=geometry['positions_angstrom'])
    atoms.info.update(charge=geometry['charge'], spin=geometry['spin_multiplicity'])
    value = float(batched_energies(calculator, [atoms], chunk=1)[0]) * protocol['UMA_eV_to_kj_per_mol']
    original = item['original_components_kj']['electronic_kj']
    difference = value - original
    assert math.isfinite(value)
    rows[key] = dict(geometry_sha256=geometry['sha256'], replay_electronic_kj=value,
                     original_electronic_kj=original, difference_kj=difference,
                     passed=abs(difference) <= protocol['UMA_replay_max_difference_kj'])
    write(HERE/'uma_replay_partial.json', dict(rows=rows, producer=calculator.metag_provenance))
    print(f'UMA replay {key}: delta={difference:.8g} kJ/mol', flush=True)
    assert rows[key]['passed'], key
write(HERE/'uma_replay.json', dict(rows=rows, producer=calculator.metag_provenance,
      runtime_profile=calculator.metag_runtime_profile,
      elapsed_seconds=time.monotonic()-start, passed=True))
