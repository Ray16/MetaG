"""One restricted-singlet DF-DFT gas/SMD pair on an unchanged stored geometry.

SMD total energy already includes electrostatic and CDS terms. Its difference
from the gas total is the complete solvent increment used by this control.
"""
import argparse
import importlib.metadata
import json
import math
import time
from control_common import HERE, authenticate, owner, sha, write

parser = argparse.ArgumentParser()
parser.add_argument('slot')
args = parser.parse_args()
inputs, protocol = authenticate()
import numpy as np
import cupy as cp
from pyscf import gto, lib
from gpu4pyscf import dft, solvent
lib.num_threads(1)
assert cp.cuda.runtime.getDeviceCount() == 1
start = time.monotonic()
output = HERE/'points'
output.mkdir(exist_ok=True)
if args.slot == 'smoke':
    reference = json.loads((HERE/'cpu_water_smoke.json').read_text())
    assert reference['source_sha256'] == sha(HERE/'cpu_smoke.py')
    geometry = dict(symbols=[x[0] for x in reference['geometry']],
                    positions_angstrom=[x[1] for x in reference['geometry']],
                    charge=0, spin_multiplicity=1)
else:
    geometry = inputs['inputs'][args.slot]['geometry']

mol = gto.M(atom=list(zip(geometry['symbols'],geometry['positions_angstrom'])),
            unit='Angstrom', basis=protocol['basis'], charge=geometry['charge'],
            spin=geometry['spin_multiplicity']-1, verbose=4,
            output=str(output/(args.slot+'.scf.log')),
            max_memory=protocol['max_memory_MB'])
overlap_min = float(np.linalg.eigvalsh(mol.intor('int1e_ovlp'))[0])
assert overlap_min >= protocol['overlap_minimum_eigenvalue'], overlap_min
if args.slot != 'smoke':
    assert mol.nao_nr() == inputs['inputs'][args.slot]['n_AO']
    assert mol.nelectron == inputs['inputs'][args.slot]['n_electrons']

def calculator():
    mf = dft.RKS(mol, xc=protocol['functional']).density_fit(auxbasis=protocol['auxbasis'])
    mf.grids.level = protocol['grid_level']
    mf.nlcgrids.level = protocol['nlc_grid_level']
    mf.conv_tol = protocol['SCF_tolerance_Eh']
    mf.conv_tol_grad = protocol['SCF_gradient_tolerance']
    mf.max_cycle = protocol['max_cycle']
    mf.chkfile = None
    assert mf.do_nlc()
    return mf

def scalar(value):
    return float(value.get() if hasattr(value, 'get') else value)

def stage(mf, **kwargs):
    energy = scalar(mf.kernel(**kwargs))
    gradient = scalar(cp.linalg.norm(mf.get_grad(mf.mo_coeff, mf.mo_occ)))
    assert mf.converged and math.isfinite(energy), 'SCF not converged'
    assert gradient <= protocol['stationary_gradient_maximum'], gradient
    occupied = mf.mo_energy[mf.mo_occ > 0]
    virtual = mf.mo_energy[mf.mo_occ == 0]
    return dict(energy_Eh=energy, converged=bool(mf.converged),
                orbital_gradient_norm=gradient, HOMO_Eh=scalar(occupied.max()),
                LUMO_Eh=scalar(virtual.min()),
                scf_summary={str(k):scalar(v) for k,v in mf.scf_summary.items()})

result = dict(slot=args.slot, geometry=geometry, owner=owner(),
              source_sha256=sha(__file__), protocol_sha256=sha(HERE/'fixed_geometry_protocol.json'),
              input_sha256=sha(HERE/'geometry_inputs.json'),
              versions={key:importlib.metadata.version(key) for key in
                        ['pyscf','gpu4pyscf-cuda12x','gpu4pyscf-libxc-cuda12x','numpy','cupy-cuda12x']},
              n_AO=mol.nao_nr(), n_electrons=mol.nelectron, overlap_minimum_eigenvalue=overlap_min,
              limitation='Converged restricted-singlet fixed-geometry point; no orbital stability proof.')
print(f'START gas {args.slot}, AO={mol.nao_nr()}', flush=True)
gas = calculator()
result['gas'] = stage(gas)
write(output/(args.slot+'.gas.json'), result)
print(f'COMPLETE gas {args.slot}', flush=True)
water = calculator()
water.with_df = gas.with_df
water = solvent.SMD(water)
water.with_solvent.solvent = 'water'
water.with_solvent.sasa_ng = protocol['SMD_sasa_ng']
water.with_solvent.equilibrium_solvation = protocol['SMD_equilibrium_solvation']
print(f'START SMD {args.slot}', flush=True)
result['water'] = stage(water, dm0=gas.make_rdm1())
result['SMD_delta_Eh'] = result['water']['energy_Eh']-result['gas']['energy_Eh']
result['CDS_Eh'] = scalar(water.with_solvent.e_cds)
result['surface_discretization_method'] = water.with_solvent.surface_discretization_method
assert abs(result['water']['scf_summary']['e_cds']-result['CDS_Eh']) < 1e-12
result['elapsed_seconds'] = time.monotonic()-start
if args.slot == 'smoke':
    checks = {key:abs(result[stage_name]['energy_Eh']-reference[key])*protocol['Hartree_to_kj_per_mol']
              for stage_name,key in [('gas','gas_Eh'),('water','water_total_Eh')]}
    result['CPU_GPU_absolute_difference_kj'] = checks
    assert max(checks.values()) <= protocol['numerical_smoke_tolerance_kj'], checks
result['passed'] = True
write(output/(args.slot+'.json'), result)
print(f'COMPLETE gas/SMD {args.slot}; delta={result["SMD_delta_Eh"]:.12g} Eh', flush=True)
