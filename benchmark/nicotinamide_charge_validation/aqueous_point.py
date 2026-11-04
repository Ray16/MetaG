"""One fixed-geometry SMD point, with the parent numerical settings unchanged."""
import argparse
import importlib.metadata
import math
import time
from control_common import HERE,authenticate,read,owner,sha,write

parser=argparse.ArgumentParser();parser.add_argument('slot');args=parser.parse_args()
inputs,protocol,numeric=authenticate()
import numpy as np
import cupy as cp
from pyscf import gto,lib
from gpu4pyscf import dft,solvent
from gpu4pyscf.scf import hf as gpu_hf
gpu_hf.overlap_zero_eigenvalue_threshold=numeric['overlap_minimum_eigenvalue']
lib.num_threads(1)
assert cp.cuda.runtime.getDeviceCount()==1
out=HERE/'points';out.mkdir(exist_ok=True)
if args.slot=='smoke':
    fixture=read(HERE/'cpu_water_smoke.json')
    assert fixture['source_sha256']==sha(HERE/'cpu_smoke.py')
    geometry=dict(symbols=[r[0] for r in fixture['geometry']],positions_angstrom=[r[1] for r in fixture['geometry']],charge=0,spin_multiplicity=1)
else:geometry=inputs[args.slot]['geometry']
start=time.monotonic()
mol=gto.M(atom=list(zip(geometry['symbols'],geometry['positions_angstrom'])),unit='Angstrom',
    basis=numeric['basis'],charge=geometry['charge'],spin=geometry['spin_multiplicity']-1,
    max_memory=numeric['max_memory_MB'],verbose=4,output=str(out/(args.slot+'.scf.log')))
smallest=float(np.linalg.eigvalsh(mol.intor('int1e_ovlp'))[0])
assert smallest>=numeric['overlap_minimum_eigenvalue'],smallest
if args.slot!='smoke':
    assert mol.nao_nr()==inputs[args.slot]['expected_n_AO']
    assert mol.nelectron==inputs[args.slot]['n_electrons']
mf=dft.RKS(mol,xc=numeric['functional']).density_fit(auxbasis=numeric['auxbasis'])
mf.grids.level=numeric['grid_level'];mf.nlcgrids.level=numeric['nlc_grid_level']
mf.conv_tol=numeric['SCF_tolerance_Eh'];mf.conv_tol_grad=numeric['SCF_gradient_tolerance']
mf.max_cycle=numeric['max_cycle'];mf.chkfile=None;assert mf.do_nlc()
mf=solvent.SMD(mf)
mf.with_solvent.solvent='water';mf.with_solvent.sasa_ng=numeric['SMD_sasa_ng']
mf.with_solvent.equilibrium_solvation=numeric['SMD_equilibrium_solvation']
scalar=lambda value:float(value.get() if hasattr(value,'get') else value)
result=dict(slot=args.slot,geometry=geometry,n_AO=mol.nao_nr(),n_electrons=mol.nelectron,
    overlap_minimum_eigenvalue=smallest,owner=owner(),source_sha256=sha(__file__),
    input_sha256=sha(HERE/'geometry_inputs.json'),protocol_sha256=sha(HERE/'protocol.json'),
    numerical_protocol_sha256=sha(HERE/'numerical_protocol.json'),
    versions={k:importlib.metadata.version(k) for k in ['pyscf','gpu4pyscf-cuda12x','gpu4pyscf-libxc-cuda12x','numpy','cupy-cuda12x']},
    physical_scope='Vertical aqueous restricted-singlet electronic energy; no nuclear optimization, RRHO, pKa, stability or ensemble proof.')
write(out/(args.slot+'.started.json'),result)
print(f'START SMD {args.slot}; charge={mol.charge}, AO={mol.nao_nr()}, electrons={mol.nelectron}',flush=True)
energy=scalar(mf.kernel())
result.update(energy_Eh=energy,converged=bool(mf.converged),n_molecular_orbitals=int(mf.mo_coeff.shape[1]),
    scf_summary={str(k):scalar(v) for k,v in mf.scf_summary.items()})
write(out/(args.slot+'.converged_stage.json'),result)
assert mf.converged and math.isfinite(energy),'SCF failed'
assert mf.mo_coeff.shape[1]==mol.nao_nr(),'Basis directions discarded'
result['orbital_gradient_norm']=scalar(cp.linalg.norm(mf.get_grad(mf.mo_coeff,mf.mo_occ,mf.get_fock())))
assert result['orbital_gradient_norm']<=numeric['stationary_gradient_maximum']
occupied=mf.mo_energy[mf.mo_occ>0];virtual=mf.mo_energy[mf.mo_occ==0]
result.update(HOMO_Eh=scalar(occupied.max()),LUMO_Eh=scalar(virtual.min()),
    CDS_Eh=scalar(mf.with_solvent.e_cds),surface_discretization_method=mf.with_solvent.surface_discretization_method,
    elapsed_seconds=time.monotonic()-start)
result['negative_HOMO_diagnostic']=result['HOMO_Eh']<0
summary=result['scf_summary']
assert abs(summary['e_cds']-result['CDS_Eh'])<1e-12
assert abs(energy-math.fsum(summary[k] for k in ['e1','e2','nuc','e_solvent','e_cds']))<1e-8
if args.slot=='smoke':
    result['CPU_absolute_difference_kj']=abs(energy-fixture['water_total_Eh'])*numeric['Hartree_to_kj_per_mol']
    assert result['CPU_absolute_difference_kj']<=numeric['numerical_smoke_tolerance_kj']
result['passed']=True
write(out/(args.slot+'.json'),result)
print(f'COMPLETE SMD {args.slot}; E={energy:.13f} Eh; HOMO={result["HOMO_Eh"]:.9g} Eh',flush=True)
