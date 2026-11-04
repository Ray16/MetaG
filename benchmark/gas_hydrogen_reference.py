"""Parent-DFT gas H2 reference for a closed thermochemical redox diagnostic.

The calculation uses no experimental reaction labels and changes no MetaG
prediction. A gas reference is kept separate from dissolved-hydrogen estimates.
"""
import argparse
import datetime
import hashlib
import json
import math
from pathlib import Path
import sys

from audit_modelseed_inputs import pin_one_cpu


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():parser.error('output directory already exists')
    resources=pin_one_cpu()
    import numpy as np
    import pyscf
    from pyscf import dft,gto,lib
    from scipy.optimize import minimize_scalar
    from ase import Atoms
    from ase.thermochemistry import IdealGasThermo
    import ase
    lib.num_threads(1)
    args.output.mkdir()
    evaluations=[]
    protocol={'schema':'metag-h2-gas-reference-v1','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'method':'restricted singlet omegaB97M-V/def2-TZVPD; gas phase',
              'charge':0,'spin_2S':0,'T_K':298.15,'pressure_Pa':[101325.,100000.],
              'orbital_grid_level':4,'nonlocal_grid_level':3,'refined_grid_levels':[5,4],
              'scf_convergence_Hartree':1e-12,'geometry_bounds_Angstrom':[.55,1.05],
              'geometry_tolerance_Angstrom':1e-5,'curvature_steps_Angstrom':[.005,.0025],
              'thermal':'harmonic vibration; linear rotor; symmetry number 2; electronic spin S=0',
              'scope':'An independently calculated primitive reference, not an experimental benchmark, fitted correction or estimator change.',
              'pyscf_version':pyscf.__version__,'ase_version':ase.__version__,'resources':resources,
              'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (args.output/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')

    def energy(distance,grid=4,nlc_grid=3):
        molecule=gto.M(atom=[['H',(0.,0.,-distance/2)],['H',(0.,0.,distance/2)]],
                       basis='def2-tzvpd',charge=0,spin=0,verbose=0)
        mf=dft.RKS(molecule);mf.xc='wb97m-v';mf.grids.level=grid;mf.nlcgrids.level=nlc_grid
        mf.conv_tol=1e-12;mf.max_cycle=100
        assert mf.do_nlc(),'Nonlocal VV10 correlation is required'
        value=float(mf.kernel())
        if not mf.converged or not math.isfinite(value):raise RuntimeError('H2 reference SCF did not converge')
        row={'distance_Angstrom':float(distance),'grid':grid,'nlc_grid':nlc_grid,'energy_Hartree':value,
             'n_AO':molecule.nao_nr(),'converged':bool(mf.converged)}
        evaluations.append(row)
        with (args.output/'evaluations.jsonl').open('a') as stream:stream.write(json.dumps(row)+'\n')
        return value

    optimum=minimize_scalar(energy,bounds=(.55,1.05),method='bounded',options={'xatol':1e-5})
    if not optimum.success:raise RuntimeError('H2 reference geometry optimization failed')
    distance=float(optimum.x);base=energy(distance);refined=energy(distance,5,4)
    curvatures=[]
    for step in [.005,.0025]:
        curvature=(energy(distance+step)-2*base+energy(distance-step))/(step*step)
        if not math.isfinite(curvature) or curvature<=0:raise RuntimeError('H2 reference has no positive bond curvature')
        curvatures.append(curvature)
    atoms=Atoms('H2',positions=[[0,0,-distance/2],[0,0,distance/2]])
    reduced_mass_kg=float(atoms.get_masses()[0])/2*1.66053906660e-27
    omega=math.sqrt(curvatures[-1]*4.3597447222071e-18/1e-20/reduced_mass_kg)
    wavenumber=omega/(2*math.pi*2.99792458e10)
    vibration_eV=omega*6.582119569e-16
    thermo=IdealGasThermo(vib_energies=[vibration_eV],potentialenergy=0.,atoms=atoms,
                         geometry='linear',symmetrynumber=2,spin=0)
    hartree_kj=2625.4996394799;ev_kj=96.485
    gas={}
    for pressure in [101325.,100000.]:
        correction=float(thermo.get_gibbs_energy(temperature=298.15,pressure=pressure,verbose=False))*ev_kj
        gas[str(pressure)]={'electronic_kj_per_mol':base*hartree_kj,
                            'harmonic_RRHO_correction_kj_per_mol':correction,
                            'G_gas_kj_per_mol':base*hartree_kj+correction}
    report={'schema':'metag-h2-gas-reference-v1','protocol':protocol,'distance_Angstrom':distance,
            'electronic_energy_Hartree':base,'refined_energy_Hartree':refined,
            'refined_grid_energy_change_kj_per_mol':(refined-base)*hartree_kj,
            'bond_curvatures_Hartree_per_Angstrom2':curvatures,
            'curvature_step_relative_difference':abs(curvatures[0]-curvatures[1])/curvatures[-1],
            'harmonic_wavenumber_cm_minus1':wavenumber,'gas_standard_states':gas,
            'n_electronic_evaluations':len(evaluations),
            'scope':protocol['scope']+' Anharmonic and nuclear-spin refinements are not included; the deployed harmonic convention is explicit.'}
    if abs(report['refined_grid_energy_change_kj_per_mol'])>.05:raise RuntimeError('Reference grid sensitivity exceeds 0.05 kJ/mol')
    if report['curvature_step_relative_difference']>.01:raise RuntimeError('Reference bond curvature is not step-stable')
    with (args.output/'result.json').open('x') as stream:json.dump(report,stream,indent=2);stream.write('\n')
    print(json.dumps({k:v for k,v in report.items() if k!='protocol'},indent=2))


if __name__=='__main__':main()
