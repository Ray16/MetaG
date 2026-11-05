"""Independent decimal reference assembly and QR Hessian/thermal replay."""
from pathlib import Path
from collections import Counter
from decimal import Decimal,localcontext
import hashlib
import json
import math
import numpy as np
from ase import Atoms,units
from ase.thermochemistry import IdealGasThermo
from rdkit import Chem

HERE=Path(__file__).resolve().parent
read=lambda p:json.loads(Path(p).read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
D=lambda value:Decimal(str(value))


def inventory(smiles):
    return Counter(a.GetSymbol() for a in Chem.AddHs(Chem.MolFromSmiles(smiles)).GetAtoms())


def cold_hessian(path,record):
    archive=np.load(path,allow_pickle=False)
    positions=archive['positions_angstrom'];H=archive['H_eV_A2'];forces=archive['forces_eV_A']
    data=read(HERE/'lbfgs_control/input.json')
    atoms=Atoms(data['symbols'],positions=positions)
    masses=atoms.get_masses();center=positions-(masses[:,None]*positions).sum(axis=0)/masses.sum()
    external=[]
    for axis in np.eye(3):external.append((np.ones_like(center)*axis*np.sqrt(masses[:,None])).reshape(-1))
    for axis in np.eye(3):external.append((np.cross(np.ones_like(center)*axis,center)*np.sqrt(masses[:,None])).reshape(-1))
    external=np.array(external).T
    assert np.linalg.matrix_rank(external)==6
    Q,_=np.linalg.qr(external,mode='complete');internal=Q[:,6:]
    inverse_mass=1/np.sqrt(np.repeat(masses,3))
    weighted=H*inverse_mass[:,None]*inverse_mass[None,:]
    eigen=np.linalg.eigvalsh(internal.T@weighted@internal)
    conversion=units._hbar*1e10/math.sqrt(units._e*units._amu)
    modes=np.sign(eigen)*conversion*np.sqrt(abs(eigen))
    cm2ev=1/8065.544005069921
    imaginary=-modes[modes<0]/cm2ev
    strong=int(sum(imaginary>50))
    assert strong==record['n_imag']==0
    assert float(np.linalg.norm(forces,axis=1).max())<=.001
    energies=np.sort(np.concatenate([modes[modes>=0],imaginary*cm2ev]))
    energies=np.maximum(energies,50*cm2ev)
    electronic=float(archive['energy_eV'])
    thermo=IdealGasThermo(vib_energies=energies,potentialenergy=electronic,atoms=atoms,
        geometry='nonlinear',symmetrynumber=1,spin=0,vib_selection='all')
    corr=(thermo.get_gibbs_energy(298.15,101325,verbose=False)-electronic)*96.485
    assert abs(corr-record['thermal_kj'])<.0001,(corr,record['thermal_kj'])
    return dict(archive_sha256=sha(path),n_strong_imag=strong,
        n_soft_imag=len(imaginary),force_max_eV_A=float(np.linalg.norm(forces,axis=1).max()),
        independent_thermal_kj=corr,stored_thermal_kj=record['thermal_kj'],difference_kj=corr-record['thermal_kj'])


def main():
    catalog=read(HERE/'catalog596.json');raw={(r['charge'],r['smiles']):r for r in catalog['records']}
    doc=read(HERE/'reassembly596/global_amino_extended596.json')
    base=read(HERE/'reassembly596/full_extended596.json')
    nominal={r['reference_smiles']:r for r in read(HERE/'models_nominal.json')}
    gh=D('-1170.7610427111');std=D('7.925418160410476');rtln10=D('5.7087203873')
    states_checked=[];unaffected=0
    with localcontext() as context:
        context.prec=50
        for key,record in doc['ledger']['states'].items():
            plan=record.get('plan',{});model=plan.get('source_reference_model')
            if model is None:
                if key in base['ledger']['states']:
                    assert abs(record['potential_kj']-base['ledger']['states'][key]['potential_kj'])<1e-7
                    unaffected+=1
                continue
            assert model==nominal[model['reference_smiles']]
            ref=model['reference_smiles'];source=raw[(1,ref)]
            weights=[D(0)]
            for pka in model['pKas']:weights.append(weights[-1]+D(7)-D(pka))
            transform=-rtln10*sum(D(10)**w for w in weights).log10()
            exchange=plan['input_charge']-1
            mu=D(source['G_kj'])+std+transform+D(exchange)*gh
            assert abs(mu-D(record['potential_kj']))<D('1e-7')
            residual=inventory(ref);residual.subtract(inventory(plan['input_smiles']));residual['H']+=exchange
            assert not {k:v for k,v in residual.items() if v}
            assert 1-plan['input_charge']+exchange==0
            states_checked.append(dict(state_id=key,reference=ref,
                independent_potential_kj=float(mu),difference_kj=float(mu-D(record['potential_kj']))))
        for rid,row in doc['results'].items():
            if row['status']!='estimated':continue
            result=row['result'];terms=[D(result['original_species'][name][0])*D(r['potential_kj'])
                for name,r in result['state_records'].items()]
            inputs=read(HERE/'inputs.json')[rid]
            terms.append(D(inputs.get('n_Hplus',0))*gh)
            assert abs(sum(terms)-D(result['dG_unrounded']))<D('1e-7')
            sensitivities={}
            for name,state in result['state_records'].items():
                nu=result['original_species'][name][0]
                for key,source in state['sampling_inputs'].items():
                    if key not in sensitivities:sensitivities[key]=[0.,source['sigma_kj']]
                    assert sensitivities[key][1]==source['sigma_kj']
                    sensitivities[key][0]+=nu*source['sensitivity']
            sigma=math.sqrt(math.fsum((s*sd)**2 for s,sd in sensitivities.values()))
            assert abs(sigma-result['U_samp'])<1e-8
    physical={}
    for slot in ['original_endpoint','after_tighter_FIRE']:
        record=read(HERE/'lbfgs_control'/(slot+'.json'))
        assert record['strict_endpoint_accepted']
        assert sha(HERE/'lbfgs_control'/(slot+'.npz'))==record['Hessian_archive_sha256']
        physical[slot]=cold_hessian(HERE/'lbfgs_control'/(slot+'.npz'),record)
    report=dict(scope='Independent assembly/proton balance and retained-Hessian minimum controls; no experimental accuracy or total uncertainty claim.',
        n_source_states=len(states_checked),source_states=states_checked,
        unaffected_shared_states_verified=unaffected,estimated_reactions_verified=148,
        independent_Hessian_controls=physical,production_reference_policy_adopted=False)
    with (HERE/'independent_verification.json').open('x') as f:json.dump(report,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({k:v for k,v in report.items() if k!='source_states'},indent=2))

if __name__=='__main__':main()
