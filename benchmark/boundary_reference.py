"""Audit the numerical boundary panel; these are model estimates, not accuracy labels."""
from decimal import Decimal,localcontext
from pathlib import Path
import argparse
import gzip
import hashlib
import json
import math
import os
import subprocess
import sys
import tempfile
import zipfile
from collections import defaultdict
from fractions import Fraction

read=lambda p:json.loads(Path(p).read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()


def manifest(packet):
    expected=read(packet/'manifest.json')['files']
    actual={str(p.relative_to(packet)) for p in packet.rglob('*') if p.is_file() and '__pycache__' not in p.parts}
    assert actual==set(expected)|{'manifest.json'}
    for name,row in expected.items():
        assert not Path(name).is_absolute() and '..' not in Path(name).parts
        path=packet/name;assert not path.is_symlink() and path.stat().st_size==row['bytes'] and sha(path)==row['sha256']


def extract(path,destination):
    with zipfile.ZipFile(path) as archive:
        entries=archive.infolist();assert len({e.filename for e in entries})==len(entries)
        assert sum(e.file_size for e in entries)<150_000_000
        for entry in entries:
            name=Path(entry.filename)
            assert not name.is_absolute() and '..' not in name.parts and (entry.external_attr>>16)&0o170000!=0o120000
            target=destination/name;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(archive.read(entry))


def audit(packet):
    from rdkit import Chem
    def atoms(smiles):
        mol=Chem.MolFromSmiles(smiles);assert mol is not None
        counts=defaultdict(Fraction)
        for a in Chem.AddHs(mol).GetAtoms():counts[(a.GetAtomicNum(),a.GetIsotope())]+=1
        return counts,Chem.GetFormalCharge(mol)
    packet=Path(packet);manifest(packet)
    records=read(packet/'raw_catalog.json')['records'];requests=read(packet/'requests.json')
    raw={};maximum_partition=0.;point_count=0
    with localcontext() as context:
        context.prec=50;D=lambda x:Decimal(str(x));rt=D('2.4789')
        with zipfile.ZipFile(packet/'raw_records.zip') as archive:
            assert set(archive.namelist())=={r['name'] for r in records}
            for record in records:
                data=archive.read(record['name']);assert hashlib.sha256(data).hexdigest()==record['sha256']
                row=json.loads(data);key=(record['charge'],record['smiles']);assert key not in raw
                assert (row['q'],row['smi'])==key
                assert atoms(row['smi'])[1]==row['q']
                assert row['G']==record['G_kj'] and row['sigma']==record['sigma_kj']
                ensemble=row['meta']['thermal']['ensemble'];points=ensemble['validated_minima']
                assert points and ensemble['n_reference_thermal_approximation']==0
                values=[]
                for point in points:
                    g=point['geometry'];body=dict(g);digest=body.pop('sha256')
                    assert hashlib.sha256(json.dumps(body,sort_keys=True,allow_nan=False).encode()).hexdigest()==digest
                    assert g['charge']==row['q']
                    value=sum(D(point[k]) for k in ['electronic_kj','solvation_kj','thermal_kj'])
                    assert abs(value-D(point['total_kj']))<D('1e-7');values.append(value)
                low=min(values);total=low-rt*sum((-(v-low)/rt).exp() for v in values).ln()
                difference=abs(total-D(row['G']));assert difference<D('1e-7')
                maximum_partition=max(maximum_partition,float(difference));point_count+=len(points);raw[key]=row
        std=D('.008314')*D('298.15')*(D('.082057')*D('298.15')).ln()
        hp=D('-26.3')-D('1104.5')-D('2.303')*D('.008314')*D('298.15')*D(7)
        scale=D('2.303')*D('.008314')*D('298.15')
        reports={};maximum_reaction=0.
        with zipfile.ZipFile(packet/'assemblies.zip') as archive:
            for phase in ['primary','augmented']:
                report=json.loads(archive.read(phase+'/report.json'));reports[phase]=report
                counts={}
                for arm in ['full','anchor14','anchor12and14']:
                    document=json.loads(archive.read(phase+'/'+arm+'.json'));rows=document['results'];assert set(rows)==set(requests)
                    counts[arm]=sum(r['status']=='estimated' for r in rows.values())
                    for rid,row in rows.items():
                        if row['status']!='estimated':continue
                        result=row['result'];rx=requests[rid]['request'];parts=[D(rx.get('n_Hplus',0))*hp];uncertainty={}
                        for name,state in result['state_records'].items():
                            plan=state['plan'];nu=D(rx['species'][name][0]);q,smi=plan['reference_charge'],plan['reference_smiles']
                            original,oq=atoms(rx['species'][name][2]);reference,rq=atoms(smi)
                            assert oq==rx['species'][name][1] and rq==q
                            exchange=plan['proton_exchange']
                            assert q+exchange==oq
                            for element in set(original)|set(reference):
                                assert reference.get(element,0)+(exchange if element==(1,0) else 0)==original.get(element,0)
                            # These fixed inputs have no hydrated or declared ring-state pools.
                            assert state['coupled_protonation'] is None
                            transfer=state.get('molecular_diagnostics',{}).get('anchored_transfer')
                            if smi=='O' and q==0:
                                molecular=D(read(packet/'water_reference.json')['water_liquid_1atm_kj'])
                            else:
                                terms=transfer['raw_terms'] if transfer else [dict(coefficient=1,charge=q,smiles=smi)]
                                composed=defaultdict(Fraction);charge=Fraction(0)
                                for term in terms:
                                    elements,tq=atoms(term['smiles']);assert tq==term['charge'];w=Fraction(str(term['coefficient']))
                                    charge+=w*tq
                                    for element,count in elements.items():composed[element]+=w*count
                                assert charge==q and all(composed[k]==reference.get(k,0) for k in set(composed)|set(reference))
                                molecular=sum(D(t['coefficient'])*D(raw[(t['charge'],t['smiles'])]['G']) for t in terms)
                            assert abs(molecular+std-D(state['solution_reference_kj']))<D('1e-7')
                            phi=D(0)
                            for site in plan['pka_sites']:
                                side,pka,kind,mult=site;assert side=='prod'
                                exponent=D(7)-D(pka) if kind=='acid' else D(pka)-D(7)
                                log10=(D(1)+(D(10)**exponent)).ln()/D(10).ln()
                                phi-=D(mult)*scale*log10
                            assert abs(phi-D(state['protonation_kj']))<D('1e-8')
                            potential=molecular+std+phi+D(plan['proton_exchange'])*hp
                            assert abs(potential-D(state['potential_kj']))<D('1e-7')
                            parts.append(nu*potential)
                            for identity,p in state['sampling_inputs'].items():
                                assert p['sigma_kj']==raw[(p['charge'],p['smiles'])]['sigma']
                                entry=uncertainty.setdefault(identity,[D(0),D(p['sigma_kj'])]);entry[0]+=nu*D(p['sensitivity'])
                        error=abs(sum(parts)-D(result['dG_unrounded']));assert error<D('1e-7'),(phase,arm,rid,error)
                        maximum_reaction=max(maximum_reaction,float(error))
                        sigma=sum((v[0]*v[1])**2 for v in uncertainty.values()).sqrt()
                        assert abs(sigma-D(result['U_samp']))<D('1e-8')
                        assert abs(sum(D(v) for v in row['components_kj'].values())-D(result['dG_unrounded']))<D('1e-7')
                    for cycle in report['cycles'][arm]:
                        stoichiometry=defaultdict(Fraction)
                        for rid,n in cycle['coefficients'].items():
                            rx=requests[rid]['request']
                            for nu,q,smi in rx['species'].values():
                                key=(q,Chem.MolToSmiles(Chem.MolFromSmiles(smi),isomericSmiles=True));stoichiometry[key]+=n*Fraction(str(nu))
                            stoichiometry['H+']+=n*Fraction(str(rx.get('n_Hplus',0)))
                        assert not any(stoichiometry.values())
                        residual=sum(D(n)*D(rows[rid]['result']['dG_unrounded']) for rid,n in cycle['coefficients'].items())
                        assert abs(residual)/sum(abs(n) for n in cycle['coefficients'].values())<D('1e-7')
                assert counts=={arm:report['assembly_status_counts'][arm].get('estimated',0) for arm in counts}
                inv=json.loads(archive.read(phase+'/invariances.json'))
                assert len(inv)==report['invariance_checks'] and all(abs(r['residual_kj'])<1e-7 for r in inv)
                assert not report['production_policy_adopted'] and report['primary_independent_experimental_labels']==0
    curvature=curvature_audit(packet)
    return dict(status='passed',raw_molecular_records=len(raw),retained_minima=point_count,
        maximum_decimal_partition_error_kj=maximum_partition,maximum_decimal_reaction_error_kj=maximum_reaction,
        availability={phase:{a:r['assembly_status_counts'][a].get('estimated',0) for a in r['assembly_status_counts']} for phase,r in reports.items()},
        augmented_invariance_checks=reports['augmented']['invariance_checks'],
        augmented_cycles_per_arm={a:len(v) for a,v in reports['augmented']['cycles'].items()},
        curvature_controls=curvature,independent_experimental_labels=0,production_policy_adopted=False)


def curvature_audit(packet):
    import numpy as np
    from ase import Atoms,units
    from ase.thermochemistry import IdealGasThermo
    differences=[];forces=[]
    for path in sorted((packet/'curvature_controls').glob('*.json')):
        row=read(path);matrix=path.with_suffix('.npz');assert sha(matrix)==row['matrix_sha256']
        with np.load(matrix) as data:H=data['hessian_eV_A2']
        g=row['geometry'];atoms=Atoms(symbols=g['symbols'],positions=g['positions_angstrom']);m=atoms.get_masses()
        xyz=atoms.positions-np.average(atoms.positions,axis=0,weights=m);root=np.sqrt(m);columns=[]
        for axis in np.eye(3):columns.append((np.broadcast_to(axis,xyz.shape)*root[:,None]).ravel())
        for axis in np.eye(3):columns.append((np.cross(axis,xyz)*root[:,None]).ravel())
        external=np.column_stack(columns);assert np.linalg.matrix_rank(external)==6
        Q,_=np.linalg.qr(external,mode='complete');B=Q[:,6:]
        weighted=H/np.outer(np.repeat(root,3),np.repeat(root,3));eig=np.linalg.eigvalsh(B.T@weighted@B)
        energies=units._hbar*1e10/math.sqrt(units._e*units._amu)*np.sqrt(np.abs(eig))
        assert not np.any(energies[eig<0]/1.23984e-4>50)
        energies=np.maximum(energies,50*1.23984e-4)
        assert row['geometry_kind']=='nonlinear' and row['symmetrynumber']==1
        thermo=IdealGasThermo(vib_energies=energies,potentialenergy=0.,atoms=atoms,geometry='nonlinear',
            symmetrynumber=1,spin=(g['spin_multiplicity']-1)/2,vib_selection='all')
        value=thermo.get_gibbs_energy(temperature=298.15,pressure=101325.,verbose=False)*96.485
        assert abs(value-row['thermal_kj'])<.0001 and abs(value-row['expected_thermal_kj'])<.02
        differences.append(abs(value-row['expected_thermal_kj']));forces.append(row['force_max_eV_A'])
    assert len(differences)==4
    return dict(n_fixed_geometries=4,maximum_retained_RRHO_difference_kj=max(differences),
        force_range_eV_A=[min(forces),max(forces)],scope='Numerical reconstruction at four retained geometries; not all-minimum replay, exact stationarity, complete aqueous ensembles or experimental validation.')


def replay(packet):
    packet=Path(packet).resolve()
    with tempfile.TemporaryDirectory(prefix='metag-boundary-replay-') as name:
        root=Path(name);extract(packet/'scoring_source.zip',root)
        extract(packet/'assemblies.zip',root/'expected')
        script=(packet/'cold_replay.py').read_text();(root/'cold_replay.py').write_text(script)
        env=dict(os.environ,OMP_NUM_THREADS='1',OMP_THREAD_LIMIT='1',MKL_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
        result=subprocess.run([sys.executable,str(root/'cold_replay.py'),str(packet)],cwd=root,env=env,capture_output=True,text=True,timeout=120)
        if result.returncode:raise RuntimeError(result.stdout+result.stderr)
        return json.loads(result.stdout.strip().splitlines()[-1])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--packet',type=Path,default=Path(__file__).with_name('boundary_reference_validation'))
    parser.add_argument('--cold-replay',action='store_true')
    args=parser.parse_args();result=audit(args.packet)
    if args.cold_replay:result['cold_replay']=replay(args.packet)
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
