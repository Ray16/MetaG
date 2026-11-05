"""Replay raw chemical proton-transfer controls, not transformed pool scores."""
from decimal import Decimal,localcontext
from collections import Counter
from pathlib import Path
import argparse,hashlib,json,math,os,subprocess,sys,tempfile,zipfile


def audit(root):
    root=Path(root);read=lambda p:json.loads(Path(p).read_text());sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    manifest=read(root/'manifest.json')
    for name,digest in manifest['files'].items():
        name=Path(name);assert not name.is_absolute() and '..' not in name.parts
        assert not (root/name).is_symlink() and sha(root/name)==digest,str(name)
    expected=read(root/'results.json');predictions={}
    with tempfile.TemporaryDirectory(prefix='metag-reference-transfer-') as temporary:
        data=Path(temporary)
        with zipfile.ZipFile(root/'calculation_records.zip') as archive:
            names=archive.namelist();assert len(names)==len(set(names))
            assert all(not Path(name).is_absolute() and '..' not in Path(name).parts for name in names)
            for name in names:
                target=data/name;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(archive.read(name))
        # Fresh process isolates module paths and uses only local retained records.
        env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',OMP_NUM_THREADS='1',OMP_THREAD_LIMIT='1',MKL_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1')
        run=subprocess.run([sys.executable,str(data/'independent_audit.py'),'--records-only'],cwd=data,env=env,capture_output=True,text=True,timeout=120)
        if run.returncode:raise ValueError('Native/curvature replay failed: '+run.stderr[-1500:])
        numeric=read(data/'independent_audit.json');assert numeric['passed'] and numeric['n_retained_geometries']==30
        inputs=read(data/'geometry_inputs.json');protocol=read(data/'protocol.json');reference=read(data/'reference.json')
        molecular=read(data/'molecular_status.json')
        expected_starts={f'{family}__{mode}__{seed:02d}' for family in inputs['species'] for mode in inputs['geometry_modes'] for seed in range(4)}
        assert len(molecular['optimized'])==48 and set(molecular['optimized'])==expected_starts
        for slot in expected_starts:
            optimized=read(data/'optimization'/slot/'result.json');assert optimized['passed'] and optimized['slot']==slot
            item=inputs['species'][optimized['species']]
            assert optimized['charge']==item['charge'] and optimized['symbols']==item['symbols']
        for coefficients in inputs['reactions'].values():
            atoms=Counter();charge=0
            for family,nu in coefficients.items():
                for element,count in Counter(inputs['species'][family]['atomic_numbers']).items():atoms[element]+=count*nu
                charge+=nu*inputs['species'][family]['charge']
            assert not any(atoms.values()) and charge==0 and sum(coefficients.values())==0
        states=read(data/'retained_geometry_index.json');requests=read(data/'scf_geometry_inputs.json')['inputs'];point_protocol=read(data/'parent_protocol.json')
        assert len(requests)==30
        points={}
        for slot,item in requests.items():
            curvature=read(data/'curvature'/slot/'result.json');optimized=read(data/'optimization'/slot/'result.json')
            assert curvature['positions_A']==optimized['positions_A']==item['geometry']['positions_angstrom']
            assert curvature['symbols']==item['geometry']['symbols'] and curvature['charge']==item['geometry']['charge']
            point=read(data/'points'/(slot+'.json'));assert point['passed'] and point['slot']==slot and point['geometry']==item['geometry']
            assert point['source_sha256']==sha(data/'reference_point.py')
            assert point['protocol_sha256']==sha(data/'parent_protocol.json') and point['input_sha256']==sha(data/'scf_geometry_inputs.json')
            assert point['n_AO']==item['n_AO'] and point['n_electrons']==item['n_electrons']
            assert point['overlap_minimum_eigenvalue']>=point_protocol['overlap_minimum_eigenvalue']
            for phase in ['gas','water']:
                result=point[phase];parts=result['scf_summary']
                assert result['converged'] and math.isfinite(result['energy_Eh']) and result['n_molecular_orbitals']==item['n_AO']
                assert 0<=result['orbital_gradient_norm']<=point_protocol['stationary_gradient_maximum']
                assert abs(parts['e2']-parts['coul']-parts['exc'])<1e-8
                summands=[parts[k] for k in ['e1','e2','nuc']]
                if phase=='water':
                    summands.extend([parts['e_solvent'],parts['e_cds']]);assert abs(parts['e_cds']-point['CDS_Eh'])<1e-12
                assert abs(math.fsum(summands)-result['energy_Eh'])<1e-8
            assert abs(point['water']['energy_Eh']-point['gas']['energy_Eh']-point['SMD_delta_Eh'])<1e-10
            points[slot]=point
        with localcontext() as context:
            context.prec=50;D=lambda v:Decimal(str(v));kt=D(protocol['pool_KT_kj']);hartree=D(protocol['units']['Hartree_to_kj']);eV=D(protocol['units']['eV_to_kj'])
            for mode in inputs['geometry_modes']:
                for electronic in ['UMA','DFT']:
                    for solvent in ['ALPB','SMD']:
                        key=f'{mode}_geometry|{electronic}|{solvent}';pools={}
                        for family in inputs['species']:
                            slots=states[f'{family}|{mode}'];values=[]
                            for slot in slots:
                                c=read(data/'curvature'/slot/'result.json');p=points[slot]
                                E=D(c['UMA_energy_eV'])*eV if electronic=='UMA' else D(p['gas']['energy_Eh'])*hartree
                                S=D(c['ALPB_delta_eV'])*eV if solvent=='ALPB' else (D(p['water']['energy_Eh'])-D(p['gas']['energy_Eh']))*hartree
                                G=E+S+D(c['RRHO_kj'])
                                if electronic=='DFT' and solvent=='SMD':G=D(p['water']['energy_Eh'])*hartree+D(c['RRHO_kj'])
                                values.append(G)
                            minimum=min(values);Z=sum((-(g-minimum)/kt).exp() for g in values);pools[family]=minimum-kt*Z.ln()
                            saved=expected['arms'][key]['species'][family]
                            assert saved['slots']==slots and abs(float(pools[family])-saved['G_kj'])<1e-7
                        rows={}
                        for reaction,coefficients in inputs['reactions'].items():
                            prediction=float(sum(D(nu)*pools[family] for family,nu in coefficients.items()))
                            saved=expected['arms'][key]['reactions'][reaction]
                            assert abs(prediction-saved['prediction_kj'])<1e-7
                            assert abs(abs(prediction-reference['chemical_reference_kj'][reaction])-saved['absolute_error_kj'])<1e-7
                            rows[reaction]=prediction
                        predictions[key]=rows
        # Reconstruct references without feeding pKas into any prediction above.
        RTln10=.00831446261815324*298.15*math.log(10)
        for key,acid in [('acetate_transfer','acetate'),('thiol_transfer','mercaptoethanol')]:
            value=RTln10*(reference['sources'][acid]['pKa']-reference['sources']['ethanolamine']['pKa'])
            assert abs(value-reference['chemical_reference_kj'][key])<1e-12
        cpu=read(data/'cpu_verification_summary.json');assert cpu['passed'] and len(cpu['tasks'])==3
        maximum=max(value for row in cpu['tasks'].values() for value in row.values());assert maximum<=.05
        scope=read(data/'quantity_scope_audit.json')
        assert all(row['transformed_pool_reaction_is_identity'] for row in scope['reactions'].values())
        assert expected['primary_independent_metabolic_references_added']==0 and not expected['production_reference_policy_adopted']
        return dict(status='passed',n_optimized_starts=48,n_retained_geometries=30,n_raw_chemical_references=2,
            native_curvature_executions_reconstructed=numeric['native_executions_reconstructed'],predictions_kj=predictions,
            CPU_GPU_maximum_difference_kj=maximum,soft_negative_mode_states=numeric['states_with_negative_modes_below_declared_threshold'],
            primary_independent_metabolic_references_added=0,production_reference_policy_adopted=False)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root',nargs='?',type=Path,default=Path(__file__).resolve().parent/'reference_transfer_validation')
    args=parser.parse_args();print(json.dumps(audit(args.root),indent=2))
