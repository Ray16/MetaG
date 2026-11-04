"""Audit a vertical phosphate-charge/redox coupling without fitting a Gibbs offset.

The four chemically distinct states form an atom/charge-balanced exchange, not
a null stoichiometric cycle. A nonzero result is allowed by thermodynamics.
"""
from collections import Counter
from decimal import Decimal, localcontext
import argparse
import hashlib
import json
import math
from pathlib import Path

def read(path):return json.loads(Path(path).read_text())
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def need(condition,message):
    if not condition:raise ValueError(message)
def close(a,b,tolerance,message):
    need(math.isfinite(a) and math.isfinite(b) and abs(a-b)<=tolerance,message)

def point_checks(point,numeric,n_AO,n_electrons):
    need(point['passed'] is True and point['converged'] is True,'Point not converged/passed')
    need(point['n_AO']==n_AO and point['n_molecular_orbitals']==n_AO,'Basis rank mismatch')
    need(point['n_electrons']==n_electrons,'Electron count mismatch')
    need(point['overlap_minimum_eigenvalue']>=numeric['overlap_minimum_eigenvalue'],'Overlap gate')
    need(0<=point['orbital_gradient_norm']<=numeric['stationary_gradient_maximum'],'Orbital gradient')
    summary=point['scf_summary']
    close(summary['e2'],summary['coul']+summary['exc'],1e-8,'Electronic component sum')
    close(point['CDS_Eh'],summary['e_cds'],1e-12,'CDS mismatch')
    close(point['energy_Eh'],math.fsum(summary[k] for k in ['e1','e2','nuc','e_solvent','e_cds']),1e-8,
          'Aqueous total must include CDS once')
    need(point['surface_discretization_method']=='SWIG','SMD surface mismatch')
    need(point['negative_HOMO_diagnostic']==(point['HOMO_Eh']<0),'HOMO diagnostic mismatch')
    need(math.isfinite(point['HOMO_Eh']) and math.isfinite(point['LUMO_Eh']),'Nonfinite orbitals')

def audit(root):
    root=Path(root);manifest=read(root/'launch_manifest.json')
    for name,digest in manifest['files_sha256'].items():need(sha(root/name)==digest,'Changed owned file: '+name)
    inputs=read(root/'geometry_inputs.json')['states'];protocol=read(root/'protocol.json')
    numeric=read(root/'numerical_protocol.json');parent=read(root/'reference_geometry_inputs.json')['inputs']
    need(set(inputs)=={'oxidized','reduced'},'Unexpected requested states')
    need(sha(root/'numerical_protocol.json')==protocol['parent_protocol_sha256'],'Parent protocol mismatch')
    need(sha(root/'reference_geometry_inputs.json')==protocol['parent_input_sha256'],'Parent input mismatch')
    smoke=read(root/'points/smoke.json');cpu=read(root/'cpu_water_smoke.json')
    point_checks(smoke,numeric,58,10)
    need(smoke['source_sha256']==sha(root/'aqueous_point.py'),'Smoke producer mismatch')
    close(smoke['CPU_absolute_difference_kj'],abs(smoke['energy_Eh']-cpu['water_total_Eh'])*numeric['Hartree_to_kj_per_mol'],1e-12,'Smoke difference mismatch')
    need(smoke['CPU_absolute_difference_kj']<=numeric['numerical_smoke_tolerance_kj'],'CPU/GPU smoke mismatch')
    states={}
    for family,item in inputs.items():
        reference=parent[item['parent_slot']];geometry=item['geometry'];original=reference['geometry']
        removed=item['removed_H_indices'];need(len(set(removed))==2,'Exactly two distinct H atoms required')
        need(all(original['symbols'][i]=='H' for i in removed),'Removed atom is not H')
        keep=[i for i in range(len(original['symbols'])) if i not in removed]
        need(geometry['symbols']==[original['symbols'][i] for i in keep],'Atom ordering changed')
        need(geometry['positions_angstrom']==[original['positions_angstrom'][i] for i in keep],'Retained coordinates changed')
        need(geometry['charge']==original['charge']-2,'Charge conversion mismatch')
        content={k:v for k,v in geometry.items() if k!='sha256'}
        need(hashlib.sha256(json.dumps(content,sort_keys=True,allow_nan=False).encode()).hexdigest()==geometry['sha256'],'Geometry content hash')
        point=read(root/'points'/(family+'.json'))
        point_checks(point,numeric,item['expected_n_AO'],item['n_electrons'])
        need(point['geometry']==geometry,'Computed geometry mismatch')
        need(point['source_sha256']==sha(root/'aqueous_point.py'),'Point producer mismatch')
        need(point['input_sha256']==sha(root/'geometry_inputs.json'),'Point input binding')
        need(point['protocol_sha256']==sha(root/'protocol.json'),'Point protocol binding')
        need(point['numerical_protocol_sha256']==sha(root/'numerical_protocol.json'),'Numerical settings binding')
        path=root/'parent_points'/(item['parent_slot']+'.json')
        need(sha(path)==protocol['parent_reference_points'][item['parent_slot']],'Parent point hash')
        previous=read(path);need(previous['geometry']==original,'Parent geometry mismatch')
        need(previous['source_sha256']==sha(root/'parent_reference_point.py'),'Parent producer binding')
        need(previous['versions']==point['versions'],'Different numerical library versions')
        # Convert the original gas/SMD result to the same single-aqueous-point schema.
        ref=dict(previous,**previous['water'],negative_HOMO_diagnostic=previous['water']['HOMO_Eh']<0)
        point_checks(ref,numeric,reference['n_AO'],reference['n_electrons'])
        states['charged_'+family]=point;states['reference_'+family]=ref
    weights=protocol['weights']
    need(weights=={'charged_reduced':1,'charged_oxidized':-1,'reference_reduced':-1,'reference_oxidized':1},'Endpoint changed')
    residual=Counter();charge=0
    for key,weight in weights.items():
        for symbol,count in Counter(states[key]['geometry']['symbols']).items():residual[symbol]+=weight*count
        charge+=weight*states[key]['geometry']['charge']
    need(not any(residual.values()) and charge==0 and sum(weights.values())==0,'Unbalanced exchange')
    with localcontext() as context:
        context.prec=40;D=lambda x:Decimal(str(x));conversion=D(numeric['Hartree_to_kj_per_mol'])
        coupling=sum(D(weight)*D(states[key]['energy_Eh']) for key,weight in weights.items())*conversion
        paired={family:float((D(states['charged_'+family]['energy_Eh'])-D(states['reference_'+family]['energy_Eh']))*conversion)
                for family in ['oxidized','reduced']}
        components={}
        groups={'electronic_and_nuclear_at_aqueous_density':['e1','e2','nuc'],
                'solvent_electrostatic':['e_solvent'],'solvent_CDS':['e_cds']}
        for name,terms in groups.items():
            components[name]=float(sum(D(weight)*sum(D(states[key]['scf_summary'][t]) for t in terms)
                                       for key,weight in weights.items())*conversion)
    close(math.fsum(components.values()),float(coupling),1e-4,'Coupling component reconstruction')
    close(paired['reduced']-paired['oxidized'],float(coupling),1e-7,'Paired deprotonation identity')
    return dict(schema='metag-nicotinamide-phosphate-charge-audit-v1',passed=True,
        scope='Vertical aqueous electronic-plus-solvent coupling; not Gibbs free energies, macroscopic pKas or a production correction.',
        weights=weights,n_new_points=2,n_parent_points=2,element_residual={},charge_residual=0,molecule_residual=0,
        coupling_kj=float(coupling),paired_vertical_deprotonation_kj=paired,components_kj=components,
        paired_difference_scope='Raw aqueous electronic differences after removing two H+ nuclei; proton reservoir omitted. Neither difference is a deprotonation Gibbs energy.',
        component_scope='Bookkeeping in separately polarized aqueous SCF densities; not independent gas/solvent free-energy components.',
        negative_HOMO_diagnostics={k:v['negative_HOMO_diagnostic'] for k,v in states.items()},
        CPU_GPU_smoke_difference_kj=smoke['CPU_absolute_difference_kj'],
        independent_experimental_labels_admitted=0,empirical_offsets_fitted=0,
        limitations=protocol['limitations']+['An atom-balanced exchange can have nonzero energy; this is not a cycle-closure failure.'])

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('root',type=Path);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();result=audit(args.root);args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
