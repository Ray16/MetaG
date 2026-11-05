"""Descriptive equal-target accuracy and exact available-network cycle checks."""
from pathlib import Path
from collections import Counter
from fractions import Fraction
import hashlib
import json
import math
import sys

HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE/'assembly_source'),str(HERE/'source_helpers')]
import sympy as sp
from rdkit import Chem
from metag.constants import CONDITIONS,PH,RT_LN10
from metag.energetics.protonation import macroscopic_transform
from metag.routing.sugar_states import prepare_declared_sugar_inputs
from coupled_source_parameters import candidate_models

read=lambda p:json.loads(Path(p).read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(name,value):
    with (HERE/name).open('x') as f:json.dump(value,f,indent=2,allow_nan=False);f.write('\n')

def metric(errors):
    if not errors:return {'n':0}
    order=sorted(map(abs,errors));n=len(order)
    return dict(n=n,MAE=math.fsum(order)/n,RMSE=math.sqrt(math.fsum(e*e for e in errors)/n),
        bias=math.fsum(errors)/n,median_abs=(order[(n-1)//2]+order[n//2])/2,
        n_abs_error_gt20=sum(e>20 for e in order))

def estimates(doc):return {rid:r['result']['dG_unrounded'] for rid,r in doc['results'].items() if r['status']=='estimated'}

def cycle_report(inputs,values):
    ids=sorted(values);columns=[];keys=set()
    for rid in ids:
        reaction=inputs[rid];col={}
        for name,(nu,q,smiles) in reaction['species'].items():
            if name in reaction.get('species_state_pools',{}):
                rx={'species':{'state':[1,q,smiles]},'species_state_pools':{'state':reaction['species_state_pools'][name]}}
                _,pools=prepare_declared_sugar_inputs(rx)
                states=tuple(sorted(s['smiles'] for s in pools['state']['source_states']))
                key=json.dumps(['pool',q,states])
            else:key=json.dumps(['molecule',q,Chem.MolToSmiles(Chem.MolFromSmiles(smiles),isomericSmiles=True)])
            col[key]=col.get(key,Fraction(0))+Fraction(str(nu));keys.add(key)
        proton=Fraction(str(reaction.get('n_Hplus',0)))
        if proton:
            col['explicit_proton_reservoir']=proton;keys.add('explicit_proton_reservoir')
        columns.append(col)
    keys=sorted(keys)
    matrix=sp.Matrix([[sp.Rational(col.get(k,0)) for col in columns] for k in keys])
    basis=[]
    for vector in matrix.nullspace():
        denominator=math.lcm(*(int(v.q) for v in vector))
        integers=[int(v*denominator) for v in vector]
        divisor=math.gcd(*(abs(v) for v in integers));integers=[v//divisor for v in integers]
        if next(v for v in integers if v)<0:integers=[-v for v in integers]
        residual=math.fsum(c*values[rid] for rid,c in zip(ids,integers))
        norm=sum(map(abs,integers));assert abs(residual)/norm<1e-7
        basis.append(dict(coefficients={rid:c for rid,c in zip(ids,integers) if c},
            residual_kj=residual,L1_normalized_residual_kj=residual/norm))
    return dict(n_available_reactions=len(ids),n_state_rows=len(keys),n_independent_cycles=len(basis),
        maximum_normalized_residual_kj=max([abs(b['L1_normalized_residual_kj']) for b in basis] or [0]),
        cycles=basis,scope='Complete nullspace of the available literal-input/pool stoichiometry with explicit proton reservoir; no experimental accuracy claim.')

def main():
    inputs=read(HERE/'inputs.json');labels=read(HERE/'reference.json')['2.0']['labels']
    paths={name:HERE/'reassembly'/(name+'.json') for name in ['full_parent589','full_extended594','global_amino_extended594']}
    paths.update({name:HERE/'reassembly596'/(name+'.json') for name in ['full_extended596','global_amino_extended596']})
    docs={name:read(path) for name,path in paths.items()};vals={name:estimates(doc) for name,doc in docs.items()}
    reports={}
    for name,doc in docs.items():
        pairs=sorted(set(vals[name])&set(labels))
        reports[name]=dict(n_requested=364,status_counts=dict(Counter(r['status'] for r in doc['results'].values())),
            n_labeled_available=len(pairs),conditional_error=metric([vals[name][rid]-labels[rid]['exp'] for rid in pairs]),
            zero_same_targets=metric([-labels[rid]['exp'] for rid in pairs]))
    base='full_extended596';new='global_amino_extended596'
    common=sorted(set(vals[base])&set(vals[new])&set(labels))
    affected=[];unaffected=[]
    for rid in common:
        touched=any(r.get('plan',{}).get('source_reference_model') for r in docs[new]['results'][rid]['result']['state_records'].values())
        (affected if touched else unaffected).append(rid)
    assert all(abs(vals[base][rid]-vals[new][rid])<1e-7 for rid in unaffected)
    def paired(ids):return dict(n=len(ids),baseline=metric([vals[base][rid]-labels[rid]['exp'] for rid in ids]),
        candidate=metric([vals[new][rid]-labels[rid]['exp'] for rid in ids]),
        zero=metric([-labels[rid]['exp'] for rid in ids]))
    by_ec={}
    for first in sorted({str(inputs[rid]['EC']).split('.')[0] for rid in common}):
        ids=[rid for rid in common if str(inputs[rid]['EC']).split('.')[0]==first];by_ec[first]=paired(ids)
    nominal={m['reference_smiles']:m for m in read(HERE/'models_nominal.json')}
    scenarios=[];ranges={rid:[] for rid in vals[new]}
    transcription=read(HERE/'source_transcription.json')
    import itertools
    for flags in itertools.product([0,1],repeat=5):
        models={m['reference_smiles']:m for m in candidate_models(transcription,CONDITIONS,flags)}
        changes={s:macroscopic_transform(m['pKas'],PH,RT_LN10,reference_deprotonations=0)['G_transform']
            -macroscopic_transform(nominal[s]['pKas'],PH,RT_LN10,reference_deprotonations=0)['G_transform'] for s,m in models.items()}
        outcomes={}
        for rid,value in vals[new].items():
            result=docs[new]['results'][rid]['result'];delta=0.
            for label,record in result['state_records'].items():
                model=record.get('plan',{}).get('source_reference_model')
                if model:delta+=inputs[rid]['species'][label][0]*changes[model['reference_smiles']]
            outcomes[rid]=value+delta;ranges[rid].append(delta)
        scored=sorted(set(outcomes)&set(labels))
        scenarios.append(dict(flags=list(flags),conditional_error=metric([outcomes[rid]-labels[rid]['exp'] for rid in scored])))
    cycles={name:cycle_report(inputs,vals[name]) for name in ['full_extended596','global_amino_extended596']}
    write('cycle_verification.json',cycles)
    report=dict(scope='Development-exposed cached reassembly with two separately fingerprinted fresh molecular prerequisites; no independent experimental validation, full fresh sweep or production policy promotion.',
        protocol_sha256=sha(HERE/'protocol.json'),arms=reports,
        equal_catalog_common_labels=paired(common),source_family_present_common=paired(affected),
        source_family_absent_common=paired(unaffected),common_by_EC=by_ec,
        paired_ids=common,source_family_present_ids=affected,
        gained_requests=sorted(set(vals[new])-set(vals[base])),lost_requests=sorted(set(vals[base])-set(vals[new])),
        source_scale_scenarios=scenarios,maximum_conditional_source_scale_range_kj=max(max(v)-min(v) for v in ranges.values()),
        source_scale_scope='32conditional scale interpretations, not statistical/total-error uncertainty; none selected by error.',
        available_network_cycles={name:{k:v for k,v in d.items() if k!='cycles'} for name,d in cycles.items()},
        source_racemase_control=dict(reaction='rxn00283',predicted_dG_kj=vals[new].get('rxn00283'),
            sampling_sigma_kj=(docs[new]['results']['rxn00283']['result']['U_samp'] if 'rxn00283' in vals[new] else None),
            exact_achiral_equilibrium_dG_kj=0.0,scope='Finite independently sampled enantiomer ensembles; sampling spread is not total uncertainty. No numerical equality was imposed.'),
        new_primary_references_admitted=0,production_reference_policy_adopted=False)
    write('comparison_report.json',report)
    print(json.dumps({k:v for k,v in report.items() if k not in ['source_scale_scenarios','paired_ids','source_family_present_ids']},indent=2))

if __name__=='__main__':main()
