"""Audit released-ensemble reaction sums and developmental risk curves without Torch."""
from collections import defaultdict
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path
import argparse
import csv
import hashlib
import json
import math
import re


def population_summary(values):
    if not values or any(isinstance(v,bool) or not isinstance(v,(float,int)) or not math.isfinite(v) for v in values):
        raise ValueError('Finite nonempty member vector required')
    with localcontext() as context:
        context.prec=50
        numbers=[Decimal(str(v)) for v in values];n=Decimal(len(numbers))
        mean=sum(numbers)/n
        sd=(sum((v-mean)**2 for v in numbers)/n).sqrt()
    return float(mean),float(sd)


def reaction_members(molecules, terms):
    """Preserve the shared checkpoint column before forming any dispersion."""
    if not terms:raise ValueError('Nonempty stoichiometry required')
    sizes={len(molecules[t['molecule']]['members_kj']) for t in terms}
    if len(sizes)!=1 or next(iter(sizes))==0:raise ValueError('Member dimensions differ')
    for term in terms:
        population_summary(molecules[term['molecule']]['members_kj'])
        if isinstance(term['coefficient'],bool) or not math.isfinite(term['coefficient']):raise ValueError('Invalid extent')
    return [math.fsum(t['coefficient']*molecules[t['molecule']]['members_kj'][j] for t in terms)
            for j in range(next(iter(sizes)))]


def audit(packet, *, identity_replay=False):
    packet=Path(packet)
    read=lambda name:json.loads((packet/name).read_text())
    sha=lambda path:hashlib.sha256(Path(path).read_bytes()).hexdigest()
    manifest=read('manifest.json')
    actual={str(p.relative_to(packet)) for p in packet.rglob('*') if p.is_file() and '__pycache__' not in p.parts}
    assert actual==set(manifest['files'])|{'manifest.json'}
    for name,entry in manifest['files'].items():
        path=packet/name
        assert not Path(name).is_absolute() and '..' not in Path(name).parts and not path.is_symlink()
        assert sha(path)==entry['sha256'] and path.stat().st_size==entry['bytes'],name
    raw=read('raw_outputs.json');inputs=read('native_inputs.json');report=read('report.json')
    protocol=read('protocol.json');development=read('development_inputs.json')
    assert raw['n_requested']==raw['n_predicted']==len(raw['raw_reactions'])==len(inputs['rows'])==69
    assert not raw['refused'] and not raw['CUDA_initialized'] and raw['device']=='cpu'
    assert raw['experimental_validation_regime']=='development_and_potential_training_exposed'
    assert raw['new_primary_references_admitted']==0 and not raw['independently_validated_risk_improvement']
    assert len(raw['members'])==100 and len({x['file'] for x in raw['members']})==100
    assert all(len(row['members_kj'])==100 for row in raw['molecular_vectors'].values())
    for member in raw['members']:
        assert member['sha256']==protocol['source_assets']['dGbyG/network/model_params/MPNN_model_params/'+member['file']]['sha256']
    originals={r['reaction']:r for r in development['rows']}
    target_rows={r['reaction']:r for r in inputs['rows']}
    assert set(originals)==set(target_rows)=={r['reaction'] for r in raw['raw_reactions']}
    max_reconstruction=0.;max_summary=0.
    for row in raw['raw_reactions']:
        rid=row['reaction'];original=originals[rid];target=target_rows[rid]
        vector=reaction_members(raw['molecular_vectors'],row['source_terms'])
        native=row['members_native_kj'];assert len(native)==100
        max_reconstruction=max(max_reconstruction,max(abs(a-b) for a,b in zip(vector,native)))
        for j,(a,b) in enumerate(zip(vector,native)):
            magnitude=math.fsum(abs(t['coefficient']*raw['molecular_vectors'][t['molecule']]['members_kj'][j]) for t in row['source_terms'])
            assert abs(a-b)<=max(1e-6,(len(row['source_terms'])+2)*2**-23*magnitude)
        mean,sd=population_summary(native)
        assert abs(mean-row['API_native_mean_kj'])<1e-4
        assert abs(sd-row['API_native_population_SD_kj'])<1e-5
        scale=row['signed_author_scale'];bridge=row['condition_bridge_kj']
        assert scale==original['signed_author_scale'] and scale!=0 and bridge==original['condition_bridge_kj']
        transformed=[v/scale+bridge for v in native]
        mu,sigma=population_summary(transformed)
        delta=max(abs(mu-target['comparator']),abs(sigma-target['spread']));max_summary=max(max_summary,delta)
        assert delta<1e-10
        assert target['reference']==original['reference'] and target['metag']==original['metag']
    rows=inputs['rows'];n=len(rows)
    spread=[r['spread'] for r in rows]
    disagreement=[abs(r['metag']-r['comparator']) for r in rows]
    rank=lambda values:[2*sum(y<x for y in values)+sum(y==x for y in values)+1 for x in values]
    scores=dict(spread=spread,disagreement=disagreement,
        rank_sum=[a+b for a,b in zip(rank(spread),rank(disagreement))],
        prediction_magnitude=[abs(r['comparator']) for r in rows])
    entries=0;aurcs={}
    for name,values in scores.items():
        risk=[]
        for k in range(1,n+1):
            threshold=sorted(values)[k-1]
            certain=[i for i,v in enumerate(values) if v<threshold]
            boundary=[i for i,v in enumerate(values) if v==threshold]
            fraction=Fraction(k-len(certain),len(boundary))
            weights={i:Fraction(1) for i in certain};weights.update({i:fraction for i in boundary})
            expected=report['curves'][name][k-1]
            assert expected['retained']==k and expected['boundary_tier_size']==len(boundary)
            for predictor in ['comparator','metag','zero']:
                error=lambda i:abs((rows[i][predictor] if predictor!='zero' else 0)-rows[i]['reference'])
                mae=math.fsum(float(w)*error(i) for i,w in weights.items())/k
                severe=float(sum(w for i,w in weights.items() if error(i)>20))
                for actual_value,stored in [(mae,expected['MAE'][predictor]),
                    (severe,expected['expected_severe_count'][predictor]),(severe/k,expected['severe_rate'][predictor])]:
                    assert abs(actual_value-stored)<1e-10;entries+=1
                if predictor=='comparator':risk.append(mae)
        aurcs[name]=math.fsum(risk)/n
        assert abs(aurcs[name]-report['comparator_AURC'][name])<1e-10
    assert report['bootstrap']['status']=='not_estimated'
    identity=None
    if identity_replay:
        from rdkit import Chem,RDLogger
        from rdkit.Chem.MolStandardize import rdMolStandardize
        RDLogger.DisableLog('rdApp.*')
        cache={}
        def normalized(smiles):
            if smiles not in cache:
                mol=Chem.AddHs(Chem.MolFromSmiles(smiles))
                for i in range(2):
                    mol=Chem.AddHs(rdMolStandardize.Uncharger().uncharge(rdMolStandardize.Normalize(mol)))
                    if i==0:mol=Chem.AddHs(Chem.MolFromSmiles(Chem.MolToSmiles(Chem.RemoveHs(mol),canonical=True)))
                cache[smiles]=Chem.MolToSmiles(Chem.RemoveHs(mol),canonical=True)
            return cache[smiles]
        def equation(text):
            result=defaultdict(float)
            for sign,side in zip([-1,1],text.split(' = ')):
                for term in side.split(' + '):
                    match=re.fullmatch(r'(?:(\d+(?:\.\d+)?) )?(.+)',term)
                    result[normalized(match[2])]+=sign*float(match[1] or 1)
            return {key:value for key,value in result.items() if value}
        training=list(csv.DictReader((packet/'reference671_training.csv').open()))
        keys=defaultdict(list)
        for index,row in enumerate(training):
            if row['type']=='reaction dg':keys[tuple(sorted(equation(row['reaction']).items()))].append(index)
        for row in raw['raw_reactions']:
            normalized_equation=equation(row['equation'])
            assert normalized_equation==row['normalized_native_equation']
            assert keys.get(tuple(sorted(normalized_equation.items())),[])==row['exact_current671_normalized_equation_rows']
        identity=dict(replayed_equations=69,exact_shipped671_state_equation_matches=sum(bool(r['exact_current671_normalized_equation_rows']) for r in raw['raw_reactions']),
            scope='Exact normalized state equations only. No exact match does not establish molecular, graph, family, label or source independence.')
    return dict(status='passed',authenticated_files=len(manifest['files']),n_reactions=69,n_members=100,
        n_molecular_vectors=len(raw['molecular_vectors']),reaction_member_sums_checked=6900,
        maximum_upstream_float32_sum_difference_kj=max_reconstruction,
        maximum_independent_summary_difference_kj=max_summary,risk_entries_independently_checked=entries,
        native_full_cohort_MAE_kj=report['curves']['spread'][-1]['MAE']['comparator'],
        comparator_AURC=aurcs,identity_replay=identity,
        experimental_validation_regime=raw['experimental_validation_regime'],independent_reliability_advantage_demonstrated=False)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--packet',type=Path,default=Path(__file__).with_name('native_uncertainty_validation'))
    parser.add_argument('--identity-replay',action='store_true')
    args=parser.parse_args()
    print(json.dumps(audit(args.packet,identity_replay=args.identity_replay),indent=2))
