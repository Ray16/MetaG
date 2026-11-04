"""Outcome-free structural feasibility audit, not independent reference admission."""
from __future__ import annotations
import argparse
from collections import Counter,defaultdict
import csv
from fractions import Fraction
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import re
from rdkit import Chem,DataStructs,rdBase
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.MolStandardize import rdMolStandardize
from metag.routing.pka_transform import _neutralize_v2
from metag.routing.sugar_states import enumerate_sugar_states

ALLOWED_COLUMNS=('description','reaction','type','T','I','pH','pMg','EC')
SOURCE_COLUMNS=('part','page','col l/r','table from top','entry nr','reaction',
                'reference','temperature','ionic_strength','p_h','p_mg')
TERM=re.compile(r'^(\d+(?:\.\d+)?(?:/\d+)?)\s+(.+)$')
GEM=Chem.MolFromSmarts('[C;X4]([O;X2;H1;+0])[O;X2;H1;+0]')
TAUTOMERS=rdMolStandardize.TautomerEnumerator()
TAUTOMERS.SetMaxTautomers(64);TAUTOMERS.SetMaxTransforms(128)

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def projected_table(path):
    with Path(path).open(newline='') as stream:
        reader=csv.DictReader(stream)
        if not set(ALLOWED_COLUMNS)<=set(reader.fieldnames or []):raise ValueError('Missing structural/condition columns')
        return [{key:row[key] for key in ALLOWED_COLUMNS} for row in reader]

def equation(text):
    sides=text.split(' = ')
    if len(sides) not in [1,2] or any(not side.strip() for side in sides):raise ValueError('Invalid equation')
    terms=[]
    for side,sign in zip(sides,[-1,1] if len(sides)==2 else [1]):
        for term in side.split(' + '):
            match=TERM.fullmatch(term.strip())
            try:
                coefficient,smiles=(Fraction(match[1]),match[2]) if match else (Fraction(1),term.strip())
            except ZeroDivisionError as error:raise ValueError('Invalid term coefficient') from error
            if coefficient<=0 or not smiles:raise ValueError('Invalid term')
            terms.append((sign*coefficient,smiles))
    return terms

def concrete(smiles):
    mol=Chem.MolFromSmiles(smiles)
    if mol is None or not mol.GetNumAtoms():raise ValueError('Invalid/empty structure')
    if any(a.GetAtomicNum()==0 for a in mol.GetAtoms()):raise ValueError('Abstract structure')
    mol=Chem.RemoveHs(mol)
    for atom in mol.GetAtoms():atom.SetAtomMapNum(0)
    return mol

def reverse_hydration(mol):
    """Conservative exposure alias; does not certify an aqueous equilibrium pool."""
    mol=Chem.Mol(mol)
    for _ in range(mol.GetNumAtoms()):
        hits=mol.GetSubstructMatches(GEM)
        if not hits:return mol
        carbon,keep,remove=min(hits);rw=Chem.RWMol(mol)
        rw.GetAtomWithIdx(carbon).SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
        oxygen=rw.GetAtomWithIdx(keep);oxygen.SetNumExplicitHs(0);oxygen.SetNoImplicit(True)
        rw.GetBondBetweenAtoms(carbon,keep).SetBondType(Chem.BondType.DOUBLE)
        rw.RemoveAtom(remove);mol=rw.GetMol();Chem.SanitizeMol(mol)
    raise ValueError('Hydration alias did not terminate')

@lru_cache(maxsize=None)
def family(smiles):
    mol=concrete(smiles);exact=Chem.MolToSmiles(mol,isomericSmiles=True)
    if exact in ['O','[H+]']:return dict(exact=exact,pool=exact,broad=exact,reservoir=True,issues=[])
    if any(a.GetIsotope() for a in mol.GetAtoms()):return dict(exact=exact,pool=None,broad=None,reservoir=False,issues=['isotope identity requires separate review'])
    if any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()) or len(Chem.GetMolFrags(mol))!=1:
        return dict(exact=exact,pool=None,broad=None,reservoir=False,issues=['radical or disconnected identity requires separate review'])
    issues=[]
    try:
        ref,*unused=_neutralize_v2(exact)
        if ref is None:raise ValueError('No global proton reference')
        mol=reverse_hydration(concrete(ref))
        state=enumerate_sugar_states(Chem.MolToSmiles(mol,isomericSmiles=True))
        if state.applicable and state.complete and state.open_smiles:
            mol=concrete(state.open_smiles)
        elif state.applicable and not state.complete:
            issues.append('pool alias unresolved: '+str(state.reason))
        pool=Chem.MolToSmiles(mol,isomericSmiles=True) if not issues else None
        flat=Chem.Mol(mol);Chem.RemoveStereochemistry(flat)
        enumerated=TAUTOMERS.Enumerate(flat)
        if issues:
            # A failed sugar alias must not become a purportedly absent family
            # merely because tautomer enumeration of its starting drawing ended.
            broad=None
        elif enumerated.status!=rdMolStandardize.TautomerEnumeratorStatus.Completed:
            broad=None;issues.append('tautomer enumeration incomplete')
        else:broad=Chem.MolToSmiles(TAUTOMERS.PickCanonical(enumerated),isomericSmiles=False)
    except (ValueError,RuntimeError) as error:
        pool=None;broad=None;issues.append('reference alias unresolved: '+str(error))
    return dict(exact=exact,pool=pool,broad=broad,reservoir=False,issues=issues)

def reaction_key(parts,tier):
    counts=defaultdict(Fraction);participants=[]
    for coefficient,smiles in parts:
        record=family(smiles)
        if record['reservoir']:continue
        key=record[tier]
        if key is None:return None
        coefficient=Fraction(str(coefficient));counts[key]+=coefficient;participants.append(key)
    values=sorted((key,value) for key,value in counts.items() if value)
    if not values:
        # Different sugar pools must not all collide into an empty signature.
        return ('within-family',tuple(sorted(set(participants))))
    denominator=1
    import math
    denominator=math.lcm(*(v.denominator for _,v in values))
    integers=[(key,int(value*denominator)) for key,value in values]
    divisor=math.gcd(*(abs(value) for _,value in integers))
    normal=tuple((key,value//divisor) for key,value in integers)
    return min(normal,tuple((key,-value) for key,value in normal))

def development_panel(path):
    data=json.loads(Path(path).read_text());reactions={};molecules={tier:set() for tier in ['exact','pool','broad']};unresolved=[]
    for rid,row in data.items():
        parts=[];complete=True
        for name,item in row.get('species',{}).items():
            try:
                coefficient=Fraction(str(item[0]))
                if not coefficient:raise ValueError('Zero coefficient')
                parts.append((coefficient,item[2]))
                info=family(item[2])
                if info['issues']:unresolved.append(dict(reaction=rid,species=name,issues=info['issues']))
                if not info['reservoir']:
                    for tier in molecules:
                        if info[tier] is not None:molecules[tier].add(info[tier])
            except (ValueError,RuntimeError,IndexError,TypeError) as error:
                complete=False;unresolved.append(dict(reaction=rid,species=name,issues=[str(error)]))
        # Retain concrete partners even when a different term is abstract.
        try:
            if not complete or not parts:raise ValueError('Incomplete development reaction')
            reactions[rid]={tier:reaction_key(parts,tier) for tier in ['exact','pool']}
        except (ValueError,RuntimeError):reactions[rid]={'exact':None,'pool':None}
    return reactions,molecules,unresolved

def audit(training,panel):
    rows=projected_table(training);development,molecules,development_unresolved=development_panel(panel)
    indexes={tier:defaultdict(list) for tier in ['exact','pool']}
    for rid,keys in development.items():
        for tier,key in keys.items():
            if key is not None:indexes[tier][key].append(rid)
    output=[];counts=Counter()
    for index,row in enumerate(rows):
        if row['type']!='reaction dg':continue
        counts['reaction_rows']+=1
        record=dict(row=index,description=row['description'],condition_metadata={k:row[k] for k in ['T','I','pH','pMg']},experimental_reference_admitted=False)
        try:
            if ' = ' not in row['reaction']:raise ValueError('Reaction target lacks two sides')
            parts=equation(row['reaction']);aliases=[dict(smiles=s,**family(s)) for _,s in parts]
            matches={tier:indexes[tier].get(reaction_key(parts,tier),[]) for tier in ['exact','pool']}
            absent={tier:sorted({a[tier] for a in aliases if not a['reservoir'] and a[tier] is not None and a[tier] not in molecules[tier]}) for tier in molecules}
            unresolved=[dict(smiles=a['smiles'],issues=a['issues']) for a in aliases if a['issues']]
            record.update(status='structural_candidate_only',matches_development=matches,absent_from_development_families=absent,unresolved_aliases=unresolved,
                equation_smiles=row['reaction'],molecules=aliases)
            counts['concrete_reaction_rows']+=1
            for tier in ['exact','pool']:
                if matches[tier]:counts[tier+'_equation_match']+=1
            if matches['exact'] or matches['pool']:counts['any_operational_equation_match']+=1
            else:
                counts['no_operational_panel_equation_match']+=1
                if absent['broad']:counts['nonmatched_with_absent_broad_family']+=1
                if unresolved:counts['nonmatched_with_unresolved_alias']+=1
        except (ValueError,RuntimeError) as error:
            counts['unresolved_reaction_rows']+=1;record.update(status='unresolved_structure',reason=str(error))
        output.append(record)
    return dict(schema='metag-benchmark-feasibility-v2',scope='Structural feasibility relative to the supplied development panel only; no historical-independence, source-ancestry, condition, novelty or accuracy certification.',
        input_sha256=dict(training=sha(training),development_panel=sha(panel)),source_sha256=sha(__file__),rdkit_version=rdBase.rdkitVersion,
        projected_training_columns=list(ALLOWED_COLUMNS),outcome_values_exported=False,predictions_performed=False,experimental_references_admitted=0,
        development_reactions=len(development),development_family_counts={k:len(v) for k,v in molecules.items()},
        development_unresolved_terms=development_unresolved,independence_certified=False,counts=dict(counts),rows=output)


def holdout_plan(training, evaluation_rows, heldout_smiles):
    """Propose structural exclusions across reaction AND formation targets.

    This does not fit a model or certify source independence. Every unresolved
    training row is excluded conservatively, including abstract wildcard labels.
    A future scaffold holdout needs its own declared grouping; these are broad
    molecular-family holdouts, not scaffold or functional-group holdouts.
    """
    rows=projected_table(training)
    evaluation_rows=set(evaluation_rows)
    if (not evaluation_rows or any(type(i) is not int or not 0<=i<len(rows) for i in evaluation_rows)
            or any(rows[i]['type']!='reaction dg' for i in evaluation_rows)):
        raise ValueError('Evaluation rows must be valid reaction target indices')
    heldout={}
    for smiles in heldout_smiles:
        info=family(smiles)
        if info['reservoir'] or info['broad'] is None or info['pool'] is None:
            raise ValueError('Held-out family must be resolved and nonreservoir')
        heldout[smiles]=info
    if not heldout:raise ValueError('At least one held-out molecular family is required')
    families={info['broad'] for info in heldout.values()}
    evaluation_keys=set();evaluation_families=set()
    for index in evaluation_rows:
        parts=equation(rows[index]['reaction'])
        if ' = ' not in rows[index]['reaction']:raise ValueError('Reaction target lacks two sides')
        for tier in ('exact','pool'):
            key=reaction_key(parts,tier)
            if key is not None:evaluation_keys.add((tier,key))
        evaluation_families.update(family(s)['broad'] for _,s in parts)
    if not families<=evaluation_families:
        raise ValueError('Held-out family does not participate in the evaluation set')
    retained=[];excluded=[];retained_molecules=set();original_other_pools=set()
    for index,row in enumerate(rows):
        reasons=[]
        if index in evaluation_rows:reasons.append('evaluation target')
        infos=[]
        try:
            if row['type'] not in ('reaction dg','formation dg'):
                raise ValueError('Unknown training target type')
            parts=equation(row['reaction'])
            if (row['type']=='reaction dg')!=(' = ' in row['reaction']):
                raise ValueError('Training target type disagrees with equation')
            # Keep concrete partners from every row for the descriptive original
            # domain comparison, even when another term later excludes the row.
            for _,smiles in parts:
                try:
                    info=family(smiles)
                    if (not info['reservoir'] and info['pool'] is not None
                            and info['broad'] not in families
                            and info['pool'] not in {h['pool'] for h in heldout.values()}):
                        original_other_pools.add(info['pool'])
                except (ValueError,RuntimeError):pass
            for _,smiles in parts:
                info=family(smiles);infos.append(info)
                if info['broad'] is None or info['pool'] is None:
                    raise ValueError('Unresolved training family')
            if any(info['broad'] in families for info in infos):
                reasons.append('held-out molecular family')
            if row['type']=='reaction dg' and any(
                    (tier,reaction_key(parts,tier)) in evaluation_keys for tier in ('exact','pool')):
                reasons.append('exact or pool-equivalent evaluation equation')
        except (ValueError,RuntimeError) as error:
            reasons.append('unresolved training row: '+str(error))
        if reasons:excluded.append(dict(row=index,target_type=row['type'],reasons=reasons))
        else:
            retained.append(index)
            retained_molecules.update(info['pool'] for info in infos if not info['reservoir'])
    if not retained or not retained_molecules:raise ValueError('No resolved training domain remains')
    generator=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048,includeChirality=True)
    training_keys=sorted(retained_molecules)
    fingerprints=[generator.GetFingerprint(concrete(s)) for s in training_keys]
    original_keys=sorted(original_other_pools)
    original_fingerprints=[generator.GetFingerprint(concrete(s)) for s in original_keys]
    neighbours=[]
    for smiles,info in sorted(heldout.items()):
        similarities=DataStructs.BulkTanimotoSimilarity(generator.GetFingerprint(concrete(info['pool'])),fingerprints)
        best=max(range(len(training_keys)),key=lambda i:similarities[i])
        original_similarities=DataStructs.BulkTanimotoSimilarity(generator.GetFingerprint(concrete(info['pool'])),original_fingerprints)
        original_best=max(range(len(original_keys)),key=lambda i:original_similarities[i])
        neighbours.append(dict(supplied_smiles=smiles,broad_family=info['broad'],pool_key=info['pool'],
            nearest_retained_pool=training_keys[best],nearest_similarity=similarities[best],
            below_declared_0_5_threshold=similarities[best]<0.5,
            original_other_family_nearest_pool=original_keys[original_best],
            original_other_family_nearest_similarity=original_similarities[original_best],
            original_domain_comparison_scope='Concrete pool keys across all original rows, excluding resolved held-out broad families. Incomplete rows may contain additional unknown chemistry; this is not a trained comparator.'))
    return dict(schema='metag-molecular-holdout-plan-v1',training_sha256=sha(training),
        source_sha256=sha(__file__),rdkit_version=rdBase.rdkitVersion,
        evaluation_rows=sorted(evaluation_rows),heldout_broad_families=sorted(families),
        retained_rows=retained,excluded_rows=excluded,
        counts=dict(training_rows=len(rows),retained_rows=len(retained),excluded_rows=len(excluded),
                    retained_by_target_type=dict(Counter(rows[i]['type'] for i in retained))),
        fingerprint=dict(kind='Morgan',radius=2,bits=2048,chirality=True),heldout_neighbours=neighbours,
        any_heldout_family_below_0_5=any(n['below_declared_0_5_threshold'] for n in neighbours),
        unresolved_training_row_policy='Exclude entire row, including every abstract/wildcard row; retain the exclusion reason.',
        scope='Molecular-family holdout feasibility only. Source ancestry, measurement duplicates and historical development exposure remain unresolved; this is not a scaffold holdout or independent reference enrollment.',
        launch_allowed=False,model_training_performed=False,predictions_performed=False,experimental_references_admitted=0)


def named_equation_key(text):
    """Strict name matching for source candidates; retain phase and stereo text.

    Only letter case, whitespace, primitive stoichiometric extent and direction
    are normalized. No synonym, molecular identity or ancestry is inferred.
    """
    if text.count(' = ')!=1:raise ValueError('Named equation must have two sides')
    counts=defaultdict(Fraction)
    for coefficient,name in equation(text):
        counts[' '.join(name.casefold().split())]+=coefficient
    values=sorted((name,c) for name,c in counts.items() if c)
    if not values:raise ValueError('Empty named equation')
    import math
    denominator=math.lcm(*(v.denominator for _,v in values))
    integers=[(name,int(value*denominator)) for name,value in values]
    divisor=math.gcd(*(abs(value) for _,value in integers))
    normal=tuple((name,value//divisor) for name,value in integers)
    return min(normal,tuple((name,-value) for name,value in normal))


def source_candidates(feasibility_report, source_csv):
    """Link candidate equations to source metadata without exporting outcomes."""
    index=defaultdict(list);unresolved=[]
    with Path(source_csv).open(newline='') as stream:
        reader=csv.DictReader(stream)
        if not set(SOURCE_COLUMNS)<=set(reader.fieldnames or []):
            raise ValueError('Missing OpenTECR metadata columns')
        for number,raw in enumerate(reader):
            row={key:raw[key] for key in SOURCE_COLUMNS}
            try:index[named_equation_key(row['reaction'])].append(dict(source_row=number,**row))
            except ValueError as error:unresolved.append(dict(source_row=number,reason=str(error)))
    output=[];counts=Counter()
    for candidate in feasibility_report['rows']:
        if (candidate['status']!='structural_candidate_only'
                or any(candidate['matches_development'].values())):continue
        counts['unmatched_candidate_equations']+=1
        try:hits=index.get(named_equation_key(candidate['description']),[])
        except ValueError:hits=[]
        counts['with_named_source_candidates']+=bool(hits)
        # These are metadata-presence counts. A filled field may still have an
        # invalid value or a different thermodynamic convention. Never infer
        # Mg-free conditions or zero ionic strength from an empty source field.
        complete=[h for h in hits if all(h[k].strip() for k in (
            'temperature','ionic_strength','p_h','p_mg'))]
        counts['with_all_four_condition_fields_present']+=bool(complete)
        output.append(dict(training_row=candidate['row'],description=candidate['description'],
            named_source_candidates=hits,source_reference_codes=sorted({h['reference'] for h in hits}),
            all_four_condition_fields_present_in_any_row=bool(complete),
            ancestry_confirmed=False,experimental_reference_admitted=False))
    return dict(schema='metag-named-source-candidates-v1',source_csv_sha256=sha(source_csv),
        training_sha256=feasibility_report['input_sha256']['training'],source_sha256=sha(__file__),
        matching_policy=named_equation_key.__doc__,source_projection=list(SOURCE_COLUMNS),
        source_rows_with_unresolved_named_equations=unresolved,counts=dict(counts),rows=output,
        scope='Named-equation source leads only; no assertion that the training target derives from any particular matching measurement. Native condition validity, bridges, primary texts and historical outcome use remain to be adjudicated.',
        outcome_values_exported=False,experimental_references_admitted=0,predictions_performed=False)


def holdout_previews(training, feasibility_report):
    """Preview every fully resolved unmatched equation with an absent family."""
    plans=[]
    for row in feasibility_report['rows']:
        if (row['status']!='structural_candidate_only' or any(row['matches_development'].values())
                or row['unresolved_aliases'] or not row['absent_from_development_families']['broad']):
            continue
        absent=set(row['absent_from_development_families']['broad'])
        heldout=sorted({m['smiles'] for m in row['molecules'] if m['broad'] in absent})
        try:
            plan=holdout_plan(training,[row['row']],heldout)
            plan['description']=row['description'];plans.append(plan)
        except ValueError as error:
            plans.append(dict(evaluation_rows=[row['row']],status='unresolved_holdout_plan',reason=str(error)))
    return dict(schema='metag-holdout-previews-v1',n_plans=len(plans),plans=plans,
        n_meeting_retained_domain_molecular_distance_screen=sum(p.get('any_heldout_family_below_0_5',False) for p in plans),
        scope='Molecular holdout feasibility under conservative training-row exclusions. Reduced training support must be controlled in any comparative experiment. These are not enrolled independent test cases.',
        launch_allowed=False,experimental_references_admitted=0)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--training',type=Path,required=True)
    p.add_argument('--panel',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--source-metadata',type=Path)
    p.add_argument('--source-output',type=Path)
    p.add_argument('--holdout-output',type=Path)
    args=p.parse_args()
    if bool(args.source_metadata)!=bool(args.source_output):
        p.error('--source-metadata and --source-output must be supplied together')
    report=audit(args.training,args.panel)
    args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    if args.source_metadata:
        args.source_output.write_text(json.dumps(source_candidates(report,args.source_metadata),indent=2,allow_nan=False)+'\n')
    if args.holdout_output:
        args.holdout_output.write_text(json.dumps(holdout_previews(args.training,report),indent=2,allow_nan=False)+'\n')
    print(json.dumps(report['counts'],indent=2))
if __name__=='__main__':main()
