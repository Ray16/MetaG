"""Independently reconstruct explicit cofactor curation accounting offline.

The default audit checks owned primary structures, derived atom accounting and
recorded input gates. --modelseed additionally authenticates upstream tables
and recomputes their graph/formula records. Neither mode computes free energies.
"""
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
import argparse
import csv
import gzip
import hashlib
import json
import importlib.util
import os
import re
import subprocess
import sys
import tempfile
import zipfile

from rdkit import Chem, rdBase
from rdkit.Chem import rdMolDescriptors

read = lambda p: json.loads(Path(p).read_text())
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()


def rows(path, key):
    values={}
    with gzip.open(path,'rt') as stream:
        for line in stream:
            row=json.loads(line)
            assert row[key] not in values
            values[row[key]]=row
    return values


def manifest(packet):
    expected=read(packet/'manifest.json')['files']
    actual={str(p.relative_to(packet)) for p in packet.rglob('*') if p.is_file() and '__pycache__' not in p.parts}
    assert actual==set(expected)|{'manifest.json'}
    for name,record in expected.items():
        assert not Path(name).is_absolute() and '..' not in Path(name).parts
        path=packet/name
        assert not path.is_symlink() and path.stat().st_size==record['bytes'] and sha(path)==record['sha256']


def graph(smiles):
    with rdBase.BlockLogs():mol=Chem.MolFromSmiles(smiles or '')
    if mol is None or not mol.GetNumAtoms() or any(a.GetAtomicNum()==0 for a in mol.GetAtoms()):return None
    formula=Counter((str(a.GetIsotope()) if a.GetIsotope() else '')+a.GetSymbol()
                    for a in Chem.AddHs(mol).GetAtoms())
    canonical=Chem.MolToSmiles(mol,isomericSmiles=True)
    charge=Chem.GetFormalCharge(mol)
    identity=hashlib.sha256(json.dumps([charge,canonical],sort_keys=True).encode()).hexdigest()
    return dict(formula=dict(formula),formal_charge=charge,canonical_identity_sha256=identity)


def column(values):
    result=defaultdict(Fraction)
    for cid,comp,nu in values:result[(cid,comp)]+=Fraction(nu)
    return {key:nu for key,nu in result.items() if nu}


def residual(terms,accounting):
    elements=defaultdict(Fraction);formal=declared=Fraction(0);known=q_known=True;signature=defaultdict(Fraction)
    for (cid,comp),nu in terms.items():
        rec=accounting.get(cid)
        if rec is None:known=q_known=False;continue
        if rec['declared_charge'] is None:q_known=False
        else:declared+=nu*rec['declared_charge']
        state=rec['graph']
        if state is None:known=False;continue
        for element,count in state['formula'].items():elements[element]+=nu*count
        formal+=nu*state['formal_charge']
        signature['explicit_proton' if cid=='cpd00067' else state['canonical_identity_sha256']]+=nu
    return dict(all_formulas_known=known,element_residual={k:str(v) for k,v in sorted(elements.items()) if v} if known else None,
        formal_charge_residual=str(formal) if known else None,declared_charge_residual=str(declared) if q_known else None,
        graph_net_empty=not any(signature.values()) if known else None)


def audit(packet, modelseed=None):
    packet=Path(packet);manifest(packet)
    summary=read(packet/'summary.json');registry=read(packet/'registry.json');reviews=read(packet/'identity_reviews.json')
    original=read(packet/'original_compound_rows.json')
    states=rows(packet/'literal_compound_admission.jsonl.gz','compound')
    literal=rows(packet/'literal_reaction_admission.jsonl.gz','reaction')
    universe=rows(packet/'source_universe.jsonl.gz','reaction')
    affected=rows(packet/'affected_reactions.jsonl.gz','reaction')
    recorded=rows(packet/'curated_outcomes.jsonl.gz','reaction')
    accounting=rows(packet/'atom_accounting.jsonl.gz','compound')
    curated_accounting=dict(accounting)
    target_states=read(packet/'curated_compound_states.json')
    families=read(packet/'anchor_families.json')
    assert set(registry)==set(target_states)==set(reviews)==set(original)=={'cpd00003','cpd00004','cpd00005','cpd00006'}
    assert set(literal)==set(recorded)=={rid for rid,row in universe.items() if row['is_obsolete']==row['is_transport']=='0'}
    assert len(literal)==summary['n_requested']==43650
    stereo_counts={};verified_primary=0
    for cid,entry in registry.items():
        target=entry['target'];source=entry['source'];review=reviews[cid]
        assert all(original[cid][k]==source[k] for k in source)
        molpath=packet/'primary_sources'/(entry['chebi_id'].split(':')[1]+'.mol')
        assert sha(molpath)==entry['source_evidence']['molfile_sha256']
        mol=Chem.MolFromMolFile(str(molpath));assert mol is not None
        smi=Chem.MolToSmiles(mol,isomericSmiles=True);q=Chem.GetFormalCharge(mol)
        assert (smi,q,Chem.MolToInchiKey(mol))==(target['smiles'],target['charge'],target['inchikey'])
        pubpath=packet/'primary_sources'/(target['inchikey']+'.json')
        assert sha(pubpath)==entry['source_evidence']['PubChem_properties_sha256']
        pub=read(pubpath)['PropertyTable']['Properties'];assert len(pub)==1;pub=pub[0]
        assert pub['CID']==entry['source_evidence']['PubChem_CID'] and pub['Charge']==q
        assert pub['InChIKey']==target['inchikey'] and pub['MolecularFormula']==rdMolDescriptors.CalcMolFormula(mol)
        assert Chem.MolToSmiles(Chem.MolFromSmiles(pub['SMILES']),isomericSmiles=True)==smi
        oldmol=Chem.MolFromSmiles(source['smiles']);assert oldmol is not None
        # The recorded mapping uses parsed source-literal/target-SMILES order,
        # independently of the Molfile's original atom numbering.
        mol=Chem.MolFromSmiles(smi)
        Chem.AssignStereochemistry(oldmol,cleanIt=True,force=True)
        Chem.AssignStereochemistry(mol,cleanIt=True,force=True)
        mapping={int(k):v for k,v in review['atom_mapping'].items()}
        assert set(mapping)==set(range(oldmol.GetNumAtoms()))
        assert set(mapping.values())==set(range(mol.GetNumAtoms()))
        count=0
        for i,j in mapping.items():
            a,b=oldmol.GetAtomWithIdx(i),mol.GetAtomWithIdx(j)
            assert (a.GetAtomicNum(),a.GetIsotope())==(b.GetAtomicNum(),b.GetIsotope())
            if a.GetAtomicNum()==6:
                ca=a.GetProp('_CIPCode') if a.HasProp('_CIPCode') else None
                cb=b.GetProp('_CIPCode') if b.HasProp('_CIPCode') else None
                assert ca==cb
                count+=ca is not None
        assert count==review['declared_carbon_stereocenters_preserved']==8
        for bond in oldmol.GetBonds():
            assert mol.GetBondBetweenAtoms(mapping[bond.GetBeginAtomIdx()],mapping[bond.GetEndAtomIdx()]) is not None
        oldgraph,newgraph=graph(source['smiles']),graph(smi)
        keys=set(oldgraph['formula'])|set(newgraph['formula'])
        differences={k:newgraph['formula'].get(k,0)-oldgraph['formula'].get(k,0) for k in keys}
        assert all(v==0 for k,v in differences.items() if k!='H')
        dh=differences.get('H',0)
        assert dh==entry['delta_H']==newgraph['formal_charge']-oldgraph['formal_charge']
        assert oldgraph==accounting[cid]['graph']
        assert not review['energy_equivalence_established']
        curated_accounting[cid]=dict(compound=cid,declared_charge=q,graph=newgraph)
        assert target_states[cid]['status']=='passes_cpu_state_gates'
        assert target_states[cid]['canonical_input']==smi
        stereo_counts[cid]=count;verified_primary+=1
    combined_states=dict(states,**target_states)
    affected_expected={rid for rid,row in literal.items() if set(registry)&set(row['compounds'])}
    assert set(affected)==affected_expected
    outcomes=dict(literal);n_known=n_proton=0
    boundary=[];passes_affected=0
    for rid,row in affected.items():
        before_terms,after_terms=column(row['original_terms']),column(row['curated_terms'])
        expected=defaultdict(Fraction,before_terms);changes=defaultdict(Fraction)
        for (cid,comp),nu in before_terms.items():
            if cid in registry:changes[comp]-=nu*registry[cid]['delta_H']
        for comp,nu in changes.items():expected[('cpd00067',comp)]+=nu
        expected={key:nu for key,nu in expected.items() if nu}
        assert after_terms==expected
        nonzero={comp:str(nu) for comp,nu in sorted(changes.items()) if nu}
        assert nonzero==row['proton_adjustments_by_compartment'];n_proton+=bool(nonzero)
        before,after=residual(before_terms,accounting),residual(after_terms,curated_accounting)
        assert before==row['balance_before'] and after==row['balance_after']
        assert before['all_formulas_known']==after['all_formulas_known']
        assert before['declared_charge_residual']==after['declared_charge_residual']
        if before['all_formulas_known']:
            assert before['element_residual']==after['element_residual'] and before['formal_charge_residual']==after['formal_charge_residual']
            n_known+=1
        old=literal[rid];assert row['literal_outcome']==old
        new=row['curated_outcome'];assert new==recorded[rid]
        assert new['source_status']==old['source_status']
        assert new['compounds']==sorted({cid for cid,comp in after_terms})
        causes=[]
        # Non-state declarations are retained; all affected columns parsed.
        if any(c['category']=='abstract_reaction' for c in old['causes']):causes.append(dict(category='abstract_reaction'))
        if len({comp for cid,comp in after_terms})>1:causes.append(dict(category='multiple_compartments_need_conditions'))
        for cid,comp in sorted((cid,comp) for cid,comp in after_terms if cid!='cpd00067'):
            if cid not in combined_states:causes.append(dict(category='missing_compound',compound=cid))
            elif combined_states[cid]['status']!='passes_cpu_state_gates':
                causes.append(dict(category=combined_states[cid]['category'],compound=cid))
        if not any('compound' in c for c in causes) and after_terms:
            assert after['all_formulas_known']
            if after['element_residual']:causes.append(dict(category='element_imbalance',residual=after['element_residual']))
            if Fraction(after['formal_charge_residual']):causes.append(dict(category='charge_imbalance',residual=after['formal_charge_residual']))
            if after['graph_net_empty']:causes.append(dict(category='identity_or_graph_collapsed_reaction'))
        expected_causes=sorted(json.dumps(c,sort_keys=True) for c in causes)
        assert expected_causes==sorted(json.dumps(c,sort_keys=True) for c in new['causes']),rid
        assert (not causes)==(new['status']=='passes_cpu_input_gates')
        coefficients={name:str(sum(nu for (cid,comp),nu in after_terms.items() if cid in {f['oxidized'],f['reduced']}))
                      for name,f in families.items()}
        assert coefficients==row['anchor_family_net_coefficients']
        if not causes:
            passes_affected+=1
            if any(Fraction(v) for v in coefficients.values()):boundary.append(rid)
        outcomes[rid]=new
    assert outcomes==recorded
    gained=sorted(rid for rid in outcomes if literal[rid]['status']!='passes_cpu_input_gates' and outcomes[rid]['status']=='passes_cpu_input_gates')
    lost=sorted(rid for rid in outcomes if literal[rid]['status']=='passes_cpu_input_gates' and outcomes[rid]['status']!='passes_cpu_input_gates')
    assert gained==summary['gained_inputs'] and lost==summary['lost_inputs']
    assert dict(Counter(r['status'] for r in literal.values()))==summary['literal_status_counts']
    assert dict(Counter(r['status'] for r in outcomes.values()))==summary['curated_status_counts']
    assert dict(Counter(r['causes'][0]['category'] for r in outcomes.values() if r['causes']))==summary['first_remaining_failure']
    warnings=Counter()
    for row in outcomes.values():
        if row['status']=='passes_cpu_input_gates':
            issues={issue for cid in row['compounds'] for issue in combined_states[cid].get('protonation_issues',[])}
            warnings['with_protonation_warnings' if issues else 'without_recorded_protonation_warnings']+=1
    assert dict(warnings)==summary['passing_warning_counts']
    assert (len(affected),len(literal)-len(affected),n_known,n_proton,passes_affected)==(
        summary['n_affected'],summary['n_unaffected'],summary['n_known_formula_residuals_preserved'],
        summary['n_proton_adjusted_reactions'],summary['affected_passing'])
    assert sorted(boundary)==summary['anchor_retaining_input_ids']
    assert summary['physical_calculations_attempted']==0 and summary['numerical_prediction_coverage'] is None
    assert not summary['production_policy_adopted'] and not summary['experimental_accuracy_assessed']
    raw_checked=0
    if modelseed is not None:
        source=read(packet/'upstream_provenance.json');tables={};reaction_tables={}
        for name,digest in source['source_files_sha256'].items():
            path=Path(modelseed)/name;assert sha(path)==digest
            if re.fullmatch(r'compound_\d+\.tsv',path.name):
                with path.open() as stream:
                    for row in csv.DictReader(stream,delimiter='\t'):tables[row['id']]=row
            elif re.fullmatch(r'reaction_\d+\.tsv',path.name):
                with path.open() as stream:
                    for row in csv.DictReader(stream,delimiter='\t'):reaction_tables[row['id']]=row
        for cid,owned in accounting.items():
            row=tables[cid]
            assert hashlib.sha256(row['smiles'].encode()).hexdigest()==owned['smiles_sha256']
            assert graph(row['smiles'])==owned['graph']
            try:q=Fraction(row['charge']);q=int(q) if q.denominator==1 else None
            except (ValueError,ZeroDivisionError):q=None
            assert q==owned['declared_charge'];raw_checked+=1
        for rid,row in affected.items():
            text=reaction_tables[rid]['stoichiometry']
            assert hashlib.sha256(text.encode()).hexdigest()==row['source_stoichiometry_sha256']
            parsed=defaultdict(Fraction)
            for token in text.split(';'):
                if not token.strip():continue
                coefficient,cid,comp,*_=token.split(':')
                parsed[(cid,comp)]+=Fraction(coefficient)
            assert {k:v for k,v in parsed.items() if v}==column(row['original_terms'])
    return dict(status='passed',primary_identities_verified=verified_primary,carbon_centers_preserved=stereo_counts,
        n_requested=43650,affected_inputs=len(affected),known_residuals_preserved=n_known,
        proton_adjusted_inputs=n_proton,newly_admitted_inputs=len(gained),lost_inputs=len(lost),
        curated_input_passes=summary['curated_status_counts']['passes_cpu_input_gates'],
        anchor_retaining_inputs=len(boundary),raw_upstream_compound_graphs_rechecked=raw_checked,
        physical_predictions_added=0,scope='Primary molecular identity plus input accounting; no Gibbs-energy accuracy or numerical coverage established.')


def cold_targets(packet):
    """Recheck all four target plans using the owned, frozen runtime source."""
    packet=Path(packet).resolve()
    with tempfile.TemporaryDirectory(prefix='metag-cofactor-targets-') as temp:
        root=Path(temp)
        freeze=read(packet/'counterfactual_freeze.json')
        with zipfile.ZipFile(packet/'runtime_source.zip') as archive:
            entries=archive.infolist()
            assert len({e.filename for e in entries})==len(entries)
            assert sum(e.file_size for e in entries)<10_000_000
            assert {e.filename for e in entries}==set(freeze['code_sha256'])
            for entry in entries:
                name=Path(entry.filename)
                assert not name.is_absolute() and '..' not in name.parts
                assert (entry.external_attr>>16)&0o170000!=0o120000
                data=archive.read(entry)
                assert hashlib.sha256(data).hexdigest()==freeze['code_sha256'][entry.filename]
                target=root/name;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(data)
        script=root/'verify_targets.py'
        script.write_text('''import importlib.abc,json,pathlib,socket,sys
root=pathlib.Path(__file__).parent
sys.path.insert(0,str(root/'src'))
class BlockQuantum(importlib.abc.MetaPathFinder):
 def find_spec(self,fullname,path=None,target=None):
  if fullname.split('.')[0] in {'torch','fairchem','pyscf','gpu4pyscf','cupy'}: raise RuntimeError('Quantum import forbidden: '+fullname)
sys.meta_path.insert(0,BlockQuantum())
def deny(*a,**k): raise RuntimeError('Network forbidden')
socket.socket.connect=deny
from metag.state_function import plan_state
from metag.energetics.conformers import spin_multiplicity
packet=pathlib.Path(sys.argv[1]);registry=json.loads((packet/'registry.json').read_text());review=json.loads((packet/'identity_reviews.json').read_text())
for cid,row in registry.items():
 target=row['target'];spin_multiplicity(target['smiles'],target['charge'])
 assert plan_state(target['charge'],target['smiles'])==review[cid]['target_plan'],cid
print(json.dumps(dict(status='passed',targets=4,original_stereo_gate_used=True,quantum_imports_forbidden=True,network_forbidden=True)))
''')
        env=dict(os.environ,OMP_NUM_THREADS='1',OMP_THREAD_LIMIT='1',MKL_NUM_THREADS='1',
            OPENBLAS_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
        result=subprocess.run([sys.executable,str(script),str(packet)],cwd=root,env=env,
            capture_output=True,text=True,timeout=60)
        if result.returncode:raise RuntimeError(result.stdout+result.stderr)
        return json.loads(result.stdout.strip().splitlines()[-1])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--packet',type=Path,default=Path(__file__).with_name('modelseed_cofactor_curation_validation'))
    parser.add_argument('--modelseed',type=Path)
    parser.add_argument('--cold-targets',action='store_true')
    args=parser.parse_args()
    result=audit(args.packet,args.modelseed)
    if args.cold_targets:result['cold_targets']=cold_targets(args.packet)
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
