"""Shared-state checks, without reading experimental errors or fitting values."""
from pathlib import Path
import copy
import hashlib
import itertools
import json
import math
import sys

HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE/'assembly_source'),str(HERE/'source_helpers')]
from rdkit import Chem
from metag.constants import G_HPLUS
from candidate import GlobalAminoEstimator
from run_reassembly import OwnedBackend,MissingRecord


def proton_drawings(reference):
    mol=Chem.MolFromSmiles(reference)
    acid_oxygen=[a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum()==8 and a.GetTotalNumHs()==1
        and any(n.GetAtomicNum()==6 and any(b.GetBondType()==Chem.BondType.DOUBLE and b.GetOtherAtom(n).GetAtomicNum()==8
             for b in n.GetBonds()) for n in a.GetNeighbors())]
    nitrogen=next(a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum()==7)
    for flags in itertools.product([0,1],repeat=len(acid_oxygen)+1):
        altered=Chem.RWMol(mol)
        for oxygen,flag in zip(acid_oxygen,flags):
            if flag:
                atom=altered.GetAtomWithIdx(oxygen);atom.SetNumExplicitHs(0)
                atom.SetNoImplicit(True);atom.SetFormalCharge(-1)
        if flags[-1]:
            atom=altered.GetAtomWithIdx(nitrogen);atom.SetNumExplicitHs(2);atom.SetFormalCharge(0)
        changed=altered.GetMol();Chem.SanitizeMol(changed)
        yield Chem.GetFormalCharge(changed),Chem.MolToSmiles(changed,isomericSmiles=True)


def main():
    catalog_path=Path(sys.argv[1]) if len(sys.argv)>1 else HERE/'catalog.json'
    result_root=Path(sys.argv[2]) if len(sys.argv)>2 else HERE/'reassembly'
    output=Path(sys.argv[3]) if len(sys.argv)>3 else HERE/'invariances594.json'
    catalog=json.loads(catalog_path.read_text());models=json.loads((HERE/'models_nominal.json').read_text())
    
    if 'source_profiles' in catalog:
        from augment_and_reassemble import ExplicitProfileBackend
        backend=ExplicitProfileBackend(catalog,catalog_path)
    else:backend=OwnedBackend(catalog)
    backend.producer['raw_catalog_sha256']=hashlib.sha256(catalog_path.read_bytes()).hexdigest()
    estimator=GlobalAminoEstimator(backend,models)
    original=json.loads((HERE/'inputs.json').read_text())
    allowed={'species','n_Hplus','species_state_pools','species_protonation_models','pka_sites'}
    candidate_file=result_root/(sys.argv[4] if len(sys.argv)>4 else 'global_amino_extended594.json')
    saved=json.loads(candidate_file.read_text())['results']
    checks=[]
    for rid,rx in original.items():
        request={k:copy.deepcopy(v) for k,v in rx.items() if k in allowed}
        try:base=estimator.score_reaction(request,key=rid)
        except MissingRecord:
            assert saved[rid]['status']=='missing_archived_molecule'
            continue
        if base['dG'] is None:
            assert saved[rid]['status']=='workflow_refusal'
            continue
        assert saved[rid]['status']=='estimated'
        assert abs(base['dG_unrounded']-saved[rid]['result']['dG_unrounded'])<1e-7
        for variant in ['reverse','double','atom_order','species_order']:
            other=copy.deepcopy(request);factor=1
            if variant in ['reverse','double']:
                factor=-1 if variant=='reverse' else 2
                for values in other['species'].values():values[0]*=factor
                other['n_Hplus']=factor*other.get('n_Hplus',0)
            if variant=='species_order':other['species']=dict(reversed(list(other['species'].items())))
            if variant=='atom_order':
                for values in other['species'].values():
                    mol=Chem.MolFromSmiles(values[2]);mol=Chem.RenumberAtoms(mol,list(reversed(range(mol.GetNumAtoms()))))
                    values[2]=Chem.MolToSmiles(mol,canonical=False,isomericSmiles=True)
            result=estimator.score_reaction(other,key=rid+'_'+variant)
            assert result['dG'] is not None,(rid,variant,result['suspect'])
            residual=result['dG_unrounded']-factor*base['dG_unrounded']
            assert abs(residual)<1e-7,(rid,variant,residual)
            assert abs(result['U_samp']-abs(factor)*base['U_samp'])<1e-8
            checks.append(dict(reaction=rid,variant=variant,residual_kj=residual))
    families=[]
    for model in models:
        reference=model['reference_smiles'];parent=None;rows=[]
        try:parent=estimator.state(1,reference)
        except MissingRecord:
            families.append(dict(reference=reference,status='parent_missing'))
            continue
        for charge,smiles in proton_drawings(reference):
            state=estimator.state(charge,smiles)
            residual=state['potential_kj']-(charge-1)*G_HPLUS-parent['potential_kj']
            assert abs(residual)<1e-7,(reference,smiles,residual)
            assert state['sampling_sigma_kj']==parent['sampling_sigma_kj']
            rows.append(dict(charge=charge,smiles=smiles,residual_kj=residual))
        families.append(dict(reference=reference,status='checked',drawings=rows))
    report=dict(scope='Shared-state reaction/drawing invariances only; not physical accuracy.',
        catalog_sha256=hashlib.sha256(catalog_path.read_bytes()).hexdigest(),
        n_estimated_reactions=len(checks)//4,n_reaction_checks=len(checks),checks=checks,
        families=families,n_full_source_families_checked=sum(r['status']=='checked' for r in families),
        n_protonation_drawings_checked=sum(len(r.get('drawings',[])) for r in families),
        maximum_reaction_residual_kj=max(abs(r['residual_kj']) for r in checks),
        no_empirical_offsets=True)
    with output.open('x') as stream:json.dump(report,stream,indent=2,allow_nan=False);stream.write('\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ['checks','families']},indent=2))

if __name__=='__main__':main()
