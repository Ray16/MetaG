"""Verify two fresh parents and reassemble with explicitly owned numerical profiles."""
from pathlib import Path
import copy
import hashlib
import json
import math
import shutil
import sys

HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE/'assembly_source'),str(HERE/'source_helpers')]
from rdkit import Chem
from metag.settings import Settings
from metag.state_pools import KT_KJ
from metag.state_function import StateFunctionEstimator
from candidate import GlobalAminoEstimator
from run_reassembly import OwnedBackend,MissingRecord

read=lambda p:json.loads(Path(p).read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(path,value):
    with Path(path).open('x') as f:json.dump(value,f,indent=2,allow_nan=False);f.write('\n')


class ExplicitProfileBackend(OwnedBackend):
    def __init__(self,catalog,path):
        super().__init__(catalog)
        self.producer={'cacheable':False,'kind':'explicit-profile-physical-reassembly',
            'catalog_sha256':sha(path),'source_profiles':catalog['source_profiles']}

    def species(self,charge,smiles,label,warnings,diagnostics):
        value,sigma=super().species(charge,smiles,label,warnings,diagnostics)
        canonical=Chem.MolToSmiles(Chem.MolFromSmiles(smiles),isomericSmiles=True)
        row=self.rows[(charge,canonical)]
        diagnostics['raw_numerical_profile']=row.get('raw_profile','original')
        return value,sigma


def verify_new_row(result,profile,catalog):
    assert result['status']=='molecular_record_complete'
    path=Path(result['cache_path']);assert sha(path)==result['cache_sha256']
    row=read(path);assert row['method']=='implicit' and row['ver']=='v1'
    assert row['q']==result['request']['charge'] and row['smi']==result['request']['smiles']
    assert row['G']==result['energy_kj'] and row['sigma']==result['sampling_sigma_kj']
    original=catalog['producer'];actual=row['settings']['producer']
    differences={k for k in set(original)|set(actual) if original.get(k)!=actual.get(k)}
    assert differences==(set() if profile=='original' else {'source_sha256'})
    assert {k:v for k,v in row['settings'].items() if k!='producer'}=={
        k:v for k,v in catalog['settings'].items() if k!='producer'}
    ensemble=row['meta']['thermal']['ensemble'];minima=ensemble['validated_minima']
    assert minima and ensemble['n_reference_thermal_approximation']==0
    values=[]
    for point in minima:
        geometry=point['geometry'];assert geometry['charge']==row['q'] and geometry['spin_multiplicity']==1
        value=point['electronic_kj']+point['solvation_kj']+point['thermal_kj']
        assert abs(value-point['total_kj'])<1e-7
        values.append(value)
    low=min(values);total=low-KT_KJ*math.log(math.fsum(math.exp(-(v-low)/KT_KJ) for v in values))
    assert abs(total-row['G'])<1e-7
    target=HERE/'verified_new_records'/path.name
    target.parent.mkdir(exist_ok=True);shutil.copyfile(path,target)
    return dict(group='fresh_parent',source=str(target.relative_to(HERE)),sha256=sha(path),name=path.name,
        charge=row['q'],smiles=row['smi'],G_kj=row['G'],sigma_kj=row['sigma'],n_minima=len(minima),
        partition_reconstruction_error_kj=total-row['G'],raw_profile=profile),actual


def main():
    original=read(HERE/'catalog.json');catalog=copy.deepcopy(original)
    for row in catalog['records']:row['raw_profile']='original'
    profiles={'original':{'producer':original['producer'],
        'optimizer_policy':'Original strict-mode protocol; retained records unchanged.'}}
    added=[]
    for file,profile in [('quantum_extension/results/aspartate_R_parent.json','original'),
                         ('recovery_calculation/results/alanine_S_parent.json','curvature_recovery')]:
        row,producer=verify_new_row(read(HERE/file),profile,original)
        assert (row['charge'],row['smiles']) not in {(r['charge'],r['smiles']) for r in catalog['records']}
        catalog['records'].append(row);added.append(row)
        if profile!='original':profiles[profile]=dict(producer=producer,
            optimizer_policy='Additional bounded final LBFGS curvature recovery,400steps/180s,fmax0.001; strict50cm gate unchanged.',
            source_files_sha256=read(HERE/'recovery_calculation/protocol.json')['source_files_sha256'])
    assert len(catalog['records'])==596
    catalog['source_profiles']=profiles
    catalog['scope']='Explicitly mixed original/recovery numerical profiles, same Hamiltonian/components; no old record relabelled.'
    write(HERE/'catalog596.json',catalog)
    write(HERE/'fresh_parent_verification.json',dict(added_records=added,profiles=profiles,
        original_record_count=594,extended_count=596,original_alanine_refusal_preserved=True,
        new_quantum_metabolite_estimates=2,new_reaction_estimates=False))
    inputs=read(HERE/'inputs.json');models=read(HERE/'models_nominal.json')
    output=HERE/'reassembly596';output.mkdir(exist_ok=False)
    allowed={'species','n_Hplus','species_state_pools','species_protonation_models','pka_sites'}
    for arm in ['full_extended596','global_amino_extended596']:
        backend=ExplicitProfileBackend(catalog,HERE/'catalog596.json')
        estimator=GlobalAminoEstimator(backend,models) if arm.startswith('global_') else StateFunctionEstimator(backend)
        estimator.config['molecular_calculation_scope']='Retained record replay; no optimization in this assembly.'
        estimator.config['minimum_optimizer_fallback']={'scope':'read actual per-record profiles',
            'source_profiles':profiles}
        results={}
        for rid,row in inputs.items():
            request={k:copy.deepcopy(v) for k,v in row.items() if k in allowed}
            try:
                result=estimator.score_reaction(request,key=rid)
                results[rid]=dict(status='estimated' if result['dG'] is not None else 'workflow_refusal',result=result)
            except MissingRecord as error:results[rid]=dict(status='missing_archived_molecule',missing=str(error),result=None)
        if arm.startswith('full_'):
            prior=read(HERE/'reassembly/full_extended594.json')['results']
            for rid,row in results.items():
                assert row['status']==prior[rid]['status']
                if row['status']=='estimated':assert abs(row['result']['dG_unrounded']-prior[rid]['result']['dG_unrounded'])<1e-7
        write(output/(arm+'.json'),dict(arm=arm,results=results,ledger=estimator.ledger(),
            used_raw_records=backend.used,catalog_sha256=sha(HERE/'catalog596.json')))
        print(arm,{s:sum(r['status']==s for r in results.values()) for s in {r['status'] for r in results.values()}},flush=True)

if __name__=='__main__':main()
