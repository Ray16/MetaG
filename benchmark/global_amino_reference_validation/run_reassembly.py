"""Run the three frozen physical reassemblies before comparing observed errors."""
from pathlib import Path
import copy
import hashlib
import json
import os
import sys
import time

HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE/'assembly_source'),str(HERE/'source_helpers')]
from rdkit import Chem
from metag.settings import Settings
from metag.state_function import StateFunctionEstimator
from candidate import GlobalAminoEstimator

read=lambda p:json.loads(Path(p).read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(path,value):
    path=Path(path);tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');tmp.replace(path)

class MissingRecord(RuntimeError):pass

class OwnedBackend:
    def __init__(self,catalog,parent_only=False):
        self.settings=Settings(**read(HERE/'raw_settings.json'))
        self.model_name=self.settings.model
        self.producer={'cacheable':False,'kind':'owned-molecule-reassembly',
            'source_producer':catalog['producer'],'raw_catalog_sha256':sha(HERE/'catalog.json'),
            'parent_only':parent_only}
        self.rows={(r['charge'],r['smiles']):r for r in catalog['records']
                   if not parent_only or r['group']=='parent'}
        self.used={};self.water=catalog['water_liquid_1atm_kj']

    def species(self,charge,smiles,label,warnings,diagnostics):
        smiles=Chem.MolToSmiles(Chem.MolFromSmiles(smiles),isomericSmiles=True)
        key=(int(charge),smiles)
        if key not in self.rows:
            raise MissingRecord(str(charge)+'|'+smiles)
        row=self.rows[key]
        self.used[row['name']]=row['sha256']
        diagnostics.update(method='authenticated-raw-record-replay',source_record=row['name'],
            source_sha256=row['sha256'],n_minima=row['n_minima'])
        return row['G_kj'],row['sigma_kj']

    def liquid_water(self):return self.water

def main():
    manifest=read(HERE/'implementation_freeze.json')
    for name,digest in manifest['files_sha256'].items():assert sha(HERE/name)==digest,name
    for name,digest in read(HERE/'protocol.json')['assembly_files_sha256'].items():
        assert sha(HERE/'assembly_source'/name)==digest,name
    assert sorted(os.sched_getaffinity(0))==[79]
    write(HERE/'reassembly_owner.json',dict(pid=os.getpid(),host=os.uname().nodename,
        start_ticks=Path(f'/proc/{os.getpid()}/stat').read_text().rsplit(')',1)[1].split()[19]))
    catalog=read(HERE/'catalog.json');models=read(HERE/'models_nominal.json')
    original=read(HERE/'inputs.json')
    allowed={'species','n_Hplus','species_state_pools','species_protonation_models','pka_sites'}
    requests={rid:{k:v for k,v in r.items() if k in allowed} for rid,r in original.items()}
    expectations=read(HERE/'baseline_expectations.json')['results']
    root=HERE/'reassembly';root.mkdir(exist_ok=False)
    started=time.monotonic();status={'state':'running','completed_arms':[]}
    write(HERE/'reassembly_status.json',status)
    for arm in ['full_parent589','full_extended594','global_amino_extended594']:
        backend=OwnedBackend(catalog,parent_only=arm=='full_parent589')
        estimator=GlobalAminoEstimator(backend,models) if arm.startswith('global_') else StateFunctionEstimator(backend)
        results={}
        for index,(rid,request) in enumerate(requests.items()):
            try:
                result=estimator.score_reaction(copy.deepcopy(request),key=rid)
                row=dict(status='estimated' if result['dG'] is not None else 'workflow_refusal',result=result)
            except MissingRecord as error:
                row=dict(status='missing_archived_molecule',missing=str(error),result=None)
            results[rid]=row
            if index%50==0:print(arm,index,'/',len(requests),flush=True)
        if arm=='full_parent589':
            for rid,row in results.items():
                prior=expectations[rid];new=row['result']
                estimated=new is not None and new['dG'] is not None
                assert estimated==(prior['dG'] is not None),('baseline mask changed',rid)
                if estimated:
                    assert abs(new['dG_unrounded']-prior['dG_unrounded'])<1e-7,('baseline changed',rid)
                    assert abs(new['U_samp']-prior['sampling_sigma_kj'])<1e-8,('baseline sampling changed',rid)
        write(root/(arm+'.json'),dict(arm=arm,results=results,ledger=estimator.ledger(),
            used_raw_records=backend.used,elapsed_seconds=time.monotonic()-started))
        status['completed_arms'].append(arm);write(HERE/'reassembly_status.json',status)
        print('ARM COMPLETE',arm,sum(r['status']=='estimated' for r in results.values()),flush=True)
    status.update(state='complete',elapsed_seconds=time.monotonic()-started)
    write(HERE/'reassembly_status.json',status)

if __name__=='__main__':main()
