"""Reproduce both boundary tables from fixed, owned records without quantum software."""
from pathlib import Path
import importlib.abc
import json
import math
import copy
import socket
import sys

HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE/'src'),str(HERE)]
class NoQuantum(importlib.abc.MetaPathFinder):
    def find_spec(self,fullname,path=None,target=None):
        if fullname.split('.')[0] in {'torch','fairchem','pyscf','gpu4pyscf','cupy'}:
            raise RuntimeError('Quantum software forbidden in boundary replay: '+fullname)
sys.meta_path.insert(0,NoQuantum())
def deny(*a,**k):raise RuntimeError('Network forbidden in boundary replay')
socket.socket.connect=deny
from rdkit import Chem
from metag.settings import Settings
from metag.state_function import StateFunctionEstimator
from benchmark.anchored_nicotinamide import AnchoredNicotinamideEstimator
from benchmark.anchored_nicotinamide_regio import RegioAnchoredEstimator

packet=Path(sys.argv[1]);read=lambda p:json.loads(Path(p).read_text())
requests=read(packet/'requests.json');water=read(packet/'water_reference.json')['water_liquid_1atm_kj']
settings=read(packet/'raw_settings.json')
class MissingMolecule(RuntimeError):pass
summary={}
for phase in ['primary','augmented']:
    expected=HERE/'expected'/phase;catalog=read(expected/'raw_catalog.json')
    class Backend:
        def __init__(self):
            self.settings=Settings(**settings);self.model_name=self.settings.model
            self.producer=dict(cacheable=False,kind='cold-physical-record-replay',source_profiles=catalog['source_profiles'])
            self.rows={(r['charge'],r['smiles']):r for r in catalog['records']}
        def species(self,q,smi,label,warnings,diagnostics):
            key=(q,Chem.MolToSmiles(Chem.MolFromSmiles(smi),isomericSmiles=True))
            if key not in self.rows:raise MissingMolecule(str(q)+'|'+key[1])
            row=self.rows[key];diagnostics.update(raw_record=row['name'],raw_sha256=row['sha256'],raw_profile=row['profile'])
            return row['G_kj'],row['sigma_kj']
        def liquid_water(self):return water
    counts={}
    for arm,cls in [('full',StateFunctionEstimator),('anchor14',AnchoredNicotinamideEstimator),('anchor12and14',RegioAnchoredEstimator)]:
        estimator=cls(Backend());saved=read(expected/(arm+'.json'))['results'];count=0
        for rid,row in requests.items():
            try:result=estimator.score_reaction(row['request'],key=rid)
            except MissingMolecule:
                assert saved[rid]['status']=='missing_validated_molecular_prerequisite';continue
            assert (result['dG'] is not None)==(saved[rid]['status']=='estimated'),(phase,arm,rid)
            if result['dG'] is not None:
                assert abs(result['dG_unrounded']-saved[rid]['result']['dG_unrounded'])<1e-7
                assert abs(result['U_samp']-saved[rid]['result']['U_samp'])<1e-8;count+=1
                for variant in ['reverse','double','atom_order','species_order']:
                    other=copy.deepcopy(row['request']);factor=-1 if variant=='reverse' else 2 if variant=='double' else 1
                    if variant in ['reverse','double']:
                        for values in other['species'].values():values[0]*=factor
                        other['n_Hplus']=factor*other.get('n_Hplus',0)
                    if variant=='species_order':other['species']=dict(reversed(list(other['species'].items())))
                    if variant=='atom_order':
                        for values in other['species'].values():
                            m=Chem.MolFromSmiles(values[2]);m=Chem.RenumberAtoms(m,list(reversed(range(m.GetNumAtoms()))))
                            values[2]=Chem.MolToSmiles(m,canonical=False,isomericSmiles=True)
                    check=estimator.score_reaction(other,key=rid+'_'+variant)
                    assert check['dG'] is not None and abs(check['dG_unrounded']-factor*result['dG_unrounded'])<1e-7
                    assert abs(check['U_samp']-abs(factor)*result['U_samp'])<1e-8
        counts[arm]=count
    summary[phase]=counts
assert not any(n.split('.')[0] in {'torch','fairchem','pyscf','gpu4pyscf','cupy'} for n in sys.modules)
print(json.dumps(dict(status='passed',requested_per_arm=12,availability=summary,quantum_imports_forbidden=True,
    network_forbidden=True,scope='Frozen policy and record replay, not quantum recomputation or experimental validation.')))
