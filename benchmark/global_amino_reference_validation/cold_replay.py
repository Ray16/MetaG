"""Reassemble the frozen research candidate in an isolated, CPU-only process."""
from pathlib import Path
import copy
import importlib.abc
import json
import sys

HERE=Path(__file__).resolve().parent
assert 'metag' not in sys.modules
sys.path[:0]=[str(HERE/'assembly_source'),str(HERE/'source_helpers')]

class BlockQuantum(importlib.abc.MetaPathFinder):
    def find_spec(self,fullname,path=None,target=None):
        if fullname.split('.')[0] in {'torch','fairchem','pyscf','gpu4pyscf','cupy'}:
            raise RuntimeError('Quantum runtime import forbidden in cold reassembly: '+fullname)
        return None

sys.meta_path.insert(0,BlockQuantum())
from candidate import GlobalAminoEstimator
from augment_and_reassemble import ExplicitProfileBackend
from run_reassembly import MissingRecord

read=lambda p:json.loads(Path(p).read_text())
catalog=read(HERE/'catalog596.json');inputs=read(HERE/'inputs.json')
models=read(HERE/'models_nominal.json')
backend=ExplicitProfileBackend(catalog,HERE/'catalog596.json')
estimator=GlobalAminoEstimator(backend,models)
expected=read(HERE/'reassembly596/global_amino_extended596.json')['results']
allowed={'species','n_Hplus','species_state_pools','species_protonation_models','pka_sites'}
estimates=0
for rid,row in inputs.items():
    try:
        result=estimator.score_reaction({k:copy.deepcopy(v) for k,v in row.items() if k in allowed},key=rid)
        status='estimated' if result['dG'] is not None else 'workflow_refusal'
    except MissingRecord:
        status='missing_archived_molecule';result=None
    assert status==expected[rid]['status'],rid
    if status=='estimated':
        estimates+=1
        assert abs(result['dG_unrounded']-expected[rid]['result']['dG_unrounded'])<1e-7,rid
        assert abs(result['U_samp']-expected[rid]['result']['U_samp'])<1e-8,rid
assert estimates==148
assert not any(name.split('.')[0] in {'torch','fairchem','pyscf','gpu4pyscf','cupy'} for name in sys.modules)
print(json.dumps(dict(status='passed',requests=364,reproduced_estimates=estimates,
    quantum_imports=0,scope='Frozen candidate scalar/availability replay; no quantum calculation or accuracy claim.')))
