"""Authenticate every unchanged raw record before reassembly or error scoring."""
from pathlib import Path
import hashlib
import json
import math
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path[:0]=[str(HERE/'assembly_source'),str(HERE/'source_helpers')]
from rdkit import Chem
from metag.state_pools import KT_KJ
from metag.constants import CONDITIONS
from coupled_source_parameters import candidate_models

read=lambda p:json.loads(Path(p).read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(name,value):
    with (HERE/name).open('x') as stream:json.dump(value,stream,indent=2,allow_nan=False);stream.write('\n')

protocol=read(HERE/'protocol.json')
for name,digest in protocol['frozen_files'].items():assert sha(HERE/name)==digest
for name,digest in protocol['assembly_files_sha256'].items():assert sha(HERE/'assembly_source'/name)==digest
old=ROOT/'analysis/metag_corrected_validation_20261030_01'
repair=ROOT/'analysis/metag_physical_repair_20261103_01'
control=read(old/'reference_control_result.json')
producer=control['runtime_producer']
indices=[('parent',old/'shared_record_index.json',old/'shared_record_ledger'),
    ('extension1',repair/'missing_raw/verified_index.json',repair/'missing_raw/verified_records'),
    ('extension2',repair/'missing_raw_wave2/verified_index.json',repair/'missing_raw_wave2/verified_records')]
catalog=[];settings=None;seen=set();parent_keys=[]
for group,index_path,folder in indices:
    index=read(index_path);entries=index if group=='parent' else index['records']
    for name,entry in sorted(entries.items()):
        path=folder/name
        assert sha(path)==entry['sha256'],path
        row=read(path)
        if settings is None:settings=row['settings']
        assert row['settings']==settings and row['settings']['producer']==producer
        assert row['method']=='implicit' and row['ver']=='v1'
        canonical=Chem.MolToSmiles(Chem.MolFromSmiles(row['smi']),isomericSmiles=True)
        key=(row['q'],canonical)
        assert key not in seen,('ambiguous raw state',key)
        seen.add(key)
        ensemble=row['meta']['thermal']['ensemble'];minima=ensemble['validated_minima']
        assert minima and ensemble['n_reference_thermal_approximation']==0
        values=[m['electronic_kj']+m['solvation_kj']+m['thermal_kj'] for m in minima]
        assert all(abs(v-m['total_kj'])<=1e-7 for v,m in zip(values,minima))
        low=min(values);energy=low-KT_KJ*math.log(math.fsum(math.exp(-(v-low)/KT_KJ) for v in values))
        assert abs(energy-row['G'])<=1e-7,path
        assert math.isfinite(row['G']) and math.isfinite(row['sigma']) and row['sigma']>=0
        item=dict(group=group,source=str(path.relative_to(ROOT)),sha256=entry['sha256'],
            name=name,charge=row['q'],smiles=canonical,G_kj=row['G'],sigma_kj=row['sigma'],
            n_minima=len(minima),partition_reconstruction_error_kj=energy-row['G'])
        catalog.append(item)
        if group=='parent':parent_keys.append(name)
assert len(catalog)==594 and len(parent_keys)==589
models=candidate_models(read(HERE/'source_transcription.json'),CONDITIONS)
assert len(models)==7
write('models_nominal.json',models)
write('catalog.json',dict(records=catalog,parent_keys=parent_keys,settings=settings,
    producer=producer,water_liquid_1atm_kj=control['water_liquid_1atm_kj'],
    source_indices=[dict(group=g,path=str(p.relative_to(ROOT)),sha256=sha(p)) for g,p,_ in indices]))
write('preparation_receipt.json',dict(records=594,parent_records=589,
    catalog_sha256=sha(HERE/'catalog.json'),models_sha256=sha(HERE/'models_nominal.json'),
    helper_files_sha256={n:sha(HERE/'source_helpers'/n) for n in protocol['helper_files_sha256']},
    inputs_sha256=sha(HERE/'inputs.json'),protocol_sha256=sha(HERE/'protocol.json'),
    production_runtime_changed=False,new_molecular_work=False))
print('Authenticated594raw records and7exact source reference models; no error scoring.')
