"""Extract the released ensemble on the fixed development cohort, using CPU only."""
from pathlib import Path
from collections import defaultdict
import csv
import datetime
import hashlib
import json
import math
import os
import socket
import sys
import time

HERE=Path(__file__).resolve().parent
sha=lambda path:hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(name):return json.loads((HERE/name).read_text())
def write(name,value):
    target=HERE/name;temp=target.with_name(target.name+f'.tmp-{os.getpid()}')
    temp.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');temp.replace(target)
protocol=read('protocol.json');start=time.monotonic()
assert sorted(os.sched_getaffinity(0))==protocol['CPU_affinity']
assert all(os.environ.get(key)=='1' for key in ['OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'])
assert max(os.getloadavg())+1<=os.cpu_count()
for path in Path('/proc').glob('[0-9]*'):
    try:
        if int(path.name)==os.getpid():continue
        argv=(path/'cmdline').read_bytes().split(b'\0')
        if argv and b'python' in Path(argv[0].decode()).name.encode():
            cpus=set(os.sched_getaffinity(int(path.name)))
            assert not(len(cpus)<os.cpu_count() and cpus&set(protocol['CPU_affinity'])), f'CPU conflict with {path.name}'
    except (OSError,UnicodeError,ValueError):pass
assert int(Path('/proc/meminfo').read_text().split('MemAvailable:')[1].split()[0])*1024>=32*1024**3
write('predictor_owner.json',dict(host=os.uname().nodename,pid=os.getpid(),
    start_ticks=Path(f'/proc/{os.getpid()}/stat').read_text().rsplit(')',1)[1].split()[19],
    cpus=sorted(os.sched_getaffinity(0)),utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
freeze=read('freeze.json')
for name,digest in freeze['files_sha256'].items():assert sha(HERE/name)==digest,name
source=Path(protocol['source_root'])
for name,asset in protocol['source_assets'].items():
    assert sha(source/name)==asset['sha256'] and (source/name).stat().st_size==asset['bytes'],name
for name,digest in protocol['source_inputs'].items():assert sha(HERE/name)==digest,name
sys.path.insert(0,str(source))
import torch
torch.set_num_threads(1);torch.set_num_interop_threads(1)
from rdkit import Chem,rdBase
from rdkit.Chem import rdChemReactions
import numpy as np

# Preserve the actual upstream directory order while binding every output column
# to its checkpoint. Only device selection changes; no source or weights change.
model_dir=source/'dGbyG/network/model_params/MPNN_model_params'
loaded=[];original_listdir=os.listdir;original_available=torch.cuda.is_available
def record_listdir(path):
    result=original_listdir(path)
    if Path(path).resolve()==model_dir.resolve():loaded.append(list(result))
    return result
os.listdir=record_listdir;torch.cuda.is_available=lambda:False
original_connect=socket.socket.connect
def no_network(*args,**kwargs):raise RuntimeError('Network lookup is forbidden in native ensemble extraction')
socket.socket.connect=no_network
try:
    from dGbyG.api.Chemistry import Compound,Reaction
    from dGbyG.api import utils
    from dGbyG.utils.ChemFunc import parse_equation,to_mol
    from dGbyG.utils.constants import default_condition
finally:
    os.listdir=original_listdir;torch.cuda.is_available=original_available
assert utils.device=='cpu' and not torch.cuda.is_initialized()
assert len(loaded)==1 and len(loaded[0])==100
assert set(loaded[0])=={f'{i}.pt' for i in range(100)}
assert default_condition==dict(T=298.15,pH=7.,I=.25,pMg=14.,e_potential=0.)
members=[dict(file=name,sha256=sha(model_dir/name)) for name in loaded[0]]
compound_cache={};object_ids={};molecular_vectors={}
def reaction(equation):
    parts=parse_equation(equation)
    for smiles in parts:
        if smiles not in compound_cache:
            compound_cache[smiles]=Compound(to_mol(smiles,'smiles'))
    return Reaction({compound_cache[s]:v for s,v in parts.items()},mol_type='compound')
def key(rxn):
    sums=defaultdict(float)
    for compound,coefficient in rxn.reaction.items():sums[compound.Smiles]+=float(coefficient)
    return tuple(sorted((smiles,coefficient) for smiles,coefficient in sums.items() if coefficient))

training=list(csv.DictReader((HERE/'reference671_training.csv').open()))
training_keys=defaultdict(list);training_normalization_failures=[]
for index,row in enumerate(training):
    if row['type']!='reaction dg':continue
    try:training_keys[key(reaction(row['reaction']))].append(index)
    except Exception as error:training_normalization_failures.append(dict(row=index,error=str(error)))
source_rows=list(csv.DictReader((HERE/'author675_structural.csv').open()))
cohort=read('development_inputs.json');assert len(cohort['rows'])==69
rows=[];raw_reactions=[];refusals=[]
for row in cohort['rows']:
    assert time.monotonic()-start<protocol['maximum_runtime_seconds']
    assert int(Path('/proc/self/status').read_text().split('VmRSS:')[1].split()[0])*1024<32*1024**3
    rid=row['reaction'];index=row['author_row'];equation=source_rows[index]['reaction']
    print(f'START {rid} author675row={index}',flush=True)
    try:
        rxn=reaction(equation)
        assert rxn.balanced, 'Native reaction balance gate refused'
        values=rxn.standard_dGr_prime_list
        assert values is not None and values.shape==(100,) and np.isfinite(values).all()
        api_mean,api_sd=map(float,rxn.standard_dGr_prime)
        contributions=[];terms=[]
        for compound,coefficient in rxn.reaction.items():
            values_i=compound.standard_dGf_prime_list
            assert values_i.shape==(100,) and np.isfinite(values_i).all()
            if id(compound) not in object_ids:
                identity=f'mol_{len(object_ids):04d}';object_ids[id(compound)]=identity
                molecular_vectors[identity]=dict(normalized_SMILES=compound.Smiles,
                    members_kj=[float(v) for v in values_i],
                    scope='Arbitrary member formation baselines; marginal mean/SD is not a calibrated molecular property.')
            identity=object_ids[id(compound)]
            contributions.append((float(coefficient),[float(v) for v in values_i]))
            terms.append(dict(molecule=identity,coefficient=float(coefficient)))
        direct=[math.fsum(c*v[j] for c,v in contributions) for j in range(100)]
        # Bound only upstream float32 summation roundoff, not prediction error.
        rounding_bounds=[max(1e-6,(len(terms)+2)*np.finfo(np.float32).eps*math.fsum(abs(c*v[j]) for c,v in contributions)) for j in range(100)]
        assert all(abs(a-float(b))<=bound for a,b,bound in zip(direct,values,rounding_bounds))
        source_vector=[float(value) for value in values]
        scale=row['signed_author_scale'];bridge=row['condition_bridge_kj']
        target_vector=[value/scale+bridge for value in source_vector]
        mean=math.fsum(target_vector)/100
        sd=math.sqrt(math.fsum((value-mean)**2 for value in target_vector)/100)
        raw_reactions.append(dict(reaction=rid,author675_row=index,equation=equation,
            normalized_native_equation=dict(key(rxn)),source_terms=terms,
            members_native_kj=source_vector,signed_author_scale=scale,condition_bridge_kj=bridge,
            API_native_mean_kj=api_mean,API_native_population_SD_kj=api_sd,
            maximum_float32_reconstruction_difference_kj=max(abs(a-b) for a,b in zip(direct,source_vector)),
            maximum_float32_roundoff_bound_kj=float(max(rounding_bounds)),
            exact_current671_normalized_equation_rows=training_keys.get(key(rxn),[])))
        rows.append(dict(reaction=rid,reference=row['reference'],metag=row['metag'],
            comparator=mean,spread=sd,n_members=100))
        print(f'COMPLETE {rid}: native transformed mean={mean:.6g},SD={sd:.6g}',flush=True)
    except Exception as error:
        refusals.append(dict(reaction=rid,error=str(error)))
        print(f'REFUSED {rid}: {error}',flush=True)
    write('partial_outputs.json',dict(rows=rows,raw_reactions=raw_reactions,molecular_vectors=molecular_vectors,
        refused=refusals,members=members,source_version=dict(torch=torch.__version__,numpy=np.__version__,rdkit=rdBase.rdkitVersion)))
assert not torch.cuda.is_initialized()
write('raw_outputs.json',dict(rows=rows,raw_reactions=raw_reactions,molecular_vectors=molecular_vectors,
    refused=refusals,members=members,n_requested=69,n_predicted=len(rows),
    training_normalization_failures=training_normalization_failures,
    source_version=dict(torch=torch.__version__,numpy=np.__version__,rdkit=rdBase.rdkitVersion),
    device=utils.device,CUDA_initialized=False,elapsed_seconds=time.monotonic()-start,
    experimental_validation_regime='development_and_potential_training_exposed',
    independently_validated_risk_improvement=False,new_primary_references_admitted=0))
if refusals:raise SystemExit('Incomplete69cohort; comparative endpoint withheld')
write('native_inputs.json',dict(schema='metag-reliability-input-v1',cohort_exposure=protocol['known_exposure'],
    spread_definition='Population SD of100released model reaction-member sums, normalized by absolute extent. Deterministic bridge uncertainty excluded.',
    severe_threshold_kj=20.,rows=rows,cluster_definition=None,protocol_sha256=sha(HERE/'protocol.json'),scope=protocol['scope']))
print('COMPLETE fixed69cohort; CUDA initialized=False',flush=True)
