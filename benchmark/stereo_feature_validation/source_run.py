"""Actual pinned graph features and archived-member contrasts, without model inference."""
from pathlib import Path
from collections import defaultdict
import ast,datetime,gzip,hashlib,importlib.util,json,math,os,sys
H=Path(__file__).resolve().parent;P=H.parents[1]
read=lambda p:json.loads(p.read_text());sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
protocol=read(H/'protocol.json');assert sorted(os.sched_getaffinity(0))==[79]
assert all(os.environ.get(k)=='1' for k in ['OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS'])
for name,digest in protocol['inputs'].items():assert sha(P/name)==digest,name
for name,digest in protocol['sources'].items():assert sha(H/'owned_source'/name)==digest,name
# Inspect the active API call, rather than the unused legacy featurizer/network.
api=ast.parse((H/'owned_source/dGbyG/api/utils.py').read_text());function=next(n for n in api.body if isinstance(n,ast.FunctionDef) and n.name=='predict_standard_dGf_prime')
calls=[n for n in ast.walk(function) if isinstance(n,ast.Call) and isinstance(n.func,ast.Name) and n.func.id=='mol_to_graph_data'];assert len(calls)==1 and len(calls[0].args)==1 and not calls[0].keywords
network_attempts=[]
def forbid_network(event,args):
 if event in {'socket.connect','socket.getaddrinfo'}:
  network_attempts.append(event);raise RuntimeError('Network forbidden in pinned-feature audit')
sys.addaudithook(forbid_network)
import torch
assert not torch.cuda.is_initialized();torch.set_num_threads(1);torch.set_num_interop_threads(1)
from rdkit import Chem,rdBase
import torch_geometric
spec=importlib.util.spec_from_file_location('pinned_dgbyg_graph_features',H/'owned_source/dGbyG/utils/NNFunc.py');features=importlib.util.module_from_spec(spec);spec.loader.exec_module(features)
defaults=features.mol_to_graph_data.__defaults__;assert defaults==(['atomic number','hybridization','aromaticity','charge'],['bond type'])

def mol(s):
 m=Chem.MolFromSmiles(s)
 if m is None:raise ValueError('Invalid source graph')
 for a in m.GetAtoms():a.SetAtomMapNum(0)
 return Chem.AddHs(m)
def stereo_info(m):
 exact=Chem.MolToSmiles(m,isomericSmiles=True)
 mirror=Chem.Mol(m)
 for a in mirror.GetAtoms():a.InvertChirality()
 Chem.AssignStereochemistry(mirror,cleanIt=True,force=True)
 flat=Chem.Mol(m);Chem.RemoveStereochemistry(flat)
 return dict(exact=exact,mirror=Chem.MolToSmiles(mirror,isomericSmiles=True),flat=Chem.MolToSmiles(flat,isomericSmiles=True),
  unspecified=[dict(type=str(s.type),center=int(s.centeredOn),specified=str(s.specified)) for s in Chem.FindPotentialStereo(m) if str(s.specified)!='Specified'],isotopes=any(a.GetIsotope() for a in m.GetAtoms()))
def arrays(m,**kwargs):
 data=features.mol_to_graph_data(m,**kwargs)
 assert data.x.device.type==data.edge_index.device.type==data.edge_attr.device.type=='cpu'
 order=sorted(range(data.edge_index.shape[1]),key=lambda j:tuple(data.edge_index[:,j].tolist()))
 return dict(x=data.x.tolist(),edge_index=data.edge_index[:,order].tolist(),edge_attr=data.edge_attr[order].tolist())
def digest_data(d):return hashlib.sha256(json.dumps(d,sort_keys=True,allow_nan=False).encode()).hexdigest()
all_arrays={}
def pair_record(left_smi,right_smi):
 left,right=mol(left_smi),mol(right_smi);a,b=stereo_info(left),stereo_info(right)
 assert a['flat']==b['flat'] and not a['isotopes'] and not b['isotopes']
 lf,rf=Chem.Mol(left),Chem.Mol(right);Chem.RemoveStereochemistry(lf);Chem.RemoveStereochemistry(rf)
 mapping=rf.GetSubstructMatch(lf,useChirality=False);assert len(mapping)==left.GetNumAtoms()==right.GetNumAtoms()
 right=Chem.RenumberAtoms(right,list(mapping));assert stereo_info(right)['exact']==b['exact']
 x,y=arrays(left),arrays(right);assert x==y
 key=digest_data(x);all_arrays.setdefault(key,x)
 category='same_canonical_state' if a['exact']==b['exact'] else 'unresolved_stereo' if a['unspecified'] or b['unspecified'] else 'enantiomer_pair' if a['mirror']==b['exact'] else 'fully_specified_non_enantiomer_pair'
 return dict(left=a,right=b,right_to_left_atom_order=list(mapping),graph_feature_sha256=key,feature_arrays_identical=True,category=category)
# Controls exercise the real function. Optional chirality exists but is not active.
control=mol('C[C@H](O)C(=O)O');inverted=Chem.Mol(control)
for atom in inverted.GetAtoms():atom.InvertChirality()
assert stereo_info(control)['exact']!=stereo_info(inverted)['exact'] and arrays(control)==arrays(inverted)
with_stereo=list(defaults[0])+['chirality'];assert arrays(control,atom_features=with_stereo)!=arrays(inverted,atom_features=with_stereo)
ez=mol('Cl/C=C/Cl');opposite=mol(r'Cl/C=C\Cl')
assert stereo_info(ez)['exact']!=stereo_info(opposite)['exact'] and arrays(ez)==arrays(opposite)
assert arrays(mol('CCO'))!=arrays(mol('CCN'))
formation=read(P/'benchmark/stereo_reference_feasibility_validation/structural_review.json');formation_rows=[]
for row in formation['pairs']:
 l,r=row['left'],row['right'];data=pair_record(l['smiles'],r['smiles']);formation_rows.append(dict(left_name=l['name'],right_name=r['name'],left_training_row=l['row'],right_training_row=r['row'],**data))
raw=read(P/'benchmark/native_uncertainty_validation/raw_outputs.json');native_groups=defaultdict(list)
for key,row in raw['molecular_vectors'].items():native_groups[stereo_info(mol(row['normalized_SMILES']))['flat']].append((key,row))
native=[]
for group in native_groups.values():
 for i,(key,a) in enumerate(group):
  for other,b in group[i+1:]:
   if stereo_info(mol(a['normalized_SMILES']))['exact']==stereo_info(mol(b['normalized_SMILES']))['exact']:continue
   data=pair_record(a['normalized_SMILES'],b['normalized_SMILES']);differences=[x-y for x,y in zip(a['members_kj'],b['members_kj'],strict=True)];assert len(differences)==100
   assert all(v==0 for v in differences)
   native.append(dict(left_molecule=key,right_molecule=other,member_difference_kj=differences,native_mean_difference_kj=0.,native_population_SD_kj=0.,**data))
# Selection uses chemistry signatures only; no reference or MetaG errors are read.
reaction_rows=[]
for row in raw['raw_reactions']:
 exact=defaultdict(float);flat=defaultdict(float)
 for term in row['source_terms']:
  source=raw['molecular_vectors'][term['molecule']];info=stereo_info(mol(source['normalized_SMILES']));nu=term['coefficient'];exact[info['exact']]+=nu;flat[info['flat']]+=nu
 exact={k:v for k,v in exact.items() if v};flat={k:v for k,v in flat.items() if v}
 if not flat and exact:
  vector=row['members_native_kj'];assert all(v==0 for v in vector)
  reaction_rows.append(dict(reaction=row['reaction'],author675_row=row['author675_row'],source_terms=row['source_terms'],native_neural_members_kj=vector,neural_mean_kj=0.,neural_population_SD_kj=0.,external_condition_bridge_kj=row['condition_bridge_kj'],scope='Neural baseline only; any external condition bridge remains separate. Existing development/training exposure retained.'))
assert not torch.cuda.is_initialized()
count=lambda rows:dict(__import__('collections').Counter(r['category'] for r in rows))
report=dict(schema='metag-stereo-feature-audit-results-v1',utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),protocol_sha256=sha(H/'protocol.json'),source_sha256=sha(Path(__file__)),
 runtime=dict(torch=torch.__version__,PyG=torch_geometric.__version__,rdkit=rdBase.rdkitVersion,CUDA_initialized=False,CPU_affinity=sorted(os.sched_getaffinity(0))),
 default_atom_features=defaults[0],default_bond_features=defaults[1],atom_feature_dimension=len(next(iter(all_arrays.values()))['x'][0]),
 controls=dict(tetrahedral_tag_inversion_default_arrays_identical=True,EZ_inversion_default_arrays_identical=True,optional_chirality_feature_distinguishes_control=True,element_change_distinguished=True),
 formation_pair_counts=count(formation_rows),formation_pairs=formation_rows,native_pair_counts=count(native),native_pairs=native,graph_cancelled_native_reactions=reaction_rows,
 network_forbidden=True,network_attempts=network_attempts,neural_forward_calls=0,weights_loaded=0,training=False,experimental_error_comparison_performed=False,independent_primary_labels_admitted=0,
 scope='Audited default feature map and archived outputs only, not all dGbyG versions or GNNs. Enantiomer equality is physically appropriate in an achiral environment. Other stereo collisions do not by themselves establish nonzero measured energies or a MetaG accuracy advantage. External condition transformations can differ from the neural baseline.')
(H/'report.json').write_text(json.dumps(report,indent=2)+'\n')
with (H/'feature_arrays.json.gz').open('wb') as file:
 with gzip.GzipFile(filename='',mode='wb',fileobj=file,mtime=0) as gz:gz.write(json.dumps(all_arrays,sort_keys=True).encode())
print(json.dumps({k:report[k] for k in ['runtime','default_atom_features','default_bond_features','formation_pair_counts','native_pair_counts','graph_cancelled_native_reactions','controls']},indent=2))
