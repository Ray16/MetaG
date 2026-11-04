"""Four owned endpoints; no coordinate optimization or experimental target use."""
from pathlib import Path
from collections import Counter
from datetime import datetime, timezone
import hashlib,json,shutil,re
import numpy as np
from ase import Atoms
from rdkit import Chem
from metag.energetics.conformers import bond_graph,same_connectivity,same_stereochemistry
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
assert (ROOT/'src/metag').is_dir() and ROOT.name=='MetaG'
OPT=ROOT/'analysis/aqueous_relaxation_pilot_20261104_02'
HESS=ROOT/'analysis/aqueous_hessian_control_20261104_02'
PARENT=ROOT/'analysis/nicotinamide_reference_control_20261104_01'
assert not (HERE/'launch_manifest.json').exists()
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
read=lambda p:json.loads(p.read_text())
def write(p,d):p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
parents=HERE/'parents';parents.mkdir(exist_ok=True)
source_records={}
def copy(source,name):
 target=parents/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,target)
 source_records[str(target.relative_to(HERE))]=dict(original_path=str(source),sha256=sha(target))
 return read(target) if target.suffix=='.json' else target
inputs=copy(OPT/'geometry_inputs.json','initial_geometry_inputs.json')['states']
diag=copy(HESS/'matched_profile_diagnostic.json','matched_profile_diagnostic.json')
copy(OPT/'geometry_audit.json','geometry_audit.json')
copy(HESS/'protocol.json','aqueous_hessian_protocol.json')
copy(HESS/'independent_hessian_audit.json','independent_hessian_audit.json')
weights=dict(oxidized_gas=1,reduced_gas=-1,oxidized_aqueous=-1,reduced_aqueous=1)
states={};trace_records={}
for slot in weights:
 candidate=copy(OPT/f'{slot}.json',f'{slot}.candidate.json')
 hessian=copy(HESS/f'{slot}.json',f'{slot}.hessian.json')
 hessian_archive=copy(HESS/f'{slot}.npz',f'{slot}.hessian.npz')
 assert hessian['candidate_sha256']==sha(OPT/f'{slot}.json')
 assert candidate['optimization_converged'] and candidate['identity_passed']
 # Last trajectory geometry, energy and charge must belong to the retained candidate.
 trace_path=OPT/f'{slot}.trace.jsonl'
 trace=json.loads(trace_path.read_text().splitlines()[-1])
 assert trace['positions_A']==candidate['positions_A'] and trace['symbols']==candidate['symbols']
 assert trace['charge']==candidate['charge']
 assert abs(trace['UMA_energy_eV']-candidate['UMA_energy_eV'])<1e-12
 trace_records[slot]=dict(original_path=str(trace_path),sha256=sha(trace_path),last_record_sha256=hashlib.sha256(json.dumps(trace,sort_keys=True).encode()).hexdigest())
 write(parents/f'{slot}.last_trace.json',trace)
 source_records[f'parents/{slot}.last_trace.json']=dict(original_path=str(trace_path)+'#last-record',sha256=sha(parents/f'{slot}.last_trace.json'))
 template_source=inputs[slot]
 m=Chem.AddHs(Chem.MolFromSmiles(template_source['parent_reference_smiles']));template=Chem.RWMol(m)
 for index in template_source['deprotonated_O_indices']:
  template.GetAtomWithIdx(index).SetFormalCharge(-1);template.GetAtomWithIdx(index).SetNoImplicit(True)
 for index in sorted(template_source['removed_H_indices'],reverse=True):template.RemoveAtom(index)
 template=template.GetMol();Chem.SanitizeMol(template)
 assert [a.GetSymbol() for a in template.GetAtoms()]==candidate['symbols']
 assert Chem.GetFormalCharge(template)==candidate['charge']
 atoms=Atoms(candidate['symbols'],positions=candidate['positions_A'])
 assert same_connectivity(atoms,bond_graph(template)) and same_stereochemistry(atoms,template)
 geometry=dict(schema='ordered-geometry-v1',symbols=candidate['symbols'],positions_angstrom=candidate['positions_A'],charge=candidate['charge'],spin_multiplicity=1)
 geometry['sha256']=hashlib.sha256(json.dumps(geometry,sort_keys=True,allow_nan=False).encode()).hexdigest()
 raw_name=f'{slot}.matched_scoring.raw.jsonl' if slot.endswith('_gas') else f'{slot}.base.xtb.jsonl'
 raw_path=copy(HESS/raw_name,f'{slot}.alpb.raw.jsonl')
 records=[json.loads(line) for line in raw_path.read_text().splitlines()]
 assert len(records)==2 and {r['phase'] for r in records}=={'gas','aqueous'}
 energies={}
 for record in records:
  assert record['returncode']==0 and record['charge']==candidate['charge']
  assert 'convergence criteria satisfied' in record['stdout']
  argv=record['argv'];assert argv[argv.index('--acc')+1]=='0.001' and '--norestart' in argv
  xyz=record['input_xyz'].splitlines();assert int(xyz[0])==len(geometry['symbols'])
  assert [line.split()[0] for line in xyz[2:]]==geometry['symbols']
  assert [[float(v) for v in line.split()[1:]] for line in xyz[2:]]==geometry['positions_angstrom']
  match=re.search(r'TOTAL ENERGY\s+([+-]?[\d.]+)\s+Eh',record['stdout']);assert match
  energies[record['phase']]=float(match[1])
 delta=(energies['aqueous']-energies['gas'])*2625.4996394799
 with np.load(hessian_archive) as archive:
  assert archive['positions_A'].tolist()==geometry['positions_angstrom']
  replayed_uma=float(archive['energy_UMA_eV'])
 assert abs((replayed_uma-candidate['UMA_energy_eV'])*96.485-hessian['UMA_energy_replay_difference_kj'])<1e-9
 uma=replayed_uma*96.485
 assert abs(uma+delta-diag['states'][slot]['combined_potential_kj'])<1e-6
 states[slot]=dict(geometry=geometry,original_candidate_sha256=sha(OPT/f'{slot}.json'),
  source_candidate=f'parents/{slot}.candidate.json',reference_smiles=template_source['deprotonated_smiles'],
  n_AO=1916 if slot.startswith('oxidized') else 1925,n_electrons=344 if slot.startswith('oxidized') else 346,
  UMA_electronic_kj=uma,original_endpoint_UMA_electronic_kj=candidate['UMA_energy_eV']*96.485,
  owned_AD_replay_minus_endpoint_kj=hessian['UMA_energy_replay_difference_kj'],
  ALPB_increment_kj=delta,matched_combined_potential_kj=uma+delta,
  inherited_single_minimum_RRHO_kj=diag['states'][slot]['RRHO_kj'],
  current_geometry_identity_verified=True)
residual=Counter();charge=0
for slot,w in weights.items():
 for element,count in Counter(states[slot]['geometry']['symbols']).items():residual[element]+=w*count
 charge+=w*states[slot]['geometry']['charge']
assert all(v==0 for v in residual.values()) and charge==0 and sum(weights.values())==0
for family in ['oxidized','reduced']:
 assert states[family+'_gas']['geometry']['symbols']==states[family+'_aqueous']['geometry']['symbols']
# Scientific settings are inherited byte-for-byte; the four new geometries are a separate control.
parent_numeric=copy(PARENT/'fixed_geometry_protocol.json','parent_fixed_geometry_protocol.json')
numeric_keys=['functional','basis','density_fitting','auxbasis','grid_level','nlc_grid_level',
 'SCF_tolerance_Eh','SCF_gradient_tolerance','stationary_gradient_maximum','max_cycle',
 'overlap_minimum_eigenvalue','max_memory_MB','solvent','SMD_sasa_ng','SMD_equilibrium_solvation',
 'Hartree_to_kj_per_mol','numerical_smoke_tolerance_kj']
numeric={key:parent_numeric[key] for key in numeric_keys}
numeric.update(schema='metag-fixed-geometry-DFT-settings-v1',
 scope='Four owned charged-NAD gas/aqueous optimized endpoints. Numerical settings inherited unchanged from the earlier fixed-geometry control.',
 parent_protocol_sha256=sha(PARENT/'fixed_geometry_protocol.json'))
write(HERE/'fixed_geometry_protocol.json',numeric)
shutil.copy2(PARENT/'reference_point.py',HERE/'reference_point.py')
shutil.copy2(PARENT/'cpu_water_smoke.json',HERE/'cpu_water_smoke.json')
shutil.copy2(PARENT/'cpu_smoke.py',HERE/'cpu_smoke.py')
assert read(HERE/'cpu_water_smoke.json')['source_sha256']==sha(HERE/'cpu_smoke.py')
write(HERE/'geometry_inputs.json',dict(inputs=states,source_records=source_records,source_trace_records=trace_records))
protocol=dict(schema='metag-aqueous-geometry-parent-control-v1',created_utc=datetime.now(timezone.utc).isoformat(),
 question='Does the parentDFT/SMD model support the direction and magnitude of the owned UMA/ALPB geometry-intervention energy effect?',
 expected_slots=list(weights),weights=weights,geometry_rule='All four previously converged identity-preserving chargedNAD endpoints, unchanged ordered coordinates; no geometry selected or altered using a new quantum result.',
 endpoint='(reduced_aqueous-oxidized_aqueous)-(reduced_gas-oxidized_gas), separately under each scoring arm',
 original_matched_UMA_ALPB_intervention_kj=sum(weights[s]*states[s]['matched_combined_potential_kj'] for s in states),
 numerical_protocol_sha256=sha(HERE/'fixed_geometry_protocol.json'),parent_worker_source_sha256=sha(HERE/'reference_point.py'),
 basis_and_functional='Restricted-singlet wb97m-v/def2-tzvpd with def2-universal-jkfit; grids4/3; SMD water590SWIG; all AO directions admitted at overlap>=1e-8 retained.',
 arms=['UMA_ALPB','DFT_ALPB','UMA_SMD','DFT_SMD'],
 UMA_energy_basis='Owned exact-autograd-Hessian replay energies used by the matched scoring profile. Original endpoint energies and their sub0.0003kJ replay differences are retained separately.',
 resource_ceiling=dict(host='lambda2',gpu_workers=1,gpu_index=7,cpus=[2,3,4,5],native_threads=1,overall_seconds=5400,attempts_per_geometry=1,consecutive_actual_failures=3),
 element_residual={k:v for k,v in residual.items() if v},charge_residual=charge,molecule_residual=sum(weights.values()),
 scope='Fixed-microgeometry electronic/solvent sensitivity, not a closed equilibrium-metabolite cycle or reactionGibbs energy. Distinct geometry constraints are distinct microstates.',
 limitations=['One geometry per optimized state; no aqueous conformational ensemble.','Gas anion HOMO and orbital-gradient diagnostics retained; no bound-gas-ion or orbital-stability proof.','SCF orbital stationarity is not nuclear-force convergence; these are not parentDFT-optimized geometries.','No new parent Hessian or RRHO evaluation; inherited thermal terms remain a separate model sensitivity only.','No pH, Mg or experimental-reference admission and no default production promotion.'],
 empirical_offsets_fitted=0,experimental_labels_admitted=0,production_runtime_changed=False,
 source_preparer_sha256=sha(Path(__file__)))
write(HERE/'protocol.json',protocol)
print(json.dumps(dict(states={s:dict(n_atoms=len(d['geometry']['symbols']),charge=d['geometry']['charge'],n_AO=d['n_AO'],n_electrons=d['n_electrons']) for s,d in states.items()},geometry_effect_kj=protocol['original_matched_UMA_ALPB_intervention_kj'],all_identity_and_raw_energy_bindings=True),indent=2))
