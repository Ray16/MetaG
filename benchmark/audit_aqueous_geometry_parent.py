"""Independent per-point energy ledger and fixed-geometry factorial reconstruction."""
from pathlib import Path
from collections import Counter
import argparse,hashlib,json,math,re
import numpy as np
SLOTS=('oxidized_gas','reduced_gas','oxidized_aqueous','reduced_aqueous')
WEIGHTS=(1,-1,-1,1)
Z={'H':1,'C':6,'N':7,'O':8,'P':15}
HARTREE=2625.4996394799
EV=96.485
BOHR_A=0.52917721092
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(path):
 def reject(value):raise ValueError('Nonfinite JSON value: '+value)
 return json.loads(Path(path).read_text(),parse_constant=reject)

def finite(value,label):
 value=float(value)
 if not math.isfinite(value):raise ValueError('Nonfinite '+label)
 return value

def nuclear_energy(geometry):
 positions=np.asarray(geometry['positions_angstrom'],dtype=float)
 if positions.shape!=(len(geometry['symbols']),3) or not np.isfinite(positions).all():raise ValueError('Invalid geometry shape/coordinates')
 result=0.0
 for i,s in enumerate(geometry['symbols']):
  for j in range(i):
   distance=float(np.linalg.norm(positions[i]-positions[j]))
   if distance<=0:raise ValueError('Coincident atoms')
   result+=Z[s]*Z[geometry['symbols'][j]]*BOHR_A/distance
 return result

def validate_point(item,result,log,numeric):
 geometry=item['geometry']
 if result['geometry']!=geometry:raise ValueError('Point geometry does not belong to declared input')
 if result['n_AO']!=item['n_AO'] or result['n_electrons']!=item['n_electrons']:raise ValueError('Basis/electron identity mismatch')
 if sum(Z[s] for s in geometry['symbols'])-geometry['charge']!=result['n_electrons']:raise ValueError('Charge/electron mismatch')
 overlap=finite(result['overlap_minimum_eigenvalue'],'overlap')
 if overlap<numeric['overlap_minimum_eigenvalue']:raise ValueError('Overlap basis quality refusal')
 matches=re.findall(r'converged SCF energy\s*=\s*([+-]?[\d.]+(?:[eEdD][+-]?\d+)?)',log)
 if len(matches)!=2:raise ValueError('Expected gas and water native SCF completion records')
 nuclear=nuclear_energy(geometry);warnings=[]
 for phase,text in zip(['gas','water'],matches):
  record=result[phase];energy=finite(record['energy_Eh'],phase+' energy')
  if not record['converged']:raise ValueError('Unconverged orbital solution')
  if record['n_molecular_orbitals']!=result['n_AO']:raise ValueError('AO directions dropped')
  if not 0<=finite(record['orbital_gradient_norm'],'orbital gradient')<=numeric['stationary_gradient_maximum']:raise ValueError('SCF orbital-gradient quality refusal')
  terms=record['scf_summary']
  # e2 duplicates coul+exc, and gap is a diagnostic; neither is added again.
  total=sum(finite(terms.get(k,0),k) for k in ['e1','coul','exc','nuc','e_solvent','e_cds'])
  if abs(total-energy)>1e-8:raise ValueError('Per-point SCF ledger mismatch')
  if abs(finite(terms['nuc'],'nuclear term')-nuclear)>1e-7:raise ValueError('Nuclear term does not belong to ordered coordinates')
  if abs(float(text.replace('D','E').replace('d','e'))-energy)>1e-8:raise ValueError('Native SCF log differs from claimed point energy')
  homo=finite(record['HOMO_Eh'],'HOMO');lumo=finite(record['LUMO_Eh'],'LUMO')
  if lumo<homo:raise ValueError('Reversed reported orbital gap')
  if phase=='gas' and homo>0:warnings.append('Positive gasHOMO; no stable isolated-anion claim')
 if abs(finite(result['SMD_delta_Eh'],'SMD increment')-(result['water']['energy_Eh']-result['gas']['energy_Eh']))>1e-10:raise ValueError('SMD increment mismatch')
 if abs(finite(result['CDS_Eh'],'CDS term')-result['water']['scf_summary']['e_cds'])>1e-12:raise ValueError('CDS convention mismatch')
 if result['surface_discretization_method']!='SWIG':raise ValueError('Surface model changed')
 return warnings

def audit(directory):
 directory=Path(directory);inputs=read(directory/'geometry_inputs.json');protocol=read(directory/'protocol.json');numeric=read(directory/'fixed_geometry_protocol.json')
 manifest=read(directory/'launch_manifest.json')
 for name,digest in manifest['files_sha256'].items():
  if sha(directory/name)!=digest:raise ValueError('Owned input/producer changed: '+name)
 if protocol['weights']!={s:w for s,w in zip(SLOTS,WEIGHTS)}:raise ValueError('Intervention orientation changed')
 if set(inputs['inputs'])!=set(SLOTS) or manifest['expected_slots']!=list(SLOTS):raise ValueError('Unreported or missing declared endpoint')
 versions=read(directory/'dependency_identity.json')['versions']
 fixture=read(directory/'cpu_water_smoke.json');smoke=read(directory/'points/smoke.json')
 expected_water=dict(symbols=[row[0] for row in fixture['geometry']],positions_angstrom=[row[1] for row in fixture['geometry']],charge=0,spin_multiplicity=1)
 if fixture['source_sha256']!=sha(directory/'cpu_smoke.py'):raise ValueError('CPU fixture source mismatch')
 if smoke['versions']!=versions or smoke.get('passed') is not True:raise ValueError('GPU fixture producer/refusal mismatch')
 validate_point(dict(geometry=expected_water,n_AO=58,n_electrons=10),smoke,(directory/'points/smoke.scf.log').read_text(),numeric)
 fixture_differences={phase:abs(smoke[phase]['energy_Eh']-fixture[key])*HARTREE for phase,key in [('gas','gas_Eh'),('water','water_total_Eh')]}
 if max(fixture_differences.values())>numeric['numerical_smoke_tolerance_kj']:raise ValueError('GPU/CPU numerical fixture disagreement')
 residual=Counter();charge=0;rows=[];nuclear_warnings=[]
 for slot,w in zip(SLOTS,WEIGHTS):
  item=inputs['inputs'][slot];geometry=item['geometry'];point=read(directory/'points'/(slot+'.json'))
  if point.get('passed') is not True:raise ValueError('Incomplete quantum point: '+slot)
  if point['versions']!=versions:raise ValueError('Point dependency versions differ from frozen producer')
  if point['source_sha256']!=sha(directory/'reference_point.py'):raise ValueError('Point producer changed')
  if point['input_sha256']!=sha(directory/'geometry_inputs.json') or point['protocol_sha256']!=sha(directory/'fixed_geometry_protocol.json'):raise ValueError('Point input/settings ownership mismatch')
  warnings=validate_point(item,point,(directory/'points'/(slot+'.scf.log')).read_text(),numeric)
  parent=read(directory/item['source_candidate'])
  if sha(directory/item['source_candidate'])!=item['original_candidate_sha256']:raise ValueError('Candidate source identity mismatch')
  if geometry['positions_angstrom']!=parent['positions_A'] or geometry['charge']!=parent['charge'] or geometry['symbols']!=parent['symbols']:raise ValueError('Candidate geometry changed')
  with np.load(directory/'parents'/(slot+'.hessian.npz')) as archive:
   if archive['positions_A'].tolist()!=geometry['positions_angstrom']:raise ValueError('Owned UMA archive geometry mismatch')
   uma=float(archive['energy_UMA_eV'])*EV
  if abs(uma-item['UMA_electronic_kj'])>1e-8:raise ValueError('Owned UMA electronic value mismatch')
  raw=[json.loads(line) for line in (directory/'parents'/(slot+'.alpb.raw.jsonl')).read_text().splitlines()]
  if len(raw)!=2 or {r['phase'] for r in raw}!={'gas','aqueous'}:raise ValueError('Incomplete native solvent pair')
  energy={}
  for r in raw:
   if r['returncode'] or r['charge']!=geometry['charge']:raise ValueError('Native solvent execution/input mismatch')
   xyz=r['input_xyz'].splitlines()
   if [x.split()[0] for x in xyz[2:]]!=geometry['symbols'] or [[float(v) for v in x.split()[1:]] for x in xyz[2:]]!=geometry['positions_angstrom']:raise ValueError('Native solvent geometry mismatch')
   match=re.search(r'TOTAL ENERGY\s+([+-]?[\d.]+)\s+Eh',r['stdout'])
   if not match or 'convergence criteria satisfied' not in r['stdout']:raise ValueError('No native solvent convergence/energy')
   energy[r['phase']]=float(match[1])
  alpb=(energy['aqueous']-energy['gas'])*HARTREE
  if abs(alpb-item['ALPB_increment_kj'])>1e-8:raise ValueError('ALPB source increment mismatch')
  dft_gas=point['gas']['energy_Eh']*HARTREE;smd=(point['water']['energy_Eh']-point['gas']['energy_Eh'])*HARTREE
  arms=dict(UMA_ALPB=uma+alpb,DFT_ALPB=dft_gas+alpb,UMA_SMD=uma+smd,DFT_SMD=point['water']['energy_Eh']*HARTREE)
  for element,count in Counter(geometry['symbols']).items():residual[element]+=w*count
  charge+=w*geometry['charge']
  rows.append(dict(slot=slot,weight=w,components_kj=dict(UMA_electronic=uma,DFT_electronic=dft_gas,ALPB=alpb,SMD=smd),arm_potentials_kj=arms,gas_HOMO_Eh=point['gas']['HOMO_Eh'],water_HOMO_Eh=point['water']['HOMO_Eh'],n_AO=point['n_AO'],n_electrons=point['n_electrons'],overlap_min=point['overlap_minimum_eigenvalue'],limitations=warnings))
 if any(residual.values()) or charge or sum(WEIGHTS):raise ValueError('Element/charge/stoichiometric intervention imbalance')
 contrasts={arm:sum(r['weight']*r['arm_potentials_kj'][arm] for r in rows) for arm in rows[0]['arm_potentials_kj']}
 per_family={}
 for family in ['oxidized','reduced']:
  gas=next(r for r in rows if r['slot']==family+'_gas');aqueous=next(r for r in rows if r['slot']==family+'_aqueous')
  per_family[family]={arm:aqueous['arm_potentials_kj'][arm]-gas['arm_potentials_kj'][arm] for arm in contrasts}
 result=dict(schema='metag-aqueous-geometry-parent-audit-v1',passed=True,source_sha256=sha(Path(__file__)),protocol_sha256=sha(directory/'protocol.json'),water_fixture_GPU_CPU_difference_kj=fixture_differences,rows=rows,geometry_intervention_kj=contrasts,per_family_aqueous_minus_gas_geometry_kj=per_family,
  parent_DFT_gas_electronic_change_to_geometry_effect_kj=contrasts['DFT_ALPB']-contrasts['UMA_ALPB'],
  parent_SMD_change_to_geometry_effect_at_UMA_kj=contrasts['UMA_SMD']-contrasts['UMA_ALPB'],
  combined_parent_DFT_SMD_minus_UMA_ALPB_change_kj=contrasts['DFT_SMD']-contrasts['UMA_ALPB'],
  source_thermal_geometry_effect_kj=read(directory/'parents/matched_profile_diagnostic.json')['RRHO_profile_shift_kj'],
  scope='Fixed-microgeometry electronic/implicit-solvent energy intervention only. No equilibrium cycle-closure requirement, parent nuclear minima/Hessians, complete conformer ensembles, biochemical standardGibbs or experimental-accuracy validation.',
  factorization_note='SMD increments are shared parent-density corrections in the two electronic arms. Single-geometry linear arithmetic makes the factorial interaction identically zero; it is not evidence of independent physical errors.',
  production_promotion=False,experimental_labels_admitted=0)
 return result

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('directory',type=Path);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
 result=audit(args.directory);args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='rows'},indent=2))
