"""Freeze six cofactor poses and a non-cubic rigid-rotation diagnostic."""
from collections import Counter
from pathlib import Path
import datetime
import hashlib
import json
import shutil
import numpy as np
from rdkit import Chem
from ase import Atoms

HERE=Path(__file__).resolve().parent
BASE=HERE.parent/'nicotinamide_ensemble_control_20261104_01'
CHARGED=HERE.parent/'nicotinamide_phosphate_charge_20261104_01'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(name,value):
    (HERE/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
assert not (HERE/'protocol.json').exists()
source=json.loads((BASE/'geometry_inputs.json').read_text())
checks=[]
def check_attachment(name,geometry,template):
    symbols=geometry['symbols'];positions=np.array(geometry['positions_angstrom'])
    assert symbols==[atom.GetSymbol() for atom in template.GetAtoms()],name
    heavy=[atom.GetIdx() for atom in template.GetAtoms() if atom.GetAtomicNum()!=1]
    margins=[];lengths=[]
    for atom in template.GetAtoms():
        if atom.GetAtomicNum()!=1:continue
        i=atom.GetIdx();expected=[n.GetIdx() for n in atom.GetNeighbors()]
        assert len(expected)==1 and expected[0] in heavy
        distances={j:float(np.linalg.norm(positions[i]-positions[j])) for j in heavy}
        nearest=min(distances,key=distances.get)
        assert nearest==expected[0],(name,i,expected,nearest,distances[nearest])
        other=min(v for j,v in distances.items() if j!=expected[0])
        margins.append(other-distances[expected[0]]);lengths.append(distances[expected[0]])
    assert Chem.GetFormalCharge(template)==geometry['charge']
    return dict(name=name,n_H=len(margins),minimum_H_attachment_margin_A=min(margins),
        maximum_H_bond_length_A=max(lengths),charge=geometry['charge'])
for key,item in source['inputs'].items():
    checks.append(check_attachment(key,item['geometry'],Chem.AddHs(Chem.MolFromSmiles(item['reference_smiles']))))
write('retained_reference_identity.json',dict(n_retained_geometries=len(checks),checks=checks,
    scope='Independent nearest-heavy-atom attachment check on all103retained geometries. It does not establish complete protonation ensembles, pKa accuracy or Gibbs accuracy.'))

poses={}
for family in ['full_ox','full_red','core_ox','core_red']:
    key=min((k for k,v in source['inputs'].items() if v['family']==family),
        key=lambda k:(source['inputs'][k]['original_components_kj']['total_kj'],source['inputs'][k]['geometry']['sha256']))
    item=source['inputs'][key]
    poses[family]=dict(geometry=item['geometry'],source_slot=key,reference_smiles=item['reference_smiles'],
        original_solvation_kj=item['original_components_kj']['solvation_kj'])
charged=json.loads((CHARGED/'geometry_inputs.json').read_text())['states']
for state,record in charged.items():
    mol=Chem.RWMol(Chem.AddHs(Chem.MolFromSmiles(record['parent_reference_smiles'])))
    for oxygen in record['deprotonated_O_indices']:mol.GetAtomWithIdx(oxygen).SetFormalCharge(-1)
    for i in sorted(record['removed_H_indices'],reverse=True):mol.RemoveAtom(i)
    template=mol.GetMol();Chem.SanitizeMol(template)
    name='charged_'+('ox' if state=='oxidized' else 'red')
    checks.append(check_attachment(name,record['geometry'],template))
    poses[name]=dict(geometry=record['geometry'],source_slot=record['parent_slot'],
        reference_smiles=record['deprotonated_smiles'],original_solvation_kj=None)
assert len(poses)==6
write('charge_reference_identity.json',dict(checks=checks[-2:]))
rotations=[dict(name='identity',matrix=np.eye(3).tolist(),determinant=1)]
rng=np.random.default_rng(20261104)
for i in range(16):
    q=rng.normal(size=4);q/=np.linalg.norm(q);w,x,y,z=q
    matrix=np.array([[1-2*(y*y+z*z),2*(x*y-z*w),2*(x*z+y*w)],
                     [2*(x*y+z*w),1-2*(x*x+z*z),2*(y*z-x*w)],
                     [2*(x*z-y*w),2*(y*z+x*w),1-2*(x*x+y*y)]])
    assert np.max(np.abs(matrix@matrix.T-np.eye(3)))<1e-14 and abs(np.linalg.det(matrix)-1)<1e-14
    rotations.append(dict(name=f'proper_{i:02d}',matrix=matrix.tolist(),determinant=1))
rotations.append(dict(name='inversion_control',matrix=(-np.eye(3)).tolist(),determinant=-1))
tasks=[]
for name,pose in poses.items():
    xyz=np.array(pose['geometry']['positions_angstrom']);center=xyz.mean(axis=0)
    for rotation in rotations:
        matrix=np.array(rotation['matrix'])
        transformed=xyz.copy() if rotation['name']=='identity' else (xyz-center)@matrix.T+center
        assert np.max(np.abs(np.linalg.norm(xyz[:,None]-xyz[None,:],axis=-1)-np.linalg.norm(transformed[:,None]-transformed[None,:],axis=-1)))<1e-12
        for grid in ['default','extreme']:
            for accuracy in ['default','0.001']:
                tasks.append(dict(key=f'{name}__{rotation["name"]}__{grid}__{accuracy}',pose=name,
                    transform=rotation['name'],determinant=rotation['determinant'],grid=grid,accuracy=accuracy,
                    symbols=pose['geometry']['symbols'],positions_angstrom=transformed.tolist(),charge=pose['geometry']['charge']))
assert len(tasks)==432
xtb='/nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/xtb/bin/xtb'
protocol=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
    scope='Fixed-geometry numerical invariance/control study; no experimental reference, prediction-error or production-policy selection.',
    n_poses=6,n_transforms=18,n_proper_orientations=17,n_tasks=432,n_expected_native_executions=864,
    rotation_rule='Identity plus16seeded Gaussian-unit-quaternion proper rotations, not only the cubic symmetry operations of a Lebedev grid. One separate coordinate-inversion control.',
    primary_endpoint='Default-production ALPB solvation range under the17proper orientations, per pose and balanced cofactor exchange.',
    comparison_arms='Default/5810point SASA grid crossed with default/.001SCF accuracy; retain allfourarms.',
    rotational_invariance_tolerance_kj=.05,native_replay_tolerance_kj=.005,
    xtb=xtb,xtb_sha256=sha(xtb),Hartree_to_kj_per_mol=2625.499639,
    coordinates='Production six-decimal XYZ serialization after each rigid transform.',
    restart_policy='Fresh isolated directory per gas/water pair. Water may read the immediately preceding gas restart, matching the production dgsolv call sequence.',
    native_flags='--gfn 2 --chrg {q} --uhf 0 --sp; the production baseline has no explicit accuracy/grid/standard-state argument. Nondefault arms add only --acc 0.001 and/or $gbsa grid=extreme.',
    fixed_pose_exchange_weights=dict(full_red=1,full_ox=-1,core_red=-1,core_ox=1),
    vertical_charge_cycle_weights=dict(charged_red=1,charged_ox=-1,full_red=-1,full_ox=1),
    source_files_sha256={str(BASE/'geometry_inputs.json'):sha(BASE/'geometry_inputs.json'),str(CHARGED/'geometry_inputs.json'):sha(CHARGED/'geometry_inputs.json')},
    native_threads=1,worker_CPUs=[60,61],controller_CPU=79,wall_seconds=3600,per_native_timeout_seconds=180,
    attempts_per_pair=1,consecutive_failure_limit=3,
    interpretation='Orientation ranges are conditional numerical diagnostics, not total-error intervals. A smaller range or denser grid alone does not establish more accurate aqueous Gibbs energies. No global reaction reassembly or method adoption follows automatically.')
write('poses.json',poses);write('transforms.json',rotations);write('tasks.json',tasks);write('protocol.json',protocol)
for directory in ['results','attempts','claims','workers']:(HERE/directory).mkdir()
print(json.dumps(dict(retained_geometries_identity_checked=103,extra_charge_poses_checked=2,poses=6,tasks=432,native_executions=864)))
