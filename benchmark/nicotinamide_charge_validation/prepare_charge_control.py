"""Prepare a separate phosphate-charge control; submit no computation."""
from collections import Counter
import datetime
import hashlib
import json
import math
from pathlib import Path
from rdkit import Chem, rdBase

HERE=Path(__file__).resolve().parent
DEST=HERE.parent/'nicotinamide_phosphate_charge_20261104_01'
assert not DEST.exists(), 'Preserve an existing control; inspect it before changing inputs'
parent=json.loads((HERE/'geometry_inputs.json').read_text())
rows={}
for family,slot in [('oxidized','full_ox__minimum_031'),('reduced','full_red__minimum_046')]:
    item=parent['inputs'][slot]
    geometry=item['geometry']
    payload={k:v for k,v in geometry.items() if k!='sha256'}
    assert hashlib.sha256(json.dumps(payload,sort_keys=True,allow_nan=False).encode()).hexdigest()==geometry['sha256']
    molecule=Chem.AddHs(Chem.MolFromSmiles(item['reference_smiles']))
    assert [atom.GetSymbol() for atom in molecule.GetAtoms()]==geometry['symbols']
    assert Chem.GetFormalCharge(molecule)==geometry['charge']
    # Verify atom indexing against actual bonded distances before deleting H+.
    for bond in molecule.GetBonds():
        a,b=bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()
        distance=math.dist(geometry['positions_angstrom'][a],geometry['positions_angstrom'][b])
        assert .6<distance<2.3,(slot,a,b,distance)
    oxygen_hydrogen=[]
    for atom in molecule.GetAtoms():
        if atom.GetSymbol()!='O' or atom.GetFormalCharge()!=0:continue
        phosphorus=[n for n in atom.GetNeighbors() if n.GetSymbol()=='P']
        hydrogens=[n for n in atom.GetNeighbors() if n.GetSymbol()=='H']
        if len(phosphorus)==len(hydrogens)==1 and atom.GetDegree()==2:
            oxygen_hydrogen.append((atom.GetIdx(),hydrogens[0].GetIdx()))
    assert len(oxygen_hydrogen)==2
    removed=sorted(h for o,h in oxygen_hydrogen)
    assert all(h>o for o,h in oxygen_hydrogen)
    revised=Chem.RWMol(molecule)
    for oxygen,hydrogen in oxygen_hydrogen:
        revised.GetAtomWithIdx(oxygen).SetFormalCharge(-1)
        revised.GetAtomWithIdx(oxygen).SetNoImplicit(True)
    for index in sorted(removed,reverse=True):revised.RemoveAtom(index)
    revised=revised.GetMol();Chem.SanitizeMol(revised)
    keep=[i for i in range(len(geometry['symbols'])) if i not in removed]
    new=dict(schema='ordered-geometry-v1',symbols=[geometry['symbols'][i] for i in keep],
             positions_angstrom=[geometry['positions_angstrom'][i] for i in keep],
             charge=geometry['charge']-2,spin_multiplicity=1)
    new['sha256']=hashlib.sha256(json.dumps(new,sort_keys=True,allow_nan=False).encode()).hexdigest()
    assert Chem.GetFormalCharge(revised)==new['charge']
    assert [atom.GetSymbol() for atom in revised.GetAtoms()]==new['symbols']
    assert Counter(new['symbols'])+Counter({'H':2})==Counter(geometry['symbols'])
    electrons=sum(atom.GetAtomicNum() for atom in revised.GetAtoms())-new['charge']
    assert electrons==item['n_electrons']
    rows[family]=dict(parent_slot=slot,parent_geometry_sha256=geometry['sha256'],
        parent_reference_smiles=item['reference_smiles'],
        deprotonated_smiles=Chem.MolToSmiles(Chem.RemoveHs(revised),isomericSmiles=True),
        removed_H_indices=removed,deprotonated_O_indices=[o for o,h in oxygen_hydrogen],
        geometry=new,n_electrons=electrons,
        expected_n_AO=item['n_AO']-18,
        expected_AO_count_note='Two def2-TZVPD H atoms removed; verify in the actual basis build.')
weights={'charged_reduced':1,'charged_oxidized':-1,'reference_reduced':-1,'reference_oxidized':1}
states={}
for family,row in rows.items():
    states['charged_'+family]=row['geometry']
    states['reference_'+family]=parent['inputs'][row['parent_slot']]['geometry']
residual=Counter();charge=0
for key,weight in weights.items():
    for symbol,count in Counter(states[key]['symbols']).items():residual[symbol]+=weight*count
    charge+=weight*states[key]['charge']
assert not any(residual.values()) and charge==0 and sum(weights.values())==0
protocol=dict(schema='metag-nicotinamide-phosphate-charge-control-v1',
    created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
    question='How does phosphate deprotonation change the vertical aqueous electronic redox contrast?',
    endpoint='(charged_reduced-charged_oxidized)-(reference_reduced-reference_oxidized)',
    weights=weights,element_residual={},charge_residual=0,molecule_residual=0,
    geometry_rule='Remove only the two verified P-O-H hydrogens from each previously retained full geometry; preserve every other coordinate and atom order.',
    method='Restricted-singlet wb97m-v/def2-tzvpd, def2-universal-jkfit; SMD water, grids4/3, sasa590, equilibrium solvent.',
    numerical_protocol_source='../nicotinamide_reference_control_20261104_01/fixed_geometry_protocol.json',
    parent_protocol_sha256=hashlib.sha256((HERE/'fixed_geometry_protocol.json').read_bytes()).hexdigest(),
    parent_input_sha256=hashlib.sha256((HERE/'geometry_inputs.json').read_bytes()).hexdigest(),
    parent_reference_points={slot:hashlib.sha256((HERE/'points'/(slot+'.json')).read_bytes()).hexdigest() for slot in [r['parent_slot'] for r in rows.values()]},
    intended_resource_ceiling='One separately reserved GPU, two serial SMD points, thirty minutes; no submission made by this preparer.',
    physical_geometry_attempts=0,experimental_labels_admitted=0,offsets_fitted=0,
    limitations=['Vertical electronic-plus-solvent coupling only; no new nuclear optimization, RRHO or conformer ensembles.',
                 'A nonzero vertical coupling does not measure a macroscopic pKa or establish the full free-energy correction.',
                 'These two charged states do not establish complete physiological protonation/Mg speciation.',
                 'No production potential or experimental accuracy is selected from this control.'],
    preparation=dict(rdkit_version=rdBase.rdkitVersion,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()))
DEST.mkdir()
(DEST/'geometry_inputs.json').write_text(json.dumps(dict(states=rows),indent=2)+'\n')
(DEST/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
print(json.dumps(dict(directory=str(DEST),charges={k:v['geometry']['charge'] for k,v in rows.items()},
                     atom_counts={k:len(v['geometry']['symbols']) for k,v in rows.items()},
                     atom_and_charge_balanced=True,calculations_submitted=0),indent=2))
