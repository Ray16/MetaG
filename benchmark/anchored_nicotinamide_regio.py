"""Research-only common-reference transfer for two dihydronicotinamide isomers.

The 1,2 and 1,4 forms receive their own methyl-core calculations and the same
exact oxidized full-scaffold reference. No reaction-dependent choice or energy
offset is used. This specifies an approximation, not its physical validation.
"""
from pathlib import Path
import copy
import hashlib
import json
import math

from rdkit import Chem

from benchmark import anchored_nicotinamide as original
from metag.chem import is_balanced
from metag.energetics.provenance import content_hash
from metag.state_function import StateFunctionEstimator, _combine_sampling

POLICY='research-ribosyl-nicotinamide-two-isomer-anchor-v2'
OX_CORE=original.OX_CORE
CORE_14=original.RED_CORE
CORE_12=Chem.MolToSmiles(Chem.MolFromSmiles('NC(=O)C1=CC=CN(C)C1'))
SUPPORTED={CORE_14:'1,4-dihydro',CORE_12:'1,2-dihydro'}


def anchor_plan(charge,smiles):
    mol=Chem.MolFromSmiles(smiles)
    if mol is None or Chem.GetFormalCharge(mol)!=charge:
        raise ValueError('Invalid molecule or declared charge')
    if any(a.GetIsotope() for a in mol.GetAtoms()):
        raise ValueError('Isotope-specific core transfer is unsupported')
    candidates=[]
    for ring_tuple in mol.GetRingInfo().AtomRings():
        ring=set(ring_tuple)
        if len(ring)!=6:continue
        atoms=[mol.GetAtomWithIdx(i) for i in ring]
        nitrogen=[a for a in atoms if a.GetAtomicNum()==7]
        if len(nitrogen)!=1 or any(a.GetAtomicNum() not in {6,7} or a.GetIsAromatic() or a.GetFormalCharge() for a in atoms):continue
        bonds=[b for b in mol.GetBonds() if b.GetBeginAtomIdx() in ring and b.GetEndAtomIdx() in ring]
        if len(bonds)!=6 or sum(b.GetBondType()==Chem.BondType.DOUBLE for b in bonds)!=2:continue
        if any(b.GetBondType() not in {Chem.BondType.SINGLE,Chem.BondType.DOUBLE} for b in bonds):continue
        n=nitrogen[0]
        if any(b.GetBondType()!=Chem.BondType.SINGLE for b in n.GetBonds()):continue
        attached=[a for a in n.GetNeighbors() if a.GetIdx() not in ring]
        if len(attached)!=1 or attached[0].GetAtomicNum()!=6:continue
        if not any(len(s)==5 and attached[0].GetIdx() in s and
            sorted(mol.GetAtomWithIdx(i).GetAtomicNum() for i in s)==[6,6,6,6,8]
            for s in mol.GetRingInfo().AtomRings()):continue
        primary=[]
        for a in atoms:
            for c in a.GetNeighbors():
                if c.GetIdx() in ring or c.GetAtomicNum()!=6 or c.GetDegree()!=3:continue
                oxygen=[x for x in c.GetNeighbors() if x.GetAtomicNum()==8 and
                    mol.GetBondBetweenAtoms(c.GetIdx(),x.GetIdx()).GetBondType()==Chem.BondType.DOUBLE]
                amide=[x for x in c.GetNeighbors() if x.GetAtomicNum()==7 and x.GetDegree()==1
                       and x.GetTotalNumHs()==2 and x.GetFormalCharge()==0]
                if len(oxygen)==len(amide)==1:
                    primary.append({c.GetIdx(),oxygen[0].GetIdx(),amide[0].GetIdx()})
        if len(primary)!=1:continue
        keep=ring|primary[0]
        external=[b for b in mol.GetBonds() if (b.GetBeginAtomIdx() in keep)!=(b.GetEndAtomIdx() in keep)]
        if len(external)!=1 or n.GetIdx() not in {external[0].GetBeginAtomIdx(),external[0].GetEndAtomIdx()}:continue
        core=Chem.RWMol();mapping={}
        for i in sorted(keep):mapping[i]=core.AddAtom(Chem.Atom(mol.GetAtomWithIdx(i)))
        for b in mol.GetBonds():
            if b.GetBeginAtomIdx() in keep and b.GetEndAtomIdx() in keep:
                core.AddBond(mapping[b.GetBeginAtomIdx()],mapping[b.GetEndAtomIdx()],b.GetBondType())
        cap=core.AddAtom(Chem.Atom(6));core.AddBond(mapping[n.GetIdx()],cap,Chem.BondType.SINGLE)
        reduced=core.GetMol();Chem.SanitizeMol(reduced)
        reduced_smiles=Chem.MolToSmiles(reduced,isomericSmiles=True)
        # Oxidize the full molecule without deleting or altering its scaffold.
        oxidized=Chem.RWMol(mol)
        for i in ring:
            a=oxidized.GetAtomWithIdx(i);a.SetIsAromatic(True);a.SetNumExplicitHs(0);a.SetNoImplicit(False)
            if a.GetAtomicNum()==7:a.SetFormalCharge(1)
        for b in oxidized.GetBonds():
            if b.GetBeginAtomIdx() in ring and b.GetEndAtomIdx() in ring:
                b.SetBondType(Chem.BondType.AROMATIC);b.SetIsAromatic(True)
        oxidized=oxidized.GetMol();Chem.SanitizeMol(oxidized)
        # Verify the exact methyl-oxidized sibling is the declared nicotinamide
        # core, excluding positional carboxamides outside this family.
        oc=Chem.RWMol(reduced)
        core_ring=set(next(r for r in reduced.GetRingInfo().AtomRings() if len(r)==6))
        for i in core_ring:
            a=oc.GetAtomWithIdx(i);a.SetIsAromatic(True);a.SetNumExplicitHs(0);a.SetNoImplicit(False)
            if a.GetAtomicNum()==7:a.SetFormalCharge(1)
        for b in oc.GetBonds():
            if b.GetBeginAtomIdx() in core_ring and b.GetEndAtomIdx() in core_ring:
                b.SetBondType(Chem.BondType.AROMATIC);b.SetIsAromatic(True)
        om=oc.GetMol();Chem.SanitizeMol(om)
        if Chem.MolToSmiles(om,isomericSmiles=True)!=OX_CORE:continue
        if reduced_smiles not in SUPPORTED:
            raise ValueError('Reduced nicotinamide isomer has no declared core-transfer policy')
        reference=Chem.MolToSmiles(oxidized,isomericSmiles=True)
        q=Chem.GetFormalCharge(oxidized)
        terms=[dict(coefficient=1.,charge=q,smiles=reference),
               dict(coefficient=1.,charge=0,smiles=reduced_smiles),
               dict(coefficient=-1.,charge=1,smiles=OX_CORE)]
        bookkeeping={'target':[-1,charge,smiles]}
        bookkeeping.update({str(i):[t['coefficient'],t['charge'],t['smiles']] for i,t in enumerate(terms)})
        if not is_balanced(bookkeeping,0):raise ValueError('Core-transfer elements/charge do not balance')
        candidates.append(dict(policy=POLICY,core_isomer=SUPPORTED[reduced_smiles],
            target_charge=charge,target_smiles=Chem.MolToSmiles(mol,isomericSmiles=True),
            oxidized_sibling_charge=q,oxidized_sibling_smiles=reference,
            raw_terms=terms,balanced_reference_transfer=True,
            physical_scope='Common full-scaffold reference and isomer-specific methyl core; locality, speciation and scaffold interactions remain unvalidated.'))
    if len(candidates)>1:raise ValueError('Multiple reduced nicotinamide rings are unsupported')
    return candidates[0] if candidates else None


class RegioAnchoredEstimator(original.AnchoredNicotinamideEstimator):
    def __init__(self,backend):
        super().__init__(backend)
        self.config.update(state_function_policy=POLICY,
            core_recipe_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            inherited_sampling_source_sha256=hashlib.sha256(Path(original.__file__).read_bytes()).hexdigest())

    def _species_energy(self,charge,smiles,label,warnings,diagnostics):
        plan=anchor_plan(charge,smiles)
        if plan is None:return StateFunctionEstimator._species_energy(self,charge,smiles,label,warnings,diagnostics)
        key=(int(charge),plan['target_smiles'])
        if key not in self._raw:
            values=[];records=[];details=[]
            for term in plan['raw_terms']:
                q,smi,w=term['charge'],term['smiles'],term['coefficient'];detail={}
                value,sigma=StateFunctionEstimator._species_energy(self,q,smi,label,warnings,detail)
                values.append(w*value);identity=content_hash(dict(charge=q,smiles=smi))
                records.append((w,dict(sampling_inputs={identity:dict(charge=q,smiles=smi,sigma_kj=sigma,sensitivity=1.)})))
                details.append(dict(term=term,raw_diagnostics=detail))
            sigma,inputs=_combine_sampling(records)
            detail=dict(anchored_transfer=copy.deepcopy(plan),primitive_records=details,primitive_sampling_inputs=inputs)
            self._raw[key]=json.dumps([math.fsum(values),sigma,detail],sort_keys=True,allow_nan=False)
            self._anchor_plans[key]=copy.deepcopy(plan)
        value,sigma,detail=json.loads(self._raw[key]);diagnostics.update(detail)
        return value,sigma

    def ledger(self):
        result=StateFunctionEstimator.ledger(self);result.pop('sha256')
        result.update(policy=POLICY,anchored_reference_plans=[copy.deepcopy(self._anchor_plans[k]) for k in sorted(self._anchor_plans)])
        return dict(result,sha256=content_hash(result))
