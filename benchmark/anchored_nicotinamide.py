"""Research-only molecule-local nicotinamide transfer approximation.

This module is not a production scoring route.  It tests whether using one
full oxidized reference per *exact* ribosyl scaffold can retain small-core
redox differences without reaction-dependent metabolite energies.  Locality
and the inherited protonation/solvation approximations require validation.
"""
import copy
import json
import math

from rdkit import Chem

from metag.chem import is_balanced
from metag.energetics.provenance import content_hash
from metag.state_function import StateFunctionEstimator, _combine_sampling

POLICY = "research-ribosyl-nicotinamide-oxidized-anchor-v1"
OX_CORE = Chem.MolToSmiles(Chem.MolFromSmiles("NC(=O)c1ccc[n+](C)c1"))
RED_CORE = Chem.MolToSmiles(Chem.MolFromSmiles("NC(=O)C1=CN(C)C=CC1"))
_RED = Chem.MolFromSmarts("[NX3H2][CX3](=O)C1=CN([#6])C=CC1")


def anchor_plan(charge, smiles):
    """Construct an exact oxidized sibling from one N-ribosyl reduced ring.

    All atoms outside the six-membered ring, their bonds, charges and stereo
    tags are retained. NAD and NADP therefore have distinct full references.
    Non-ribosyl molecules, including the two methyl cores, use the full model.
    No reaction, enzyme class, raw energy or experimental label is consulted.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or Chem.GetFormalCharge(mol) != charge:
        raise ValueError("invalid nicotinamide input identity/charge")
    matches = mol.GetSubstructMatches(_RED)
    if not matches:
        return None
    rings = set()
    all_rings = mol.GetRingInfo().AtomRings()
    for match in matches:
        matched = set(match)
        for ring in all_rings:
            if len(ring) != 6 or not set(ring) <= matched:
                continue
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            nitrogens = [a for a in atoms if a.GetAtomicNum() == 7]
            if len(nitrogens) != 1 or any(a.GetAtomicNum() not in (6, 7) for a in atoms):
                continue
            attached = [a for a in nitrogens[0].GetNeighbors() if a.GetIdx() not in ring]
            if len(attached) != 1 or attached[0].GetAtomicNum() != 6:
                continue
            # Scope is the N-ribosyl scaffold, not arbitrary N-alkyl analogues.
            if any(len(sugar) == 5 and attached[0].GetIdx() in sugar
                   and sorted(mol.GetAtomWithIdx(i).GetAtomicNum() for i in sugar) == [6, 6, 6, 6, 8]
                   for sugar in all_rings):
                rings.add(tuple(sorted(ring)))
    if not rings:
        return None
    if len(rings) != 1:
        raise ValueError("multiple reduced nicotinamide rings are outside the anchored policy")
    ring = set(next(iter(rings)))
    edited = Chem.RWMol(mol)
    for i in ring:
        atom = edited.GetAtomWithIdx(i)
        atom.SetIsAromatic(True)
        atom.SetNumExplicitHs(0)
        atom.SetNoImplicit(False)
        if atom.GetAtomicNum() == 7:
            atom.SetFormalCharge(1)
    for bond in edited.GetBonds():
        if bond.GetBeginAtomIdx() in ring and bond.GetEndAtomIdx() in ring:
            bond.SetBondType(Chem.BondType.AROMATIC)
            bond.SetIsAromatic(True)
    oxidized = edited.GetMol()
    Chem.SanitizeMol(oxidized)
    target = Chem.MolToSmiles(mol, isomericSmiles=True)
    reference = Chem.MolToSmiles(oxidized, isomericSmiles=True)
    reference_charge = Chem.GetFormalCharge(oxidized)
    terms = [(1., reference_charge, reference), (1., 0, RED_CORE), (-1., 1, OX_CORE)]
    accounting = {"target": [-1, charge, target]}
    accounting.update({str(i): list(term) for i, term in enumerate(terms)})
    if not is_balanced(accounting, 0):
        raise ValueError("anchored reference/core transfer fails atom or charge accounting")
    # Oxidation may alter only the six ring atoms and their internal bonds.
    for atom in mol.GetAtoms():
        if atom.GetIdx() not in ring:
            other = oxidized.GetAtomWithIdx(atom.GetIdx())
            if (atom.GetAtomicNum(), atom.GetFormalCharge(), atom.GetChiralTag()) != (
                    other.GetAtomicNum(), other.GetFormalCharge(), other.GetChiralTag()):
                raise ValueError("reference conversion changed spectator identity")
    return {"policy": POLICY, "target_charge": int(charge), "target_smiles": target,
            "oxidized_sibling_charge": reference_charge, "oxidized_sibling_smiles": reference,
            "raw_terms": [{"coefficient": w, "charge": q, "smiles": s} for w, q, s in terms],
            "balanced_reference_transfer": True,
            "physical_scope": "Fixed-scaffold locality approximation; no measured offset or interaction correction."}


class AnchoredNicotinamideEstimator(StateFunctionEstimator):
    """One molecular ledger, with exact shared-primitive sampling derivatives."""
    def __init__(self, backend):
        super().__init__(backend)
        self.config.update(state_function_policy=POLICY,
                           physical_validation="development study; not a production policy")
        self._anchor_plans = {}

    def _species_energy(self, charge, smiles, label, warnings, diagnostics):
        plan = anchor_plan(charge, smiles)
        if plan is None:
            return super()._species_energy(charge, smiles, label, warnings, diagnostics)
        key = (int(charge), plan["target_smiles"])
        if key not in self._raw:
            values, primitive_records, details = [], [], []
            for term in plan["raw_terms"]:
                q, smi, weight = term["charge"], term["smiles"], term["coefficient"]
                detail = {}
                value, sigma = super()._species_energy(q, smi, label, warnings, detail)
                values.append(weight * value)
                identity = content_hash({"charge": q, "smiles": smi})
                primitive_records.append((weight, {"sampling_inputs": {identity: {
                    "charge": q, "smiles": smi, "sigma_kj": sigma, "sensitivity": 1.}}}))
                details.append(dict(term=term, raw_diagnostics=detail))
            sigma, inputs = _combine_sampling(primitive_records)
            detail = dict(anchored_transfer=copy.deepcopy(plan), primitive_records=details,
                          primitive_sampling_inputs=inputs)
            self._raw[key] = json.dumps([math.fsum(values), sigma, detail], sort_keys=True, allow_nan=False)
            self._anchor_plans[key] = copy.deepcopy(plan)
        value, sigma, detail = json.loads(self._raw[key])
        diagnostics.update(detail)
        return value, sigma

    def _sampling_inputs(self, charge, smiles, routes):
        inherited = super()._sampling_inputs(charge, smiles, routes)
        expanded = []
        for row in inherited.values():
            key = (row["charge"], row["smiles"])
            if key in self._anchor_plans:
                _, _, detail = json.loads(self._raw[key])
                expanded.append((row["sensitivity"], {
                    "sampling_inputs": detail["primitive_sampling_inputs"]}))
            else:
                identity = content_hash({"charge": row["charge"], "smiles": row["smiles"]})
                expanded.append((1., {"sampling_inputs": {identity: row}}))
        return _combine_sampling(expanded)[1]

    def ledger(self):
        ledger = super().ledger()
        ledger.pop("sha256")
        ledger["policy"] = POLICY
        ledger["anchored_reference_plans"] = [copy.deepcopy(self._anchor_plans[k])
                                             for k in sorted(self._anchor_plans)]
        return dict(ledger, sha256=content_hash(ledger))
