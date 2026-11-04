"""Narrow, experimental alpha-amino-acid transamination routing hypothesis.

Research diagnostic only. This rule has not been validated against physical
energies or molecular protonation data and is not a MetaG production option.
"""
from rdkit import Chem

from metag.routing import pka_transform as pk


KETO_ACID = Chem.MolFromSmarts('[C;X3](=O)[C;X3](=O)[O]')


def candidate(species):
    """Exactly one primary alpha-ammonium and one alpha-keto acid on each side.

    Deliberately excludes cofactors, side-chain nitrogens, amides, phosphorus,
    sulfur, neutral input amines and non-unit stoichiometry. This is a structural
    screen, not a claim that a pKa approximation applies to every match.
    """
    if len(species) != 4:
        return False
    sides = {-1: [], 1: []}
    for coeff, charge, smiles in species.values():
        if coeff not in (-1, 1):
            return False
        mol = Chem.MolFromSmiles(smiles)
        if mol is None or any(a.GetAtomicNum() not in (1, 6, 7, 8) for a in mol.GetAtoms()):
            return False
        n_atoms = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 7]
        aa = mol.GetSubstructMatches(pk._ALPHA_AMINO_ACID)
        keto = mol.GetSubstructMatches(KETO_ACID)
        if len(n_atoms) == 1 and len(aa) == 1 and not keto:
            n = n_atoms[0]
            if n.GetIdx() != aa[0][0] or n.GetFormalCharge() != 1 or n.GetTotalNumHs() != 3:
                return False
            kind = 'amino'
        elif not n_atoms and len(keto) == 1 and not aa:
            kind = 'keto'
        else:
            return False
        sides[int(coeff)].append(kind)
    return all(sorted(side) == ['amino', 'keto'] for side in sides.values())


def wrap(original):
    def routed(species, *args, **kwargs):
        if candidate(species):
            kwargs['force_base'] = True
        return original(species, *args, **kwargs)
    return routed
