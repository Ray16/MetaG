"""Rotational symmetry number + linearity detection (metag.symmetry). The high-impact species -- water
(sigma=2, the hydrolysis rotational-entropy term) and the linear CO2/O2/H2/N2 (ubiquitous in ModelSEED) --
MUST be exact. Rare high-symmetry floppy species may undercount conservatively, never overcount."""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
from metag.symmetry import geometry_and_sigma
import numpy as np
import pytest

RDLogger.DisableLog("rdApp.*")

CASES = [
    ("O", 2, "nonlinear"), ("O=C=O", 2, "linear"), ("O=O", 2, "linear"),
    ("[H][H]", 2, "linear"), ("N#N", 2, "linear"), ("[C-]#[O+]", 1, "linear"),
    ("C#N", 1, "linear"), ("N", 3, "nonlinear"), ("C", 12, "nonlinear"),
    ("CC", 6, "nonlinear"), ("CO", 1, "nonlinear"),
    ("OC[C@@H](O)[C@@H](O)[C@H](O)[C@H](O)C=O", 1, "nonlinear"),
]


def _geom(smi):
    m = Chem.AddHs(Chem.MolFromSmiles(smi))
    AllChem.EmbedMolecule(m, randomSeed=3)
    try:
        AllChem.MMFFOptimizeMolecule(m)
    except Exception:
        pass
    c = m.GetConformer()
    return ([a.GetSymbol() for a in m.GetAtoms()],
            [[c.GetAtomPosition(i).x, c.GetAtomPosition(i).y, c.GetAtomPosition(i).z]
             for i in range(m.GetNumAtoms())])


def test_symmetry_and_linearity():
    for smi, want_sigma, want_geo in CASES:
        s, crd = _geom(smi)
        geo, sigma, ndrop = geometry_and_sigma(s, crd)
        assert (sigma, geo) == (want_sigma, want_geo), f"{smi}: got ({geo},{sigma}) want ({want_geo},{want_sigma})"
        assert ndrop == {"monatomic": 3, "linear": 5, "nonlinear": 6}[geo]


@pytest.mark.parametrize('scale', [.1, 1., 10.])
def test_bent_light_geometry_is_not_linear_under_scale_or_rotation(scale):
    angle = np.deg2rad(134.)/2
    xyz = scale*np.array([[0,0,0],[1.08*np.cos(angle),1.08*np.sin(angle),0],
                         [1.08*np.cos(angle),-1.08*np.sin(angle),0]])
    # Symmetry tolerance is in Angstrom; the classifier itself also tests rank.
    assert geometry_and_sigma(['C','H','H'], xyz)[0::2] == ('nonlinear',6)
    if scale == 1.:
        assert geometry_and_sigma(['C','H','H'], xyz) == ('nonlinear',2,6)
    xyz = xyz @ np.array([[0,1,0],[0,0,1],[1,0,0]]) + [4.,-8.,12.]
    assert geometry_and_sigma(['C','H','H'], xyz)[0::2] == ('nonlinear',6)


def test_geometric_rank_resolves_nearly_linear_triatomics():
    from metag.symmetry import is_linear
    assert is_linear(['O','C','O'], [[-1.16,0,0],[0,0,0],[1.16,1e-5,0]])
    assert not is_linear(['O','C','O'], [[-1.16,0,0],[0,0,0],[1.16,.02,0]])


if __name__ == "__main__":
    test_symmetry_and_linearity()
    print(f"PASS: {len(CASES)}/{len(CASES)} symmetry cases exact")
