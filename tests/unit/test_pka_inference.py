"""Regression checks for pKa domain, provenance, and thermodynamic bookkeeping."""
import json
import math

import pytest
from rdkit import Chem

from metag.energetics.protonation import reaction_transform, site_free_energy
from metag.reactions import validate_reaction
from metag.routing import pka_transform as pka


@pytest.mark.parametrize("neutral,anion", [
    ("CCS", "CC[S-]"), ("Oc1ccccc1", "[O-]c1ccccc1"),
    ("CC(=O)O", "CC(=O)[O-]"), ("OP(=O)(O)O", "[O-]P(=O)([O-])[O-]"),
])
def test_acid_reference_and_transform_independent_of_input_charge(monkeypatch, neutral, anion):
    left, right = pka._neutralize(neutral), pka._neutralize(anion)
    assert left[0] == right[0]
    assert sorted(left[1]) == sorted(right[1])
    assert left[2] == right[2] == 0


@pytest.mark.parametrize("smiles", ["CCO", "CC(=O)S", "COc1ccccc1", "CSSC"])
def test_weak_acid_extension_leaves_other_groups_alone(smiles):
    assert pka._canonicalize_maxanion(smiles) == smiles


def test_monoprotic_equal_population_limit_and_mirror():
    rt_ln10 = 5.7
    assert site_free_energy(7, "acid", 7, rt_ln10) == pytest.approx(-rt_ln10 * math.log10(2))
    assert site_free_energy(4, "acid", 7, rt_ln10) == site_free_energy(10, "base", 7, rt_ln10)
    # Population-dominance limit must remain finite far beyond 10**308.
    assert site_free_energy(-1000, "acid", 7, rt_ln10) == pytest.approx(-rt_ln10 * 1007)
    assert site_free_energy(1000, "acid", 7, rt_ln10) == 0


def test_reversal_scaling_and_spectator_cancellation():
    sites = [["react", 4, "acid", 0.5], ["prod", 9, "base", 2]]
    value = reaction_transform(sites, 7, 5.7)
    assert reaction_transform([["prod" if s == "react" else "react", p, k, m]
                               for s, p, k, m in sites], 7, 5.7) == -value
    assert reaction_transform([[s, p, k, 3*m] for s, p, k, m in sites], 7, 5.7) == pytest.approx(3*value)
    assert reaction_transform([["react", 4], ["prod", 4]], 7, 5.7) == 0


@pytest.mark.parametrize("site", [[], ["wrong", 4], ["react", 4, "wrong"],
    ["prod", float("nan")], ["react", 4, "acid", -1], ["react", 4, "acid", 0],
    ["react", 4, "acid", True], ["react", 4, "acid", 1, "extra"]])
def test_bad_user_terms_rejected_before_scoring(site):
    with pytest.raises(ValueError):
        validate_reaction({"species": {"ethanol": [-1, 0, "CCO"], "ether": [1, 0, "COC"]},
                           "pka_sites": [site]})
