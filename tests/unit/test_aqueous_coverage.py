"""Incomplete state enumeration must not masquerade as an equilibrium estimate."""
from unittest.mock import Mock

import pytest

from metag.energetics.aqueous_states import hydrate_all_sites
from metag.routing.aldehyde_hydration import n_hydration_sites


class _Backend:
    model_name = "fake-uma"
    producer = {"cacheable": False}

    def __init__(self, species_value, water_value):
        self.species = Mock(return_value=species_value)
        self.liquid_water = Mock(return_value=water_value)


@pytest.mark.parametrize("smiles,count,computed", [
    ("O=CCC(=O)CC=O", 3, 5),  # seven subsets, five distinct molecular hydrates
    ("O=CCC(=O)CC(=O)CC=O", 4, 0),
])
def test_hydration_limit_requires_complete_ensemble(smiles, count, computed):
    assert n_hydration_sites(smiles) == count
    energies, sigmas = {"x": -20.0}, {"x": 0.1}
    routes = {"warnings": []}
    backend = _Backend((-25.0, 0.1), -10.0)
    hydrate_all_sites(backend, {}, "x", 0, smiles, energies, sigmas, 0.0, lambda _: None, routes,
                      acid_transform=lambda _: 0.0)
    assert backend.species.call_count == computed
    assert bool(routes.get("incomplete_hydration")) == (count > 3)
    if count > 3:
        assert energies["x"] == -20.0
        backend.liquid_water.assert_not_called()


def test_hydration_limit_fails_closed_at_reaction_level():
    from metag import pipeline as p
    reaction = {"species": {
        "four-carbonyl": [-1, 0, "O=CCC(=O)CC(=O)CC=O"],
        "hydrogen": [-1, 0, "[H][H]"],
        "three-carbonyl": [1, 0, "OCCC(=O)CC(=O)CC=O"],
    }}
    result = p.score_reaction(reaction, _Backend((-50.0, 0.1), -10.0), log=lambda _: None)
    assert result["dG"] is None
    assert result["routes"]["input_balanced"]
    assert "incomplete hydration" in result["suspect"]
