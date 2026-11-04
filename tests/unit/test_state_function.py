"""State-function invariants with synthetic molecular energies, not accuracy tests."""
import copy
import hashlib
import json
from collections import Counter
from pathlib import Path

import pytest
from rdkit import Chem

from metag.settings import Settings
from metag.state_function import StateFunctionEstimator, plan_state


class Backend:
    model_name = "synthetic-state-function-control"
    producer = {"cacheable": False, "kind": "synthetic"}

    def __init__(self):
        self.settings = Settings()
        self.calls = Counter()
        self.water_calls = 0

    def species(self, q, smiles, label, warnings, diagnostics):
        self.calls[(q, smiles)] += 1
        molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
        # Include nonadditive structure information and a changing repeat term;
        # exact elemental additivity would make cycle checks uninformative.
        h = int(hashlib.sha256(smiles.encode()).hexdigest()[:6], 16) % 97
        value = -100. * sum(a.GetAtomicNum() for a in molecule.GetAtoms()) + h + self.calls[(q, smiles)]
        diagnostics["method"] = "synthetic"
        return value, .5

    def liquid_water(self):
        self.water_calls += 1
        return -1000. + self.water_calls


def ester():
    return {"species": {"ester": [-1, 0, "COC(C)=O"], "water": [-1, 0, "O"],
                        "acetate": [1, -1, "CC(=O)[O-]"], "methanol": [1, 0, "CO"]},
            "n_Hplus": 1}


def reverse(reaction, scale=-1):
    out = copy.deepcopy(reaction)
    out["species"] = {n: [scale*c, q, s] for n, (c, q, s) in out["species"].items()}
    out["n_Hplus"] = scale * out.get("n_Hplus", 0)
    return out


def test_reversal_scaling_and_aliases_use_identical_molecular_records():
    backend = Backend()
    estimator = StateFunctionEstimator(backend)
    rx = ester()
    a = estimator.score_reaction(rx)
    b = estimator.score_reaction(reverse(rx))
    c = estimator.score_reaction(reverse(rx, 2.5))
    renamed = dict(rx, species={"alias_"+n: v for n, v in rx["species"].items()})
    d = estimator.score_reaction(renamed)
    assert a["dG"] is not None, a["suspect"]
    assert b["dG_unrounded"] == pytest.approx(-a["dG_unrounded"], abs=1e-10)
    assert c["dG_unrounded"] == pytest.approx(2.5*a["dG_unrounded"], abs=1e-10)
    assert d["dG_unrounded"] == a["dG_unrounded"]
    assert max(backend.calls.values()) == 1
    assert backend.water_calls == 1


def test_shared_species_cancels_even_with_sampling_noise_and_hydration():
    estimator = StateFunctionEstimator(Backend())
    reaction = {"species": {"left": [-1, -1, "CC(=O)C(=O)[O-]"],
                            "right": [1, -1, "CC(=O)C(=O)[O-]"]}}
    result = estimator.score_reaction(reaction)
    assert result["dG"] == 0, result["suspect"]
    assert result["U_samp"] == 0
    assert result["state_records"]["left"] == result["state_records"]["right"]


def test_proton_drawings_share_reference_and_sampling_uncertainty():
    estimator = StateFunctionEstimator(Backend())
    reaction = {"species": {"acid": [-1, 0, "CC(=O)O"],
                            "anion": [1, -1, "CC(=O)[O-]"]}, "n_Hplus": 1}
    result = estimator.score_reaction(reaction)
    assert result["dG"] == 0, result["suspect"]
    assert result["U_samp"] == 0


def test_records_and_ledger_exports_are_copies():
    estimator = StateFunctionEstimator(Backend())
    a = estimator.state(-1, "CC(=O)[O-]")
    original = copy.deepcopy(a)
    a["potential_kj"] = 0
    a["plan"]["pka_sites"].clear()
    export = estimator.ledger()
    export["config"]["state_function_policy"] = "changed"
    assert estimator.state(-1, "CC(=O)[O-]") == original
    assert estimator.ledger()["config"]["state_function_policy"] != "changed"


def test_molecular_reference_is_independent_of_smiles_drawing():
    assert plan_state(-1, "CC(=O)[O-]") == plan_state(-1, "[O-]C(C)=O")


def test_backend_change_cannot_mix_producers():
    backend = Backend()
    estimator = StateFunctionEstimator(backend)
    estimator.state(0, "CO")
    backend.producer = {"kind": "different"}
    with pytest.raises(ValueError, match="identity/settings changed"):
        estimator.state(0, "CO")


@pytest.mark.parametrize("field,value", [("pka_sites", [["prod", 4., "acid"]]),
                                         ("species_state_pools", {"methanol": {}})])
def test_unowned_or_unsupported_state_declarations_fail_before_molecular_work(field, value):
    backend = Backend()
    reaction = dict(ester(), **{field: value})
    result = StateFunctionEstimator(backend).score_reaction(reaction)
    assert result["dG"] is None and field in result["suspect"]
    assert not backend.calls


def test_unbalanced_input_fails_before_molecular_work():
    backend = Backend()
    reaction = dict(ester(), n_Hplus=0)
    result = StateFunctionEstimator(backend).score_reaction(reaction)
    assert result["dG"] is None and "balanced" in result["suspect"]
    assert not backend.calls


def test_known_amino_acid_cycle_closes_with_one_shared_ledger():
    panel = json.loads((Path(__file__).parents[2]/"src/metag/data/reactions_opentecr_std.json").read_text())
    cycle = {"rxn00182": 1, "rxn00248": -1, "rxn00260": 1, "rxn00347": -1, "rxn00799": 1}
    estimator = StateFunctionEstimator(Backend())
    values = []
    for rid, weight in cycle.items():
        result = estimator.score_reaction(panel[rid], key=rid)
        assert result["dG"] is not None, (rid, result["suspect"])
        values.append(weight * result["dG_unrounded"])
    assert sum(values) == pytest.approx(0, abs=1e-8)


@pytest.mark.parametrize("value,sigma", [(float("nan"), .5), (1., -1.), (True, .5)])
def test_invalid_backend_outputs_never_enter_ledger(value, sigma):
    backend = Backend()
    backend.species = lambda *args: (value, sigma)
    estimator = StateFunctionEstimator(backend)
    result = estimator.score_reaction(ester())
    assert result["dG"] is None
    assert not estimator.ledger()["states"]
