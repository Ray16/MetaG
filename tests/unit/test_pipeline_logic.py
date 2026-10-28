"""CPU tests of metag.pipeline.score_reaction bookkeeping with a deterministic fake backend (no GPU):
  * 1 atm -> 1 M standard state: every species (solutes AND liquid water) gets +RT ln 24.46, so a
    reaction's species sum shifts by exactly STD_STATE_KJ·Δn and Δn = 0 reactions are untouched;
  * the truncated-core failure path re-scores the ORIGINAL dict with full molecules (was KeyError);
  * the sampling spread U_samp follows the backend sigma;
  * route provenance is returned."""
import pytest

P = pytest.importorskip("metag.pipeline")
from metag.constants import STD_STATE_KJ
from metag.energetics.errors import SpeciesRearranged

G_FAKE = {"CCO": -100.0, "CC=O": -50.0, "[HH]": -5.0, "CC(=O)O": -300.0, "CCCC": -400.0, "CC": -190.0,
          "CO": -80.0, "C": -30.0, "O=C=O": -350.0, "CCOC(C)=O": -350.0, "CCCCCC": -600.0,
          "CC(O)O": -260.0}
# every fixture reaction is element- and charge-balanced: an unbalanced one (correctly) returns no estimate
G_WATER = -200.0


class FakeBackend:
    model_name = "fake-uma"
    producer = {"cacheable": False, "kind": "test"}

    def __init__(self, energy=None, sigma=0.5):
        self.energy = energy or (lambda q, smi: G_FAKE[smi])
        self.sigma = sigma
        self.calls = []

    def species(self, q, smi, label, warnings, diagnostics):
        self.calls.append((q, smi))
        value = self.energy(q, smi)
        return (None, None) if value is None else (value, self.sigma)

    def liquid_water(self):
        return G_WATER


def _rx(species, n_h=0):
    return {"note": "test", "n_Hplus": n_h, "species": species}


def _score(rx, backend=None):
    return P.score_reaction(rx, backend or FakeBackend(), log=lambda *_: None)


def test_std_state_constant():
    assert abs(STD_STATE_KJ - 7.93) < 0.01


def test_std_state_shifts_by_delta_n():
    # acetic acid -> methane + CO2 (Δn = +1, 1 reactant -> 2 products)
    lyase = _rx({"A": [-1, 0, "CC(=O)O"], "B": [1, 0, "C"], "C": [1, 0, "O=C=O"]})
    # esterification acid + ethanol -> ester + water (Δn = 0 only if liquid water is counted too)
    iso = _rx({"A": [-1, 0, "CC(=O)O"], "B": [-1, 0, "CCO"], "E": [1, 0, "CCOC(C)=O"], "W": [1, 0, "O"]})
    on_l, on_i = _score(lyase)["stages"]["species_sum"], _score(iso)["stages"]["species_sum"]
    assert abs(on_l - (-30 - 350 + 300) - STD_STATE_KJ) < 0.06     # Δn = +1
    assert abs(on_i - (-350 - 200 + 300 + 100)) < 0.06             # Δn = 0 (water counted) -> unchanged


def test_truncation_failure_retries_full_molecules(monkeypatch):
    import metag.routing.truncate as T
    # a (balanced) truncated core whose species QM fails -> the ORIGINAL reaction is re-scored in full
    monkeypatch.setattr(T, "build_truncated_reaction", lambda sp, radius=2, **kw: ({"X_t": [-1, 0, "CCCCCC"],
                                                                              "Y_t": [1, 0, "CCCCCC"]}, 0))
    backend = FakeBackend(lambda q, smi: None if smi == "CCCCCC" else G_FAKE[smi])
    r = _score(_rx({"A": [-1, 0, "CCO"], "H": [-1, 0, "[HH]"], "B": [1, 0, "CO"], "C": [1, 0, "C"]}), backend)
    assert r is not None                                     # was: KeyError('rxn') via run_reaction
    assert r["routes"]["truncated"] is False
    assert abs(r["dG_raw"] - (-80 - 30 + 100 + 5)) < 0.06    # Δn = 0


def test_sampling_spread_follows_backend_sigma():
    rx = _rx({"A": [-1, 0, "CCO"], "B": [1, 0, "CC=O"], "H2": [1, 0, "[HH]"]})
    small = _score(rx, FakeBackend(sigma=0.5))
    big = _score(rx, FakeBackend(sigma=30.0))
    assert big["dG"] == small["dG"]
    assert big["U_samp"] > small["U_samp"] + 20
    assert "sigma_pred" not in big and "ci95" not in big


def test_routes_and_scope_reported():
    r = _score(_rx({"A": [-1, 0, "CCO"], "B": [1, 0, "CC=O"], "H2": [1, 0, "[HH]"]}))
    assert set(r["routes"]) >= {"cofactor_ring", "truncated", "ph0", "errors"}
    assert r["std_state_kJ"] == round(STD_STATE_KJ, 2)


def _rearranging_backend():
    def energy(q, smi):
        raise SpeciesRearranged("A", smi, 3)
    return FakeBackend(energy)


def test_rearrangement_without_zwitterion_is_not_rerouted():
    backend = _rearranging_backend()
    rx = _rx({"A": [-1, 0, "CCO"], "B": [1, 0, "CCO"]})
    result = _score(rx, backend)
    assert result["dG"] is None
    assert len(backend.calls) == 1
    assert result["config"]["force_neutral"] is False
    attempt = result["routes"]["species_rearrangement"]
    assert attempt["neutral_reference_permitted"] is False
    assert attempt["requested_species"] == rx["species"]


def test_rearrangement_does_not_replace_supplied_pka_sites():
    backend = _rearranging_backend()
    rx = _rx({"A": [-1, 0, "CCO"], "B": [1, 0, "CCO"]})
    rx["pka_sites"] = [["react", 15.0]]
    result = _score(rx, backend)
    assert result["dG"] is None and len(backend.calls) == 1
    assert result["routes"]["species_rearrangement"]["requested_pka_sites"] == rx["pka_sites"]


@pytest.mark.parametrize("smi,charge", [("C[NH3+]", 1), ("CC(=O)[O-]", -1)])
def test_rearrangement_cannot_expand_the_declared_reference_policy(smi, charge):
    backend = _rearranging_backend()
    result = _score(_rx({"A": [-1, charge, smi], "B": [1, charge, smi]}), backend)
    assert result['dG'] is None and len(backend.calls) == 1
    assert result['config']['force_neutral'] is False
    assert not result['routes']['species_rearrangement']['neutral_reference_permitted']


def test_protonation_diagnostic_matches_point_correction():
    rx = _rx({"acid": [-1, 0, "CC(=O)O"], "methane": [1, 0, "C"],
              "carbon_dioxide": [1, 0, "O=C=O"]})
    transformed = _score(rx)
    diagnostic = transformed["routes"]["protonation_diagnostics"]
    st = transformed["stages"]
    assert transformed["dG_raw"] - st["species_sum"] - st["proton"] == pytest.approx(
        diagnostic["correction_kj"], abs=.11)
    assert st["pka_transform"] == pytest.approx(diagnostic["correction_kj"], abs=.005)
    assert diagnostic["structural_coverage"]["acid"]["issues"]
    assert "uncertainty" in diagnostic["scope"]
