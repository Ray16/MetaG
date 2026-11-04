"""score_reaction with a deterministic CPU fake backend (no UMA, no xtb)."""
import hashlib
import math

import pytest

from metag.pipeline import effective_config, no_estimate, score_reaction
from metag.settings import Settings

AMP_HYDROLYSIS = {
    "n_Hplus": 0, "note": "AMP + H2O -> adenosine + Pi",
    "species": {
        "H2O": [-1, 0, "O"],
        "AMP": [-1, -2, "Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)([O-])[O-])[C@@H](O)[C@H]1O"],
        "Phosphate": [1, -2, "O=P([O-])([O-])O"],
        "Adenosine": [1, 0, "Nc1ncnc2c1ncn2[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O"],
    },
}
ESTER_HYDROLYSIS = {
    "n_Hplus": 0, "note": "methyl acetate + H2O -> acetic acid + methanol",
    "species": {"H2O": [-1, 0, "O"], "MeOAc": [-1, 0, "COC(C)=O"],
                "AcOH": [1, 0, "CC(=O)O"], "MeOH": [1, 0, "CO"]},
}


class FakeBackend:
    """Deterministic energies keyed by (charge, SMILES)."""
    model_name = "fake-model"

    def __init__(self, fail=None):   # fail = SMILES whose energy is unavailable
        self.settings = Settings()
        self.producer = {"cacheable": False, "kind": "test-fake"}
        self.fail = fail or set()
        self.calls = []

    def species(self, q, smiles, label, warnings, diagnostics):
        self.calls.append((q, smiles, label))
        if smiles in self.fail:
            return None, None
        h = int(hashlib.sha256(f"{q}|{smiles}".encode()).hexdigest()[:8], 16)
        return -1000.0 - (h % 100000) / 10.0, 0.5

    def liquid_water(self):
        return -1000.0


@pytest.mark.parametrize("reaction", [AMP_HYDROLYSIS, ESTER_HYDROLYSIS], ids=["amp", "ester"])
def test_fake_backend_gives_finite_dG_whose_stages_sum(reaction):
    r = score_reaction(reaction, FakeBackend(), key="t", log=lambda *_: None)
    assert r["dG"] is not None and math.isfinite(r["dG"]), r.get("suspect")
    st = r["stages"]
    assert abs(st["species_sum"] + st["proton"] + st["pka_transform"] - r["dG"]) <= 0.06
    assert r["dG"] == r["dG_raw"]
    assert r["config"]["model"] == "fake-model" and r["config"]["solv_model"] == "alpb"


def test_score_is_deterministic():
    a = score_reaction(ESTER_HYDROLYSIS, FakeBackend(), key="t", log=lambda *_: None)
    b = score_reaction(ESTER_HYDROLYSIS, FakeBackend(), key="t", log=lambda *_: None)
    assert a["dG"] == b["dG"] and a["stages"] == b["stages"]


@pytest.mark.parametrize("allow_truncate", [True, False])
def test_failed_species_gives_no_estimate(allow_truncate):
    backend = FakeBackend(fail={"COC(C)=O"})
    r = score_reaction(ESTER_HYDROLYSIS, backend, key="t", log=lambda *_: None, allow_truncate=allow_truncate)
    assert r["dG"] is None and r["dG_raw"] is None
    assert "MeOAc" in r["suspect"] and r["routes"]["failed_species"]


@pytest.mark.parametrize("energy,sigma", [(float('nan'), 0.1), (float('inf'), 0.1),
                                         (1.0, None), (1.0, -1.0), (1.0, float('nan')),
                                         (True, 0.1), (1.0, True)])
def test_invalid_backend_species_values_fail_closed(energy, sigma):
    backend = FakeBackend()
    backend.species = lambda *args: (energy, sigma)
    result = score_reaction(ESTER_HYDROLYSIS, backend, key="t", log=lambda *_: None,
                            allow_truncate=False)
    assert result['dG'] is None
    assert result['routes']['failed_species']
    assert any('invalid species energy or sampling diagnostic' in w for w in result['routes']['warnings'])


@pytest.mark.parametrize("value", [None, float('nan'), float('inf'), True])
def test_invalid_water_reference_fails_closed(value):
    backend = FakeBackend()
    backend.liquid_water = lambda: value
    result = score_reaction(ESTER_HYDROLYSIS, backend, key="t", log=lambda *_: None,
                            allow_truncate=False)
    assert result['dG'] is None
    assert 'H2O' in result['routes']['failed_species']


def test_explicit_water_request_raises():
    with pytest.raises(ValueError):
        score_reaction(dict(ESTER_HYDROLYSIS, explicit=True), FakeBackend(), key="t", log=lambda *_: None)


def test_missing_backend_and_model_is_a_type_error():
    with pytest.raises(TypeError):
        score_reaction(ESTER_HYDROLYSIS, log=lambda *_: None)


@pytest.mark.parametrize('bad_reference', [None, '1.2', [True], [float('nan')], [float('inf')], ['1.2']])
def test_invalid_reference_metadata_is_rejected_before_backend_work(bad_reference):
    backend = FakeBackend()
    with pytest.raises(ValueError, match='reaction exp'):
        score_reaction(dict(ESTER_HYDROLYSIS, exp=bad_reference), backend, log=lambda *_: None)
    assert backend.calls == []


def test_no_estimate_and_effective_config():
    r = no_estimate("k", ESTER_HYDROLYSIS, "boom", log=lambda *_: None)
    assert r["dG"] is None and r["suspect"] == "boom" and r["reaction"] == "k"
    cfg = effective_config(FakeBackend(), allow_truncate=False)
    assert cfg["allow_truncate"] is False and cfg["settings"] == Settings().as_dict()
