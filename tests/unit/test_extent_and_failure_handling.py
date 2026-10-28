"""Extent invariance, optional inputs and failure handling (mocked QM)."""
import json

import pytest

from rdkit import Chem

from metag import pipeline as p
from metag.reactions import load_reactions, stoichiometric_extent
from metag.routing.planner import route_reaction


def equation(scale=1):
    return {"species": {"ethanol": [-scale, 0, "CCO"], "ether": [scale, 0, "COC"]}}


class Backend:
    """Deterministic species energies keyed by canonical SMILES (no QM)."""
    model_name = "fake"
    producer = {"cacheable": False}

    def __init__(self, energy, water=-10.0):
        self.energy, self.water = energy, water

    def species(self, q, smi, label, warnings, diagnostics):
        return self.energy(Chem.MolToSmiles(Chem.MolFromSmiles(smi)))

    def liquid_water(self):
        return self.water


def test_benchmark_equations_have_unit_extent():
    assert all(stoichiometric_extent(rx["species"]) == 1 for rx in load_reactions().values())


@pytest.mark.parametrize("scale", [0.5, 2, 10])
def test_routing_and_point_estimate_are_extent_invariant(monkeypatch, scale):
    # Without normalization, the doubled equation bypasses a core builder
    # restricted to unit stoichiometry and gets a different species energy.
    from metag.routing import truncate
    def build(species, **kw):
        if any(abs(v[0]) != 1 for v in species.values()):
            return None
        return {"a": [-1, 0, "CCCO"], "b": [1, 0, "CCOC"]}, 0
    monkeypatch.setattr(truncate, "build_truncated_reaction", build)
    energies = {"CCO": -10, "COC": -8, "CCCO": -30, "CCOC": -25}
    backend = Backend(lambda smi: (energies[smi], 0.5))
    base = p.score_reaction(equation(), backend, log=lambda *_: None)
    changed = p.score_reaction(equation(scale), backend, log=lambda *_: None)
    assert changed["routes"]["truncated"]
    assert changed["dG"] == pytest.approx(scale * base["dG"], abs=0.11)
    assert changed["U_samp"] == pytest.approx(scale * base["U_samp"], abs=0.6)


def test_routing_preserves_supplied_pka_site_extent():
    rx = equation(2)
    rx["pka_sites"] = [["react", 4.5], ["prod", 6.0, "base", 2]]
    routed, _, _ = route_reaction(rx, allow_truncate=False, log=lambda *_: None)
    assert routed["pka_sites"] == [["react", 4.5, "acid", 1.0], ["prod", 6.0, "base", 2.0]]


def test_score_without_optional_note():
    r = p.score_reaction(equation(), Backend(lambda smi: (-50.0, 0.2)), log=lambda *_: None,
                         allow_truncate=False)
    assert r["dG"] == 0 and r["note"] == ""


def test_cli_records_all_failures_and_keeps_completed_work(monkeypatch, tmp_path):
    from metag.cli import run_batch
    reactions = {k: equation() for k in ("ok", "invalid", "exception")}
    good = {"reaction": "ok", "dG": 1.0, "dG_raw": 1.0, "U_samp": 0.0, "exp": None, "err": []}
    def run(rx, backend, *, key, log):
        if key == "ok":
            return good
        # Results must already be persisted before the next computation begins.
        assert json.loads((tmp_path / "metag_all.json").read_text())[0] == good
        if key == "exception":
            raise ValueError("unsupported species")
        return p.no_estimate(key, equation(), "unbalanced", lambda *_: None)
    monkeypatch.setattr(p, "score_reaction", run)
    run_batch(reactions, list(reactions), object(), str(tmp_path), log=lambda *_: None)
    records = json.loads((tmp_path / "metag_all.json").read_text())
    assert [r["reaction"] for r in records] == ["ok", "invalid", "exception"]
    assert all(r["dG"] is None for r in records[1:])
    assert "unsupported species" in records[2]["suspect"]


def test_failed_hydration_does_not_return_a_point():
    rx = {"species": {"ethanol": [-1, 0, "CCO"], "acetaldehyde": [1, 0, "CC=O"],
                       "hydrogen": [1, 0, "[H][H]"]}}
    backend = Backend(lambda smi: (None, None) if smi == "CC(O)O" else (-50.0, 0.1))
    result = p.score_reaction(rx, backend, log=lambda *_: None)
    assert result["dG"] is None
    assert "incomplete hydration" in result["suspect"]


@pytest.mark.parametrize("energy,sigma", [(float("nan"), 1), (float("inf"), 1), (-10, -1), (-10, float("nan"))])
def test_cache_rejects_invalid_values(monkeypatch, tmp_path, energy, sigma):
    from metag.energetics import species_cache as c
    monkeypatch.setattr(c, "CACHE_DIR", str(tmp_path))
    monkeypatch.setattr(c, "_ENABLED", True)
    c.put("CCO", 0, "implicit", {}, energy, sigma)
    assert c.get("CCO", 0, "implicit", {}) is None
    # Also reject bad files produced by previous versions or damaged outside this API.
    path, smi = c._path("CCO", 0, "implicit", {})
    from pathlib import Path
    Path(path).write_text(json.dumps({"smi": smi, "q": 0, "method": "implicit", "settings": {},
                                    "ver": c.CACHE_VERSION, "G": energy, "sigma": sigma}))
    assert c.get("CCO", 0, "implicit", {}) is None


def test_cache_checks_identity_and_cleans_failed_publication(monkeypatch, tmp_path):
    from metag.energetics import species_cache as c
    from pathlib import Path
    monkeypatch.setattr(c, "CACHE_DIR", str(tmp_path))
    monkeypatch.setattr(c, "_ENABLED", True)
    c.put("CCO", 0, "implicit", {}, -50, 1)
    path, _ = c._path("CCO", 0, "implicit", {})
    record = json.loads(Path(path).read_text())
    record["settings"] = {"wrong_model": True}
    Path(path).write_text(json.dumps(record))
    assert c.get("CCO", 0, "implicit", {}) is None
    c.put("CO", 0, "implicit", {}, -50, 1, meta={"bad": object()})
    assert not list(tmp_path.glob("*.tmp"))
