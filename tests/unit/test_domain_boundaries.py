"""Fail-closed boundaries and reusable preflight, independently of reaction labels."""
import pytest
from rdkit import Chem

from metag.routing.planner import route_reaction
from metag.routing import pka_transform as pk
from metag.routing.pka_domain import UnsupportedSiteError
from metag.tools.preflight import audit_reactions, main
from metag.validation.protonation import inspect_species


def reaction():
    return {"species": {"a": [-1, 0, "CCO"], "b": [1, 0, "COC"]}}


@pytest.mark.parametrize("smiles", ["CP(=O)(O)O", "CP(=O)([O-])[O-]", "OP(O)O"])
def test_unsupported_p_not_misclassified_as_phosphate(smiles):
    with pytest.raises(UnsupportedSiteError, match="phosphoryl"):
        pk._neutralize_v2(smiles)


def test_ladder_overflow_does_not_invent_or_drop_values():
    with pytest.raises(UnsupportedSiteError, match="exceed"):
        pk._assign_ladder([0, 1, 2], [1.5, 6.5], "test")


@pytest.mark.parametrize("smiles, expected", [
    ("OP(=O)(O)O", [2.15, 7.2, 12.35]),
    ("COP(=O)(O)O", [1.5, 6.5]),
    ("COP(=O)(O)OC", [1.5]),
])
def test_supported_phosphate_values_preserved(smiles, expected):
    assert pk._neutralize_v2(smiles)[1] == expected


def test_preflight_reports_neutral_reference_does_not_prove_coverage():
    row = inspect_species("n1ccccc1")
    assert row["reference_charge"] == 0
    assert "aromatic_n_speciation_incomplete" in row["issues"]
    assert row["status"] == "review"
    assert inspect_species("CP(=O)(O)O")["status"] == "unsupported"


def test_preflight_checks_input_and_routed_species_without_labels():
    report = audit_reactions({"arbitrary": reaction(), "bad": {}})
    assert report["reactions"]["bad"]["status"] == "invalid"
    assert report["n_unique_structures"] == 2
    assert "routed_species" in report["reactions"]["arbitrary"]
    assert report["source_sha256"]
    renamed = audit_reactions({"another": dict(reaction(), exp=[99999], note="anything")})
    assert report["reactions"]["arbitrary"] == renamed["reactions"]["another"]


def test_missing_preflight_input_is_not_an_empty_success(tmp_path):
    with pytest.raises(FileNotFoundError):
        main([str(tmp_path / "missing.json"), "-o", str(tmp_path / "report.json")])


@pytest.mark.parametrize("phosphoryl", ["C", "N"])
def test_unsupported_route_returns_no_estimate_before_species_energy(monkeypatch, phosphoryl):
    from metag import pipeline
    def forbidden(*args, **kwargs):
        raise AssertionError("species energy must not run for unsupported pKa chemistry")
    monkeypatch.setattr(pipeline, "_species_free_energies", forbidden)
    rx = {"species": {"ester": [-1, 0, f"COP({phosphoryl})(=O)O"], "water": [-1, 0, "O"],
                      "acid": [1, 0, f"{phosphoryl}P(=O)(O)O"], "methanol": [1, 0, "CO"]}}
    backend = type("Backend", (), {"model_name": "fake", "producer": {"cacheable": False}})()
    result = pipeline.score_reaction(rx, backend, log=lambda _: None, allow_truncate=False)
    assert result["dG"] is None
    assert any(("unsupported phosphoryl" if phosphoryl == "C" else "quarantined") in error
               for error in result["routes"]["errors"])
