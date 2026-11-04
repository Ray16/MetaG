"""Checks for the historical benchmark reconstruction policy."""
import importlib.util
import math
from collections import defaultdict
from pathlib import Path


spec = importlib.util.spec_from_file_location(
    "mg_bounded_reference", Path(__file__).resolve().parents[2] / "benchmark/rebuild_mg_bounded.py"
)
mg = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mg)


def row(**changes):
    data = {"K_prime": "10", "p_h": "7", "p_mg": "", "temperature": "298.15",
            "ionic_strength": "0.1", "_source_row": 2, "_reference_metadata": {}}
    data.update(changes)
    return data


def test_mg_bounds_distinguish_missing_recipe_and_explicit_free_activity():
    assert mg.mg_interval(row(), {}) == (4.0, 14.0)
    assert mg.mg_interval(row(), {"Cofactor": "MgCl2"}) == (2.0, 14.0)
    assert mg.mg_interval(row(p_mg="3"), {}) is None
    lo, hi = mg.mg_interval(row(), {"Cofactor": "c(MgCl2)=5 mM"})
    assert math.isclose(lo, -math.log10(0.005)) and hi == 14.0


def test_missing_magnesium_does_not_hide_invalid_ionic_strength():
    assert mg.base_exclusion(row(ionic_strength="-0.1"), {}) == "invalid_ionic_strength"
    assert mg.base_exclusion(row(ionic_strength="0.1"), {}) is None
    assert mg.base_exclusion(row(p_mg="3", ionic_strength="-0.1"), {}) == "invalid_ionic_strength"


def test_reconstruction_excludes_a_mg_sensitive_measurement():
    class Transformer:
        def reaction(self, rid):
            return rid

        def correction(self, rid, ph, ionic, pmg, temperature):
            return pmg / 10

    ms = {"r": {"a": -1, "b": 1}}
    index = defaultdict(list)
    index[mg.B.signature(ms["r"], {mg.B.PROTON, mg.B.WATER})].append(("r", 1))
    rows = [(row(), "C00001 = C00002", "fixture")]
    args = ({"r": {}}, 0.9, Transformer(), ms, {"C00001": "a", "C00002": "b"}, index, rows)
    records = []
    labels, excluded, used = mg.rebuild(*args, excluded_records=records)
    assert labels == {} and not used
    assert excluded == {"Mg_sensitive_beyond_threshold": 1}
    assert records == [{"reaction": "r", "source_row": 2, "equation": "C00001 = C00002",
                        "mapping": "fixture", "orientation": 1,
                        "mapping_candidates": [("r", 1)], "reason": "Mg_sensitive_beyond_threshold"}]
    accepted_records = []
    labels, excluded, used = mg.rebuild(args[0], 1.1, *args[2:], excluded_records=accepted_records)
    assert labels["r"]["n"] == 1
    assert labels["r"]["max_mg_spread"] == 1.0
    assert used["r"][0]["pMg_interval"] == (4.0, 14.0)
    assert accepted_records == [] and excluded == {}


def test_exclusion_audit_matches_counts_without_changing_labels():
    class Transformer:
        def reaction(self, rid):
            return rid

        def correction(self, rid, ph, ionic, pmg, temperature):
            return pmg / 10

    ms = {"r": {"a": -1, "b": 1}}
    index = {mg.B.signature(ms["r"], {mg.B.PROTON, mg.B.WATER}): [("r", 1)]}
    rows = [(row(_source_row=2), "C00001 = C00002", "fixture"),
            (row(_source_row=3, K_prime=""), "C00001 = C00002", "fixture"),
            (row(_source_row=4, p_mg="5"), "C00001 = C00002", "fixture"),
            (row(_source_row=5), "C00001 = C00003", "fixture")]
    args = ({"r": {}}, 0.9, Transformer(), ms, {"C00001": "a", "C00002": "b"}, index, rows)
    without_audit = mg.rebuild(*args)
    records = []
    with_audit = mg.rebuild(*args, excluded_records=records)
    assert with_audit == without_audit
    assert with_audit[1] == {"Mg_sensitive_beyond_threshold": 1, "missing_or_invalid_Kprime": 1}
    assert [(r["source_row"], r["reason"]) for r in records] == [
        (2, "Mg_sensitive_beyond_threshold"), (3, "missing_or_invalid_Kprime")]
    assert with_audit[0]["r"]["n"] == 1
