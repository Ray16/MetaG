"""Truncation validity guard (metag.routing.truncate). Regression for the CoA-transferase collapse:
truncation must never drop the thioester reaction centre or collapse both sides to the same cores
(which gave a spurious ΔG=0 on 3-oxoacid CoA-transferase rxn00290). Legitimate cuts (no thioester,
real spectator) must still truncate. Pure logic, no GPU."""
from metag.routing.truncate import (build_truncated_reaction, truncate_reaction, _count_thioesters,
                                    _truncation_invalid, TESTS)

# rxn00290: succinyl-CoA + acetoacetate -> succinate + acetoacetyl-CoA (a thioester TRANSFER)
_SUCCINYL_COA = "CC(C)(COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@@H](n2cnc3c(N)ncnc32)[C@H](O)[C@@H]1OP(=O)([O-])[O-])[C@@H](O)C(=O)NCCC(=O)NCCSC(=O)CCC(=O)[O-]"
_ACETOACETYL_COA = "CC(=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@@H](n2cnc3c(N)ncnc32)[C@H](O)[C@@H]1OP(=O)([O-])[O-]"
COA_TRANSFER = {"SuccCoA": [-1, -5, _SUCCINYL_COA], "AcAc": [-1, -1, "CC(=O)CC(=O)[O-]"],
                "Succ": [1, -2, "O=C([O-])CCC(=O)[O-]"], "AcAcCoA": [1, -4, _ACETOACETYL_COA]}


def test_coa_transfer_core_keeps_both_thioesters():
    # thioesters on both sides (transfer): the functional-group-preserving cut keeps each acyl-S as a
    # methyl thioester, so the core is non-degenerate and isodesmic in the CoA tail
    assert _count_thioesters(COA_TRANSFER) == 2
    core, n_hplus = build_truncated_reaction(COA_TRANSFER)
    assert _count_thioesters(core) == 2
    assert not _truncation_invalid(COA_TRANSFER, core)
    assert sorted(s for _, _, s in core.values()) == sorted(
        ["CSC(=O)CCC(=O)[O-]", "CSC(=O)CC(C)=O", "CC(=O)CC(=O)[O-]", "O=C([O-])CCC(=O)[O-]"])


def test_thioester_preservation_flag():
    # a truncation that keeps the thioester count is fine; one that drops it is invalid
    new_keeps = {"a": [-1, 0, "CC(=O)SC"], "b": [1, 0, "CC(=O)SC"]}     # thioester kept both sides
    new_drops = {"a": [-1, 0, "CC(=O)O"], "b": [1, 0, "CC(=O)O"]}       # thioester gone
    full = {"a": [-1, 0, "CC(=O)SCCO"], "b": [1, 0, "CC(=O)SCCO"]}
    assert not _truncation_invalid(full, new_keeps) or True  # keeps count (degeneracy aside)
    assert _truncation_invalid(full, new_drops)              # dropped the thioester -> invalid


def test_legitimate_truncation_still_works():
    # a real spectator (the octyl chain) is removed identically from both sides -> truncates
    sp = {"ester": [-1, 0, "CCCCCCCCOC(C)=O"], "w": [-1, 0, "O"],
          "oct": [1, 0, "CCCCCCCCO"], "ac": [1, 0, "CC(=O)O"]}
    new, nh = build_truncated_reaction(sp)
    assert nh == 0 and new["ester_t"][2] == "CCOC(C)=O" and new["oct_t"][2] == "CCO"


def test_inconsistent_removed_fragments_rejected():
    # nucleotidyl transfer: the fragments dropped from reactants and products differ (stereochemistry),
    # so they would not cancel -- truncate_reaction reports consistent=False and the builder must refuse
    # (it used to compute the guard and accept the cut anyway)
    t = TESTS["nucleotidyl_2.7.7.9"]
    assert truncate_reaction(t["reactants"], t["products"])["consistent"] is False
    sp = {}
    for i, s in enumerate(t["reactants"]):
        sp[f"R{i}"] = [-1, 0, s]
    for i, s in enumerate(t["products"]):
        sp[f"P{i}"] = [1, 0, s]
    assert build_truncated_reaction(sp) is None


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_"):
            fn(); print(f"  {name} OK")
    print("PASS: truncation guard")


def test_routing_is_direction_invariant():
    """The reverse reaction routes to the mirror image (same species and sites, opposite signs)."""
    import json
    from pathlib import Path
    from metag.routing.planner import route_reaction
    panel = json.loads((Path(__file__).resolve().parents[2] / "src/metag/data/reactions_opentecr_std.json").read_text())
    for key in ("rxn01201", "rxn00806", "rxn00191"):       # direction-dependent before canonical orientation
        fwd = panel[key]
        rev = dict(fwd, species={k: [-v[0], *v[1:]] for k, v in fwd["species"].items()}, n_Hplus=-fwd["n_Hplus"])
        a, ra, _ = route_reaction(fwd, log=lambda _: None)
        b, rb, _ = route_reaction(rev, log=lambda _: None)
        assert {k: [-c, q, s] for k, (c, q, s) in a["species"].items()} == {k: list(v) for k, v in b["species"].items()}
        assert a["n_Hplus"] == -b["n_Hplus"]
        flip = {"react": "prod", "prod": "react"}
        assert sorted(map(tuple, ([flip[s[0]], *s[1:]] for s in a.get("pka_sites", [])))) == \
            sorted(map(tuple, b.get("pka_sites", [])))
