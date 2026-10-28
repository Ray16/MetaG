"""Regression tests for the 2026-09-30 code-review fixes (pure logic / mocked QM, no GPU).

Each test pins one invariant that was silently violated before: routed reactions stay balanced, invalid
reactions return no estimate, water is water however it is written, fractional coefficients are not truncated, the pKa ladders of phosphoryl-N /
acyl-phosphate / carbonic acid, the activated-aldehyde hydration, and the
float32-aware conformer dedup."""
import numpy as np
import pytest
from rdkit import Chem

from metag.chem import is_water, is_balanced, reaction_residual
from metag.routing import pka_transform as PK
from metag.routing.planner import route_reaction
from metag.settings import Settings

P = pytest.importorskip("metag.pipeline")


# ---------------------------------------------------------------- shared helpers (metag.chem)
def test_water_any_spelling():
    for s in ("O", "[OH2]", "[H]O[H]"):
        assert is_water(s, 0)
    assert not is_water("O", -1) and not is_water("[OH-]", -1) and not is_water("OO", 0)


def test_balance_with_protons_and_fractions():
    # 2 H2O2 -> 2 H2O + O2 written with fractional coefficients
    assert is_balanced({"a": [-1, 0, "OO"], "b": [1, 0, "O"], "c": [0.5, 0, "O=O"]})
    assert not is_balanced({"a": [-1, 0, "OO"], "b": [1, 0, "O"]})
    # lactate + NAD-ring model -> pyruvate + dihydro ring + H+  (proton closes H and charge)
    lac = {"lac": [-1, -1, "CC(O)C(=O)[O-]"], "nad": [-1, 1, "C[n+]1cccc(C(N)=O)c1"],
           "pyr": [1, -1, "CC(=O)C(=O)[O-]"], "nadh": [1, 0, "CN1C=CCC(C(N)=O)=C1"]}
    assert is_balanced(lac, n_hplus=1) and not is_balanced(lac, n_hplus=0)
    res, _ = reaction_residual({"x": [-1, -1, "CCO"]})          # declared charge != SMILES charge
    assert "declared_charge_mismatch" in res


# ---------------------------------------------------------------- input validation + fail-closed
def test_validate_reaction_rejects_charge_mismatch():
    with pytest.raises(ValueError):
        P.validate_reaction({"species": {"a": [-1, -1, "CCO"], "b": [1, 0, "CC=O"]}})
    with pytest.raises(ValueError):
        P.validate_reaction({"species": {"a": [-1, 0, "CC(=O)[O-]"]}})
    with pytest.raises(ValueError):
        P.validate_reaction({"species": {"a": [-1, 0, "not a smiles"]}})
    P.validate_reaction({"species": {"a": [-1, -1, "CC(=O)[O-]"]}})   # consistent -> no error


def test_validate_reaction_rejects_disconnected_species():
    with pytest.raises(ValueError, match="disconnected SMILES"):
        P.validate_reaction({"species": {"a": [-1, 0, "CCO.CC"], "b": [1, 0, "CCCCO"]}})


@pytest.mark.parametrize("record", [
    [-1, 0],
    [0, 0, "CCO"],
    [float("nan"), 0, "CCO"],
    [-1, 0.5, "CCO"],
    [-1, 0, ""],
])
def test_validate_reaction_rejects_malformed_species_records(record):
    with pytest.raises(ValueError):
        P.validate_reaction({"species": {"a": record}})
    with pytest.raises(ValueError):
        P.validate_reaction({"n_Hplus": float("inf"), "species": {"a": [-1, 0, "CCO"]}})


G = {"CCO": -100.0, "CC=O": -50.0, "[HH]": -5.0, "CC(=O)O": -300.0, "CCOC(C)=O": -350.0}


class Backend:
    """Deterministic species energies keyed by canonical SMILES (no QM); unlisted states (hydrates) get
    a prohibitive energy, i.e. zero population."""
    model_name = "fake"
    producer = {"cacheable": False}

    def __init__(self, settings=None):
        self.settings = settings

    def species(self, q, smi, label, warnings, diagnostics):
        return G.get(Chem.MolToSmiles(Chem.MolFromSmiles(smi)), 1e4), 0.5

    def liquid_water(self):
        return -200.0


def _score(sp, n_h=0, backend=None):
    return P.score_reaction({"note": "test", "n_Hplus": n_h, "species": sp}, backend or Backend(),
                            log=lambda *_: None, allow_truncate=False)


def test_unbalanced_reaction_returns_no_estimate():
    # Proton count wrong -> reject before QM, with no number.
    class NoQM(Backend):
        def species(self, *args, **kwargs):
            raise AssertionError("unbalanced input reached species QM")
        liquid_water = species

    r = _score({"A": [-1, 0, "CCO"], "B": [1, 0, "CC=O"], "H2": [1, 0, "[HH]"]}, n_h=1, backend=NoQM())
    assert r["suspect"] and r["dG"] is None and r["dG_raw"] is None
    assert r["dG_raw_unreliable"] is None
    assert r["routes"]["input_balanced"] is False
    assert "residual" in r["suspect"]


def test_water_spelling_does_not_change_dG():
    base = {"A": [-1, 0, "CC(=O)O"], "B": [-1, 0, "CCO"], "E": [1, 0, "CCOC(C)=O"]}
    a = _score(dict(base, W=[1, 0, "O"]))["dG_raw"]
    b = _score(dict(base, W=[1, 0, "[H]O[H]"]))["dG_raw"]
    assert a is not None and a == b                          # "[H]O[H]" used to be scored as a 1 M solute


def test_balance_guard_reverts_unbalanced_rewrite(monkeypatch):
    # a cofactor core that drops atoms (the mixed-disulfide failure) must be reverted, not scored
    import metag.routing.cofactor_cores as CC
    rx = {"note": "t", "n_Hplus": 0, "species": {"A": [-1, 0, "CCO"], "B": [1, 0, "CC=O"], "H2": [1, 0, "[HH]"]}}
    monkeypatch.setattr(CC, "cofactor_ring", lambda sp: {"A": [-1, 0, "CCCO"], "B": [1, 0, "CC=O"], "H2": [1, 0, "[HH]"]})
    new, routes, _ = route_reaction(rx, allow_truncate=False, log=lambda *_: None)
    assert new["species"] == rx["species"] and routes["cofactor_ring"] is False
    assert any("cofactor_ring" in e and "reverted" in e for e in routes["reverted"]) and not routes["errors"]


# ---------------------------------------------------------------- configuration fingerprint
def test_config_fingerprints_numeric_acceptance_controls():
    base = P.effective_config(Backend(Settings()))
    changed = P.effective_config(Backend(Settings.from_env({"DG_SANITY_KJ": "321"})))
    assert changed["settings"]["dg_sanity_kj"] == 321 != base["settings"]["dg_sanity_kj"]


# ---------------------------------------------------------------- pKa ladders
def _pkas(smi):
    return sorted(p for _, p in PK._classify_species(PK._canonicalize_maxanion(smi))[1])


def test_phosphoramidate_acylphosphate_carbonate_ladders():
    with pytest.raises(ValueError, match="phosphoramidate.*quarantined"):
        _pkas("CN(CC(=O)[O-])C(N)=[NH+]P(=O)([O-])[O-]")
    assert _pkas("CC(=O)OP(=O)([O-])[O-]") == PK.ACYL_P_LADDER               # acetyl phosphate
    assert _pkas("O=C([O-])O") == PK.CARBONATE_LADDER                         # bicarbonate
    assert _pkas("COP(=O)([O-])[O-]") == PK.P_LADDER[1]                       # plain monoester unchanged
    assert _pkas("O=P([O-])([O-])O") == PK.P_LADDER[0][:3]                     # free Pi unchanged
    assert _pkas("CC(=O)[O-]") == [PK.CARBOXYL_PKA]                          # carboxylate unchanged


def test_ph0_fractional_coefficients_carry_multiplicity():
    # 1/2 acetate + 1/4 O2 -> 1/2 glycolate (balanced, fractional)
    sp = {"a": [-0.5, -1, "CC(=O)[O-]"], "o2": [-0.25, 0, "O=O"], "b": [0.5, -1, "OCC(=O)[O-]"]}
    assert is_balanced(sp)
    ns, sites, nh = PK.build_ph0_reaction(sp, 0)
    assert all(len(s) == 4 and s[3] == 0.5 for s in sites) and len(sites) == 2


def test_redox_proton_is_carried_not_refused():
    # lactate + NAD+ -> pyruvate + NADH + H+ (ring model): neutralized reaction is short one H;
    # carried as n_H+ instead of refusing pH-0 (the same metabolites are pH-0 routed elsewhere)
    lac = {"lac": [-1, -1, "CC(O)C(=O)[O-]"], "nad": [-1, 1, "C[n+]1cccc(C(N)=O)c1"],
           "pyr": [1, -1, "CC(=O)C(=O)[O-]"], "nadh": [1, 0, "CN1C=CCC(C(N)=O)=C1"]}
    why = []
    out = PK.build_ph0_reaction(lac, 1, base=False, why=why)
    assert out is not None, why
    ns, sites, nh = out
    assert is_balanced(ns, nh)


# ---------------------------------------------------------------- aldehyde hydration
def test_gem_diol_hydrates_the_activated_aldehyde():
    from metag.routing.aldehyde_hydration import gem_diol
    assert gem_diol("O=CCCC(=O)C=O") == "O=CCCC(=O)C(O)O"       # alpha-keto CHO, not the alkyl CHO
    assert gem_diol("O=CC(=O)[O-]") == "O=C([O-])C(O)O"          # glyoxylate
    assert gem_diol("O=CCCC") is None                            # non-activated: gated out


# ---------------------------------------------------------------- conformer dedup + connectivity
def test_dedup_energy_tolerance_exceeds_float32_step():
    from metag.energetics.conformers import UniqueMinima, EV2KJ
    u = UniqueMinima()
    e_kj = -72000.0 * EV2KJ                                        # ATP-size total energy
    step = float(np.spacing(np.float32(72000.0))) * EV2KJ
    assert step > 0.5                                              # the old 0.5 kJ tolerance was below it
    assert u._e_tol(e_kj) >= 2 * step and u._e_tol(-10.0) == u.e_tol
    assert UniqueMinima(e_tol=-1.0)._e_tol(e_kj) < 0               # legacy never-merge preserved


def test_connectivity_check_flags_proton_transfer():
    from rdkit import Chem
    from metag.energetics.conformers import pool_confs, bond_graph, same_connectivity
    smi = "[NH3+]CC(=O)[O-]"
    m = Chem.AddHs(Chem.MolFromSmiles(smi)); ref = bond_graph(m)
    a = pool_confs(smi, 0, 1, 1)[0]
    assert same_connectivity(a, ref)
    on = next(x.GetIdx() for x in m.GetAtoms() if x.GetSymbol() == "O" and x.GetFormalCharge() == -1)
    hn = next(n.GetIdx() for n in m.GetAtomWithIdx(0).GetNeighbors() if n.GetSymbol() == "H")
    x = a.get_positions(); x[hn] = x[on] + np.array([0.97, 0.0, 0.0]); a.set_positions(x)
    assert not same_connectivity(a, ref)


def test_rescore_source_fingerprint_changes_with_scoring_source(tmp_path, monkeypatch):
    from metag.validation.provenance import source_tree_fingerprint
    package = tmp_path / "metag"
    package.mkdir(parents=True)
    scoring = package / "pipeline.py"
    scoring.write_text("VERSION = 1\n")
    runner = tmp_path / "runner.py"
    runner.write_text("runner = 1\n")
    before = source_tree_fingerprint(package, [runner])
    scoring.write_text("VERSION = 2\n")
    assert source_tree_fingerprint(package, [runner]) != before
