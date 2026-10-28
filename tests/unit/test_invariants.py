"""Thermodynamic invariants of the assembled estimator + failure injection (mocked QM, no GPU).

  * reversal:      score(-R) = -score(R)
  * extensivity:   score(k R) = k score(R)                  (no fixed magnitude ceiling)
  * cancellation:  adding the same spectator to both sides changes nothing
  * failure injection at the external boundaries: an xtb failure on a low-energy conformer, or a species
    whose every conformer rearranges, never yields a normal result under the requested species' identity."""
import numpy as np
import pytest
from rdkit import Chem

P = pytest.importorskip("metag.pipeline")

G = {"CCO": -100.0, "CC=O": -50.0, "[HH]": -5.0, "CC(=O)O": -300.0, "CCOC(C)=O": -350.0, "C": -30.0,
     "O=C=O": -350.0, "CCCCCC": -600.0}


class _Backend:
    """Deterministic species energies keyed by canonical SMILES (no QM). States absent from the table
    (e.g. carbonyl hydrates) get a prohibitively high energy, i.e. zero population."""
    model_name = "fake"
    producer = {"cacheable": False}

    def species(self, q, smi, label, warnings, diagnostics):
        canon = Chem.MolToSmiles(Chem.MolFromSmiles(smi))
        return G.get(canon, 1e4), 0.5

    def liquid_water(self):
        return -200.0


@pytest.fixture
def mocked():
    return _Backend()


def _score(sp, n_h=0, backend=None):
    return P.score_reaction({"note": "t", "n_Hplus": n_h, "species": sp}, backend or _Backend(),
                            log=lambda *_: None, allow_truncate=False)


EST = {"A": [-1, 0, "CC(=O)O"], "B": [-1, 0, "CCO"], "E": [1, 0, "CCOC(C)=O"], "W": [1, 0, "O"]}
DECARB = {"A": [-1, 0, "CC(=O)O"], "B": [1, 0, "C"], "C": [1, 0, "O=C=O"]}


def _rev(sp):
    return {k: [-v[0], v[1], v[2]] for k, v in sp.items()}


def _mul(sp, k):
    return {n: [k * v[0], v[1], v[2]] for n, v in sp.items()}


@pytest.mark.parametrize("sp", [EST, DECARB])
def test_reversal_antisymmetric(mocked, sp):
    f, r = _score(sp), _score(_rev(sp))
    assert f["dG"] is not None and abs(f["dG"] + r["dG"]) < 0.11 and abs(f["dG_raw"] + r["dG_raw"]) < 0.11


@pytest.mark.parametrize("k", [0.5, 2, 5])
def test_extensive_in_stoichiometry(mocked, k):
    base = _score(DECARB)["dG_raw"]
    big = _score(_mul(DECARB, k))
    assert big["suspect"] is None and abs(big["dG_raw"] - k * base) < 0.11 * max(1, k)


def test_large_extensive_reaction_not_rejected(mocked, monkeypatch):
    # a physically valid ~ -1000 kJ reaction must be scored, not failed by a magnitude ceiling
    monkeypatch.setitem(G, "CCOC(C)=O", -1350.0)
    r = _score(EST)
    assert r["suspect"] is None and r["dG_raw"] < -900


def test_spectator_cancels(mocked):
    base = _score(DECARB)["dG_raw"]
    with_spec = _score(dict(DECARB, S1=[-1, 0, "CCCCCC"], S2=[1, 0, "CCCCCC"]))
    assert abs(with_spec["dG_raw"] - base) < 1e-6


def test_stages_add_up(mocked):
    r = _score(EST)
    st = r["stages"]
    assert abs(st["species_sum"] + st["proton"] + st["pka_transform"] - st["unanchored"]) < 0.02
    assert abs(st["unanchored"] - st["reported"]) < 0.02
    assert r["conditions"]["pH"] == 7.0 and r["conditions"]["ionic_strength_M"] == 0.0


# ---------------------------------------------------------------- failure injection (implicit_G internals)
S = pytest.importorskip("metag.energetics.species")


class _PU:
    pass


def _species(smi, warnings):
    from metag.settings import Settings
    return S.species_free_energy(_PU(), 0, smi, "x", settings=Settings(), model_name="fake",
                                 producer={"cacheable": False}, log=lambda *_: None, warnings=warnings)


def _fake_backend(monkeypatch, solv=lambda i: -10.0, rearrange=False):
    from ase import Atoms
    calls = {"n": 0}
    def pool(smi, q, seed, pool, spin=1):
        m = Chem.AddHs(Chem.MolFromSmiles(smi))
        from rdkit.Chem import AllChem
        AllChem.EmbedMolecule(m, randomSeed=seed)
        x = m.GetConformer().GetPositions()
        return [Atoms([a.GetSymbol() for a in m.GetAtoms()], x + 0.01 * i, info={"charge": q, "spin": spin})
                for i in range(3)]
    monkeypatch.setattr(S, "pool_confs", pool)
    monkeypatch.setattr(S, "batched_energies", lambda pu, c: np.arange(len(c), dtype=float))
    monkeypatch.setattr(S, "batched_fire", lambda pu, sel, **k: (sel, np.arange(len(sel)) * 0.1 - 1000.0,
                                                                np.ones(len(sel), bool)))
    def dg(sym, pos, q, mult=1):
        calls["n"] += 1
        return solv(calls["n"])
    monkeypatch.setattr(S, "dgsolv", dg)
    monkeypatch.setattr(S, "same_connectivity", lambda a, g: not rearrange)
    monkeypatch.setattr(S, "uma_gibbs_corr", lambda *a, **k: (5.0, {"n_imag": 0, "max_imag_cm": 0.0}))
    monkeypatch.setattr(S.cache, "_ENABLED", False)


def test_xtb_failure_on_low_energy_conformer_fails_species(monkeypatch):
    _fake_backend(monkeypatch, solv=lambda i: None)          # every solvation call fails (incl. retry)
    w = []
    g, s = _species("CCO", w)
    assert g is None and s is None


def test_all_conformers_rearranged_raises(monkeypatch):
    _fake_backend(monkeypatch, rearrange=True)
    with pytest.raises(P.SpeciesRearranged):
        _species("CCO", [])


def test_rearranged_species_without_zwitterion_fails_closed():
    # The neutral-reference retry is permitted only for a remaining zwitterion; otherwise one attempt,
    # then no estimate.
    seen = []

    class Rearranging(_Backend):
        def species(self, q, smi, *a):
            seen.append(smi); raise P.SpeciesRearranged("x", smi, 3)
    r = _score(DECARB, backend=Rearranging())
    assert r["dG"] is None and len(seen) == 1
    assert "species identity unresolved" in r["suspect"]
    assert r["routes"]["rearrangement_reroute"] is False


def test_routing_exception_fails_closed(mocked, monkeypatch):
    import metag.routing.cofactor_cores as CC
    def boom(sp):
        raise RuntimeError("injected")
    monkeypatch.setattr(CC, "cofactor_ring", boom)
    r = _score(DECARB)
    assert r["dG"] is None and "routing error" in r["suspect"]


# ---------------------------------------------------------------- spin, thermal projection, calibration
def test_spin_from_radicals_and_metal_refusal():
    from metag.energetics.conformers import spin_multiplicity
    assert spin_multiplicity("O=O", 0) == 3 and spin_multiplicity("[CH3]", 0) == 2
    assert spin_multiplicity("[CH2]", 0) == 3                 # triplet methylene (2 radical electrons)
    with pytest.raises(ValueError):
        spin_multiplicity("[Fe+2]", 2)


def test_imaginary_mode_is_detected_not_floored():
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from ase import Atoms
    from metag.energetics.thermal import internal_vib_energies
    m = Chem.AddHs(Chem.MolFromSmiles("CC")); AllChem.EmbedMolecule(m, randomSeed=1)
    mp = AllChem.MMFFGetMoleculeProperties(m); ff = AllChem.MMFFGetMoleculeForceField(m, mp)
    ff.MMFFAddTorsionConstraint(2, 0, 1, 5, False, 0.0, 0.0, 1e4); ff.Minimize(maxIts=5000)
    pos = m.GetConformer().GetPositions()
    def grad(x):
        c = m.GetConformer()
        for i, p in enumerate(x.reshape(-1, 3)):
            c.SetAtomPosition(i, p.tolist())
        return np.array(AllChem.MMFFGetMoleculeForceField(m, mp).CalcGrad()) * 0.0433641
    x0 = pos.ravel(); H = np.zeros((x0.size,) * 2)
    for k in range(x0.size):
        xp = x0.copy(); xp[k] += 1e-3; xm = x0.copy(); xm[k] -= 1e-3
        H[k] = (grad(xp) - grad(xm)) / 2e-3
    vib, imag = internal_vib_energies(Atoms([a.GetSymbol() for a in m.GetAtoms()], pos), 0.5 * (H + H.T),
                                      "nonlinear")
    assert len(vib) + len(imag) == 3 * m.GetNumAtoms() - 6 and len(imag) == 1 and imag[0] > 200
