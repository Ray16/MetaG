"""pKa and hydration rules (no QM): weak-acid anion neutralisation, alpha,beta-unsaturated
carboxyl pKa class, aryl-amine basicity in the base-path gate, general carbonyl hydration sites, the
multi-site hydration mixture, and the extra-force hook of the batched optimiser."""
import math
import os

import numpy as np

from metag.routing import aldehyde_hydration as ah
from metag.routing import pka_transform as pk


def test_thiolate_and_phenolate_are_neutralised():
    gsh = "[NH3+][C@@H](CCC(=O)N[C@@H](C[S-])C(=O)NCC(=O)[O-])C(=O)[O-]"
    smi, pkas, q = pk._neutralize(gsh)
    assert "[S-]" not in smi and "CS" in smi
    assert pk.THIOL_PKA in pkas
    tyr = "N[C@@H](Cc1ccc([O-])cc1)C(=O)[O-]"
    smi, pkas, q = pk._neutralize(tyr)
    assert q == 0 and pk.PHENOL_PKA in pkas
    assert pk._neutralize("CC(=O)[S-]")[2] == -1          # thiocarboxylate is not a thiolate


def test_unsaturated_carboxyl_pka():
    # class rule: mesaconate, acrylate; malate rule unchanged
    assert pk._neutralize("C/C(=C\\C(=O)[O-])C(=O)[O-]")[1] == [3.75, 3.75]     # default: per-site class rule
    assert pk._neutralize("C=CC(=O)[O-]")[1] == [4.35]
    assert sorted(pk._neutralize("O=C([O-])C[C@H](O)C(=O)[O-]")[1]) == [3.8, 4.75]


def test_alpha_keto_carboxyl_uses_keto_microstate_pka():
    assert pk._neutralize("CC(=O)C(=O)[O-]")[1] == [1.8]                       # Lopalco et al. 2016


def test_arylamine_is_not_a_basic_amine():
    adenylosuccinate_synthase = {
        "Asp": (-1, -1, "[NH3+][C@@H](CC(=O)[O-])C(=O)[O-]"),
        "IMP": (-1, -2, "O=c1[nH]cnc2c1ncn2[C@@H]1O[C@H](COP(=O)([O-])[O-])[C@@H](O)[C@H]1O"),
        "S-AMP": (1, -4, "O=C([O-])C[C@H](Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)([O-])[O-])[C@@H](O)"
                         "[C@H]1O)C(=O)[O-]"),
        "H2O": (1, 0, "O"),
    }
    assert pk._amine_cn_change(adenylosuccinate_synthase) == -1         # the Asp ammonium is destroyed


def test_hydration_sites():
    assert ah.hydration_sites("CC(=O)C(=O)O") == []                     # alpha-keto acid: excluded
    assert ah.hydration_sites("O=C(O)CC(=O)C(=O)O") == []               # oxaloacetate: excluded
    assert [d for _, d in ah.hydration_sites("O=CC(=O)O")] == ["O=C(O)C(O)O"]           # glyoxylate kept
    assert [d for _, d in ah.hydration_sites("O=C(CO)COP(=O)(O)O")] == ["O=P(O)(O)OCC(O)(O)CO"]   # DHAP
    assert len(ah.hydration_sites("CC(=O)C=O")) == 2                                    # methylglyoxal
    assert ah.hydration_sites("CC(=O)SC") == [] and ah.hydration_sites("CC(N)=O") == []  # thioester, amide
    assert ah.hydration_sites("CC(=O)O") == []                                           # acid


def test_hydration_states_exact_enumeration():
    # every subset of sites is an explicit, computed state; no uncomputed (implicitly multiplied) states
    states = ah.hydration_states("CC(=O)C=O")                               # methylglyoxal: 2 sites
    assert sorted(n for n, _ in states) == [1, 1, 2]
    assert "CC(O)(O)C(O)O" in {s for _, s in states}
    for dg in (-12.0, -1.0, 0.0, 3.0, 20.0):                               # one state = two-state formula
        assert abs(ah.mixture_G_states(0.0, [(1, dg)], 0.0) - ah.mixture_G(0.0, dg, 0.0)) < 1e-9


def test_hydration_sites_independent_of_atom_order():
    a = sorted(s for _, s in ah.hydration_states("O=CC(C)CC(=O)CC=O"))       # three carbonyls
    b = sorted(s for _, s in ah.hydration_states("O=CCC(=O)CC(C)C=O"))       # same molecule, other order
    assert a == b and len(a) == 7                                          # 2^3 - 1 hydrated states


def test_hydration_calibration_reproduces_khyd_fit():
    # log K_exp = 0.639 log K_calc - 0.411  <=>  ΔG_exp = 0.639 ΔG_calc + 2.35 kJ (per event)
    rt_ln10 = 2.303 * 8.314e-3 * 298.15
    assert abs(ah.calibrated_dg_hyd(0.0) - 2.35) < 1e-9
    assert abs(ah.calibrated_dg_hyd(-10.0) - (-6.39 + 2.35)) < 1e-9


def test_constant_edit_changes_fingerprint(monkeypatch):
    import metag.pipeline as P
    from metag.routing import pka_transform as pk
    from types import SimpleNamespace
    backend = SimpleNamespace(model_name="fake-uma", producer={"cacheable": False})
    base = P.effective_config(backend)["pka_constants"]
    monkeypatch.setitem(pk.CARBOXYL_PKA_ALPHA, "oxo", 2.5)
    assert P.effective_config(backend)["pka_constants"] != base


def test_deployed_hydration_calibration_matches_reproducible_fit():
    # The deployed constants must match the archived, reproducible validation fit.
    import json
    from pathlib import Path
    path = Path(__file__).resolve().parents[1] / "regression" / "fixtures" / "khyd_validation.json"
    fit = json.load(open(path))["calibration_domain_fit"]
    a, b = ah.HYDRATION_CAL
    assert abs(a - fit["a"]) < 0.0015 and abs(b - fit["b_kJ_per_event"]) < 0.01


def _ethanol_minima():
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from ase import Atoms
    from metag.energetics.conformers import UniqueMinima
    m = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    cids = AllChem.EmbedMultipleConfs(m, 6, randomSeed=3)
    syms = [a.GetSymbol() for a in m.GetAtoms()]
    ats = [Atoms(symbols=syms, positions=m.GetConformer(c).GetPositions()) for c in cids]
    return UniqueMinima(template=Chem.Mol(m), e_tol=-1.0), ats   # e_tol < 0: never merge (distinct minima)


def test_dedup_replaces_geometry_with_lower_g_duplicate():
    import numpy as np
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from ase import Atoms
    from metag.energetics.conformers import UniqueMinima
    m = Chem.AddHs(Chem.MolFromSmiles("C"))
    AllChem.EmbedMolecule(m, randomSeed=1)
    syms = [a.GetSymbol() for a in m.GetAtoms()]
    p = m.GetConformer().GetPositions()
    u = UniqueMinima(template=Chem.Mol(m))
    u.add(Atoms(symbols=syms, positions=p), -100.0, -110.0)
    distorted = p.copy(); distorted[1] += [0.03, 0.0, 0.0]   # same basin (RMSD << 0.25 Å), lower G
    assert u.add(Atoms(symbols=syms, positions=distorted), -100.5, -111.0) is False
    assert u.G == [-111.0] and u.E == [-100.5]
    # the stored representative is the NEW structure (RDKit aligns it rigidly onto the old frame, so compare
    # interatomic distances, which are invariant to that alignment)
    dist = lambda x: np.linalg.norm(x[:, None] - x[None], axis=-1)
    stored = u.ref[0].GetConformer().GetPositions()
    assert np.allclose(dist(stored), dist(distorted), atol=1e-6)
    assert not np.allclose(dist(stored), dist(p), atol=1e-3)


def test_thermal_track_applies_validated_correction_to_accepted_minimum():
    from metag.energetics.conformers import boltz
    from metag.energetics.thermal_ensemble import ThermalTrack
    from metag.settings import Settings
    u, ats = _ethanol_minima()
    for i, a in enumerate(ats[:3]):
        u.add(a, -100.0 + i, -110.0 + i)                 # minimum 0 lowest, then 1, 2 (all in the window)
    tr = ThermalTrack(None, u.template, 0, 1, correction=lambda *a, **k: (50.0, {}), combine=boltz,
                      window=Settings().thermal_window_kj, maximum=Settings().thermal_max)
    tot_valid_on_1 = tr.final(u, {}, therm_valid=40.0, j_valid=1,
                                validated={j: 40.0 if j == 1 else 50.0 for j in range(len(u.G))})
    # the validated value belongs to minimum 1 (not the lowest): lowest keeps its own tracked correction
    _, c, order, _ = tr._totals(u, override=(1, 40.0))
    assert c[1] == 40.0 and c[0] == 50.0 and order[0] == 0
    assert tot_valid_on_1[2]["validated_index_is_lowest"] is False


def test_unique_minima_drop():
    u, ats = _ethanol_minima()
    for i, a in enumerate(ats[:3]):
        u.add(a, float(i), float(i))
    keep = u.drop({1})
    assert keep == [0, 2] and u.G == [0.0, 2.0] and len(u.ref) == 2
