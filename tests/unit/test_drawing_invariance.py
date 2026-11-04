"""Routing must not depend on whether a basic amine/guanidine is drawn neutral or protonated."""
import pytest
from metag.chem import reaction_residual
from metag.routing.pka_transform import canonicalize_basic_amines
from metag.routing.planner import route_reaction

ATP = "Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)([O-])OP(=O)([O-])OP(=O)([O-])[O-])[C@@H](O)[C@H]1O"
ADP = "Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)([O-])OP(=O)([O-])[O-])[C@@H](O)[C@H]1O"


@pytest.mark.parametrize("smiles, expected", [
    ("N", "[NH4+]"),
    ("CN(C)C", "C[NH+](C)C"),
    ("NCCO", "[NH3+]CCO"),
    ("N=C(N)NCCC[C@H](N)C(=O)O", "NC(=[NH2+])NCCC[C@H]([NH3+])C(=O)O"),   # arginine
    ("CN(CC(=O)O)C(N)=N", "CN(CC(=O)O)C(N)=[NH2+]"),                      # creatine
    ("CC(N)=N", "CC(N)=[NH2+]"),                                            # acetamidine
])
def test_strong_bases_protonated(smiles, expected):
    assert canonicalize_basic_amines(smiles) == (expected, expected.count("+"))


@pytest.mark.parametrize("smiles", [
    "Nc1ccccc1",                          # aniline, pKa 4.6
    "CC(N)=O", "NC(N)=O",                 # amide, urea
    "Nc1ncnc2[nH]cnc12",                  # adenine
    "NO",                                 # hydroxylamine, pKa 6
    "CN1CC(=O)N=C1N",                     # creatinine (cyclic amidine), pKa 4.8
    "CC(=O)NC(N)=N",                      # acylguanidine
    "CN(CC(=O)O)C(=N)NP(=O)(O)O",         # phosphocreatine (phosphagen)
    "NC(=O)C1=CN(C)C=CC1",                # dihydronicotinamide ring N
    "[NH3+]CC(=O)[O-]",                   # already protonated
])
def test_weak_or_nonbases_untouched(smiles):
    assert canonicalize_basic_amines(smiles) == (smiles, 0)


def _route(species):
    _, dq = reaction_residual(species, 0)
    rx, routes, _ = route_reaction({"species": species, "n_Hplus": -dq}, log=lambda *a: None)
    return (sorted(tuple(v) for v in rx["species"].values()), rx["n_Hplus"],
            sorted(map(str, rx.get("pka_sites") or [])), routes["errors"])


@pytest.mark.parametrize("neutral, cation, q", [("N", "[NH4+]", 1)])
def test_amide_hydrolysis_drawing_invariant(neutral, cation, q):
    base = {"am": [-1, 0, "CC(N)=O"], "w": [-1, 0, "O"], "ac": [1, -1, "CC(=O)[O-]"]}
    a = _route(dict(base, nh=[1, 0, neutral]))
    b = _route(dict(base, nh=[1, q, cation]))
    assert a == b
    assert a[1] == -1          # the ammonium proton is taken from the bath


def test_ethanolamine_kinase_drawing_invariant():
    def rxn(etn, q):
        return {"etn": [-1, q, etn], "ATP": [-1, -4, ATP],
                "pe": [1, -1, "[NH3+]CCOP(=O)([O-])[O-]"], "ADP": [1, -3, ADP]}
    assert _route(rxn("NCCO", 0)) == _route(rxn("[NH3+]CCO", 1))


def test_canonicalization_recorded_in_routes():
    _, dq = reaction_residual({"am": [-1, 0, "CC(N)=O"], "w": [-1, 0, "O"],
                               "ac": [1, -1, "CC(=O)[O-]"], "nh": [1, 0, "N"]}, 0)
    _, routes, _ = route_reaction({"species": {"am": [-1, 0, "CC(N)=O"], "w": [-1, 0, "O"],
                                               "ac": [1, -1, "CC(=O)[O-]"], "nh": [1, 0, "N"]},
                                   "n_Hplus": -dq}, log=lambda *a: None)
    assert routes["amine_canonicalized"] == ["nh"]
