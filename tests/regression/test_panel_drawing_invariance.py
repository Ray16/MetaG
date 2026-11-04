"""Every packaged panel reaction redrawn fully neutral (PubChem-style) must route to the same chemistry.

Known exceptions are listed explicitly so the set can only shrink. They are acid-drawing effects on
spectator truncation and partially protonated phosphonates (rxn00182/184, 02305, 02403, 02930, 00462,
02685, 02900), a deprotonated aromatic N (02902), and phosphagens that are refused under both drawings
and differ only in the refusal message (00396, 01371, 04785)."""
import json
from importlib import resources

import pytest
from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize

from metag.chem import reaction_residual
from metag.routing.planner import route_reaction

KNOWN = {"rxn00182", "rxn00184", "rxn00396", "rxn00462", "rxn01371", "rxn02305", "rxn02403",
         "rxn02685", "rxn02900", "rxn02902", "rxn02930", "rxn04785"}


def _snap(reaction):
    try:
        rx, routes, _ = route_reaction(reaction, log=lambda *a: None)
    except Exception as exc:                       # same exception class under both drawings is fine
        return "EXC " + type(exc).__name__
    return json.dumps({"sp": sorted(tuple(v) for v in rx["species"].values()),
                       "nh": round(float(rx.get("n_Hplus") or 0), 9),
                       "pk": sorted(map(str, rx.get("pka_sites") or [])),
                       "err": routes["errors"], "ph0": routes.get("ph0")}, sort_keys=True)


@pytest.mark.integration
def test_panel_routing_invariant_to_neutral_redrawing():
    RDLogger.DisableLog("rdApp.*")
    panel = json.loads(resources.files("metag.data").joinpath("reactions_opentecr_std.json").read_text())
    uncharger = rdMolStandardize.Uncharger()
    differ = set()
    for rid, reaction in panel.items():
        species = {}
        for name, (c, _, smi) in reaction["species"].items():
            mol = uncharger.uncharge(Chem.MolFromSmiles(smi))
            species[name] = [c, Chem.GetFormalCharge(mol), Chem.MolToSmiles(mol)]
        _, dq = reaction_residual(species, 0)
        redrawn = dict(reaction, species=species, n_Hplus=-dq)
        if reaction_residual(species, redrawn["n_Hplus"])[0]:
            continue
        if _snap(reaction) != _snap(redrawn):
            differ.add(rid)
    assert differ <= KNOWN, sorted(differ - KNOWN)
