"""base_reference='global' gives every species one protonation reference state in every reaction."""
import json
from collections import Counter
from fractions import Fraction
from importlib import resources

import pytest

from metag.routing.planner import route_reaction

# Glutamate dehydrogenase, malate dehydrogenase, aspartate transaminase, aspartate ammonia-lyase and
# fumarate hydratase form a closed input cycle. Under the reaction-class base gate, glutamate and
# aspartate are routed neutral in some of these reactions and as ammonium in others.
WITNESS = {"rxn00182": 1, "rxn00248": -1, "rxn00260": 1, "rxn00347": -1, "rxn00799": 1}


def _panel():
    return json.loads(resources.files("metag.data").joinpath("reactions_opentecr_std.json").read_text())


def _net(cycle, base_reference, allow_truncate=True):
    panel, species, sites, n_h = _panel(), Counter(), Counter(), Fraction(0)
    for rid, weight in cycle.items():
        rx, routes, _ = route_reaction(panel[rid], allow_truncate=allow_truncate, log=lambda *a: None,
                                       base_reference=base_reference)
        assert not routes["errors"] and routes["base_reference"] == base_reference
        w = Fraction(weight)
        for c, q, s in rx["species"].values():
            species[(s, q)] += Fraction(c).limit_denominator() * w
        n_h += Fraction(rx["n_Hplus"]).limit_denominator() * w
        for site in rx.get("pka_sites") or []:
            sign = 1 if site[0] == "prod" else -1
            sites[(site[1], site[2])] += sign * Fraction(site[3] if len(site) > 3 else 1).limit_denominator() * w
    return ({k: v for k, v in species.items() if v}, {k: v for k, v in sites.items() if v}, n_h)


@pytest.mark.parametrize("allow_truncate", [True, False])
def test_global_reference_closes_amino_acid_witness_cycle(allow_truncate):
    assert _net(WITNESS, "global", allow_truncate) == ({}, {}, 0)


def test_gated_reference_leaves_mixed_amine_states():
    species, _, _ = _net(WITNESS, "gated")
    charges = sorted(q for _, q in species)
    assert charges == [0, 0, 1, 1]          # neutral and ammonium forms of Glu and Asp do not cancel


def test_invalid_base_reference_rejected():
    with pytest.raises(ValueError, match="base_reference"):
        route_reaction(_panel()["rxn00182"], log=lambda *a: None, base_reference="sometimes")
