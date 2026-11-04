"""Scientific identity controls from the independent audit; no molecular energies mocked as evidence."""
import copy
import json
import random
from pathlib import Path

import pytest
from rdkit import Chem
from metag.routing.planner import route_reaction
from metag.chem import reaction_residual

CASES = json.loads((Path(__file__).parents[1] / 'fixtures/audit_representation_cases.json').read_text())


def signature(rx):
    return {
        'species': sorted((name, c, q, Chem.MolToSmiles(Chem.MolFromSmiles(s)))
                          for name, (c, q, s) in rx['species'].items()),
        'pka_sites': sorted(tuple(x) for x in rx.get('pka_sites', [])),
        'n_Hplus': rx.get('n_Hplus', 0),
    }


@pytest.mark.parametrize('case', CASES, ids=lambda x: x['reaction'])
def test_equivalent_drawings_and_dictionary_order_share_chemical_route(case):
    original = copy.deepcopy(case['original_input'])
    expected = signature(route_reaction(original, log=lambda *_: None)[0])
    for trial in range(12):
        alternate = copy.deepcopy(case['atom_order_only_input'])
        rng = random.Random(trial)
        rows = list(alternate['species'].items()); rng.shuffle(rows)
        for name, (c, q, s) in rows:
            mol = Chem.MolFromSmiles(s)
            order = list(range(mol.GetNumAtoms())); rng.shuffle(order)
            alternate['species'][name] = [c, q, Chem.MolToSmiles(
                Chem.RenumberAtoms(mol, order), canonical=False, isomericSmiles=True)]
        alternate['species'] = {name: alternate['species'][name] for name, _ in rows}
        routed = route_reaction(alternate, log=lambda *_: None)[0]
        assert signature(routed) == expected
        assert reaction_residual(routed['species'], routed.get('n_Hplus', 0)) == ({}, 0)
    assert original == case['original_input']


@pytest.mark.parametrize('case', CASES, ids=lambda x: x['reaction'])
@pytest.mark.parametrize('extent', [-1, 0.5, 2])
def test_routing_preserves_direction_and_extent_after_canonicalization(case, extent):
    rx = copy.deepcopy(case['original_input'])
    baseline = route_reaction(rx, log=lambda *_: None)[0]
    rx['species'] = {n: [extent*c, q, s] for n, (c, q, s) in rx['species'].items()}
    rx['n_Hplus'] = extent*rx.get('n_Hplus', 0)
    if rx.get('pka_sites'):
        rx['pka_sites'] = [[('prod' if x[0]=='react' else 'react') if extent<0 else x[0],
                            x[1], x[2] if len(x)>2 else 'acid',
                            abs(extent)*(x[3] if len(x)>3 else 1.)] for x in rx['pka_sites']]
    routed = route_reaction(rx, log=lambda *_: None)[0]
    routed['species'] = {n: [c/extent,q,s] for n,(c,q,s) in routed['species'].items()}
    routed['n_Hplus'] /= extent
    if routed.get('pka_sites'):
        routed['pka_sites'] = [[('prod' if x[0]=='react' else 'react') if extent<0 else x[0],
                               x[1], x[2] if len(x)>2 else 'acid',
                               (x[3] if len(x)>3 else 1.)/abs(extent)] for x in routed['pka_sites']]
    # Normalize implicit unit multiplicities in the baseline too.
    for obj in (baseline, routed):
        obj['pka_sites'] = [list(x[:2])+[x[2] if len(x)>2 else 'acid',x[3] if len(x)>3 else 1.]
                            for x in obj.get('pka_sites',[])]
    assert signature(routed) == signature(baseline)


@pytest.mark.parametrize('case', CASES, ids=lambda x: x['reaction'])
def test_species_aliases_do_not_change_chemical_route(case):
    rx = copy.deepcopy(case['original_input'])
    renamed = dict(rx, species={f'alias_{i}':v for i,v in enumerate(reversed(list(rx['species'].values())))})
    def chemistry(obj):
        sig = signature(obj)
        return sorted(tuple(x[1:]) for x in sig['species']), sig['pka_sites'], sig['n_Hplus']
    assert chemistry(route_reaction(rx, log=lambda *_:None)[0]) == chemistry(
        route_reaction(renamed, log=lambda *_:None)[0])
