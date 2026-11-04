"""Complete solution-state partitions and shared raw sampling derivatives."""
import math

import pytest
from rdkit import Chem

from metag.constants import STD_STATE_KJ
from metag.routing.sugar_states import enumerate_sugar_states
from metag.settings import Settings
from metag.state_function import StateFunctionEstimator
from metag.state_pools import KT_KJ

GLC = 'OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O'
FRU = 'OC[C@H]1OC(O)(CO)[C@@H](O)[C@@H]1O'
G6P = 'O=P([O-])([O-])OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O'
F6P = 'O=P([O-])([O-])OC[C@H]1OC(O)(CO)[C@@H](O)[C@@H]1O'
NANA = 'CC(=O)N[C@H]1[C@H]([C@H](O)[C@H](O)CO)OC(O)(C(=O)[O-])C[C@@H]1O'


class Backend:
    model_name = 'controlled-pool-test'
    producer = {'cacheable': False, 'kind': 'synthetic'}

    def __init__(self):
        self.settings = Settings()
        self.values, self.calls = {}, []
        from metag.state_function import plan_state
        for parent, charge, offset in [(GLC, 0, 0.), (FRU, 0, 5.),
                                        (G6P, -2, 0.), (F6P, -2, 5.), (NANA, -1, 0.)]:
            for i, row in enumerate(enumerate_sugar_states(parent, charge).states):
                reference = plan_state(charge, row['smiles'])['reference_smiles']
                self.values[reference] = offset+i

    def species(self, q, smi, label, warnings, diagnostics):
        self.calls.append((q, smi))
        diagnostics['method'] = 'synthetic'
        return self.values.get(smi, 1e6), .1

    def liquid_water(self):
        return 0.


def rx(left=GLC, right=FRU, charge=0, extent=1.):
    return {'species': {'left': [-extent, charge, left], 'right': [extent, charge, right]},
            'species_state_pools': {n: {'kind': 'reducing_sugar_equilibrium',
                                      'source': 'Controlled pool declaration'} for n in ('left', 'right')}}


def test_complete_pool_reverses_and_scales_and_keeps_all_states():
    backend = Backend()
    estimator = StateFunctionEstimator(backend)
    values = [estimator.score_reaction(rx(extent=e))['dG'] for e in [1., -1., 2., .5]]
    assert values == [5., -5., 10., 2.5]
    result = estimator.score_reaction(rx())
    assert len(result['state_records']['left']['states']) == 5
    assert len(result['state_records']['right']['states']) == 5
    assert len(backend.calls) == len(set(backend.calls))


@pytest.mark.parametrize('parent,charge', [(GLC, 0), (G6P, -2), (NANA, -1)])
def test_open_and_ring_references_share_one_pool_potential(parent, charge):
    estimator = StateFunctionEstimator(Backend())
    source = enumerate_sugar_states(parent, charge)
    reaction = rx(source.open_smiles, source.states[-1]['smiles'], charge)
    result = estimator.score_reaction(reaction)
    assert result['dG'] == 0, result['suspect']
    assert result['U_samp'] == 0
    assert result['state_records']['left'] == result['state_records']['right']


def test_pool_partition_includes_the_standard_state_once_per_molecule():
    backend = Backend()
    estimator = StateFunctionEstimator(backend)
    result = estimator.score_reaction(rx())
    left = result['state_records']['left']
    # Hydrates have negligible weight in this controlled backend.
    energies = [backend.values[row['smiles']] for row in left['states']]
    low = min(energies)
    expected = low-KT_KJ*math.log(sum(math.exp(-(e-low)/KT_KJ) for e in energies))+STD_STATE_KJ
    assert left['potential_kj'] == pytest.approx(expected)
    assert sum(row['conditional_population'] for row in left['states']) == pytest.approx(1.)


def test_failed_required_state_refuses_pool_without_renormalizing():
    backend = Backend()
    opened = enumerate_sugar_states(GLC).open_smiles
    backend.values[opened] = None
    estimator = StateFunctionEstimator(backend)
    result = estimator.score_reaction(rx())
    assert result['dG'] is None
    assert 'incomplete declared sugar ensemble' in result['suspect']
    assert not any(row['plan'].get('kind') == 'reducing_sugar_equilibrium'
                   for row in estimator.ledger()['states'].values())


def test_undeclared_anomer_ambiguity_is_not_a_pool():
    backend = Backend()
    estimator = StateFunctionEstimator(backend)
    reaction = rx()
    reaction.pop('species_state_pools')
    result = estimator.score_reaction(reaction)
    assert result['dG'] is None and not backend.calls


def test_shared_constituent_covariance_is_combined_before_quadrature():
    estimator = StateFunctionEstimator(Backend())
    opened = enumerate_sugar_states(GLC).open_smiles
    reaction = rx(opened, opened)
    reaction['species_state_pools'].pop('left')
    result = estimator.score_reaction(reaction)
    assert result['dG'] is not None, result['suspect']
    pool = result['state_records']['right']
    probabilities = {r['smiles']: r['conditional_population'] for r in pool['states']}
    expected = .1*math.sqrt((probabilities[opened]-1)**2
                            +sum(p*p for s, p in probabilities.items() if s != opened))
    assert result['U_samp'] == pytest.approx(expected)


def test_pool_source_labels_do_not_select_numerical_values():
    estimator = StateFunctionEstimator(Backend())
    a = estimator.score_reaction(rx())
    other = rx()
    other['species_state_pools']['left']['source'] = 'Another declaration of the same state basis'
    b = estimator.score_reaction(other)
    assert a['state_records']['left'] == b['state_records']['left']
    assert a['dG'] == b['dG']
    assert a['declared_species_state_pools'] != b['declared_species_state_pools']


def test_nested_hydration_and_pool_sampling_derivatives_match_finite_differences():
    from metag.routing.aldehyde_hydration import enumerate_hydration_states
    opened = enumerate_sugar_states(GLC).open_smiles
    hydrated = enumerate_hydration_states(opened).states[0][1]
    values = Backend().values
    values[hydrated] = -8.
    declaration = {'kind': 'reducing_sugar_equilibrium', 'source': 'Controlled pool'}

    def evaluate(changed):
        backend = Backend()
        backend.values = changed
        return StateFunctionEstimator(backend).state(0, GLC, pool=declaration)

    baseline = evaluate(values)
    eps = 1e-4
    for row in baseline['sampling_inputs'].values():
        plus, minus = dict(values), dict(values)
        default = values.get(row['smiles'], 1e6)
        plus[row['smiles']], minus[row['smiles']] = default+eps, default-eps
        derivative = (evaluate(plus)['potential_kj']-evaluate(minus)['potential_kj'])/(2*eps)
        assert row['sensitivity'] == pytest.approx(derivative, abs=1e-7, rel=2e-5)
