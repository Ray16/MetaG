"""An explicit experimental-condition request must never become a fixed-target score."""
import copy

import pytest

from metag.constants import CONDITIONS
from metag.pipeline import score_reaction
from metag.state_function import StateFunctionEstimator
from metag.validation.calculation_conditions import calculation_condition_issues


class Backend:
    model_name = 'synthetic-condition-contract-control'
    producer = {'cacheable': False, 'kind': 'synthetic'}

    def __init__(self):
        self.calls = 0

    def species(self, *args):
        self.calls += 1
        return -100., 0.1

    def liquid_water(self):
        self.calls += 1
        return -200.


def reaction():
    return {'species': {'left': [-1, 0, 'CCO'], 'right': [1, 0, 'CCO']}}


def score(profile, rx, backend):
    if profile == 'reduced':
        return score_reaction(rx, backend, log=lambda *_: None)
    return StateFunctionEstimator(backend).score_reaction(rx)


@pytest.mark.parametrize('profile', ['reduced', 'state_function'])
@pytest.mark.parametrize('requested', [
    {'pH': 7.5, 'ionic_strength_M': 0.12},  # native GPDH reference conditions
    {'temperature_K': 310.15},
    {'pMg': 3.0},
    {'quantity': 'standard_reaction_gibbs'},
])
def test_unsupported_request_has_no_energy_or_molecular_work(profile, requested):
    backend = Backend()
    rx = dict(reaction(), conditions=requested)
    result = score(profile, rx, backend)
    assert result['dG'] is None
    assert result.get('dG_unrounded') is None
    assert 'unsupported calculation conditions' in result['suspect']
    assert backend.calls == 0
    assert result['conditions'] == CONDITIONS


@pytest.mark.parametrize('profile', ['reduced', 'state_function'])
def test_matching_and_omitted_conditions_have_identical_physical_results(profile):
    rx = reaction()
    baseline = score(profile, rx, Backend())
    request = dict(rx, conditions=copy.deepcopy(CONDITIONS))
    before = copy.deepcopy(request)
    matching = score(profile, request, Backend())
    assert baseline['dG'] is not None
    assert matching['dG'] == baseline['dG']
    assert matching['U_samp'] == baseline['U_samp']
    assert request == before


@pytest.mark.parametrize('conditions', [None, {'pH': float('nan')},
    {'T_K': True}, {'T_K': '298.15'}, {'temperature_C': 25}])
def test_invalid_or_unknown_requirements_are_not_silently_accepted(conditions):
    assert calculation_condition_issues(dict(reaction(), conditions=conditions))


def test_matching_partial_publication_schema_and_conflicting_aliases():
    assert not calculation_condition_issues(dict(reaction(), conditions={
        'temperature_K': 298.15, 'pH': 7, 'ionic_strength_M': 0,
        'magnesium_binding': 'absent', 'solute_standard_state_M': 1,
        'water_standard_state': 'pure_liquid', 'quantity': 'standard_transformed_reaction_gibbs',
        'energy_units': 'kJ/mol'}))
    issues = calculation_condition_issues(dict(reaction(), conditions={
        'T_K': 298.15, 'temperature_K': 310.15}))
    assert len(issues) == 1 and '310.15' in issues[0]
