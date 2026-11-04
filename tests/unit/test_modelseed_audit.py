"""Protect literal structural accounting from false passes and silent curation."""
import importlib.util
from pathlib import Path
from fractions import Fraction

import pytest

BENCHMARK = Path(__file__).resolve().parents[2]/'benchmark'


def load(name):
    spec = importlib.util.spec_from_file_location(name, BENCHMARK/(name+'.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


A = load('audit_modelseed_inputs')
V = load('verify_modelseed_inputs')


def test_exact_stoichiometry_preserves_compartments_and_cancels_duplicates():
    text = '-0.5:cpdA:0:name;0.5:cpdA:0:name;-1:cpdB:0:B;1:cpdB:1:B'
    expected = {('cpdB', '0'): Fraction(-1), ('cpdB', '1'): Fraction(1)}
    assert A.parse_equation(text) == V.terms(text) == expected


@pytest.mark.parametrize('text', ['-1:cpdA', '-1::0:name', '-1:cpdA::name', 'bad:cpdA:0:name'])
def test_incomplete_or_nonnumeric_terms_fail_closed(text):
    with pytest.raises(ValueError):
        A.parse_equation(text)
    with pytest.raises(ValueError):
        V.terms(text)


@pytest.mark.parametrize('smiles', ['', ' ', '*', 'CC*'])
def test_empty_or_abstract_graph_does_not_become_a_known_zero_formula(smiles):
    assert V.graph(smiles) is None


def test_negative_element_residual_is_not_dropped():
    compounds = {'A': dict(charge='0', smiles='O'), 'B': dict(charge='0', smiles='[H][H]')}
    result = V.balance(V.terms('-1:A:0:A;1:B:0:B'), compounds)
    assert result['element_residual'] == {'O': '-1'}
    assert V.reservoir_completion(result) is None  # Cannot add oxygen alone as water.


def test_isotope_defects_are_not_repaired_with_ordinary_water_or_protons():
    compounds = {'A': dict(charge='0', smiles='[13CH3]O'), 'B': dict(charge='0', smiles='CO')}
    result = V.balance(V.terms('-1:A:0:A;1:B:0:B'), compounds)
    assert result['element_residual'] == {'13C': '-1', 'C': '1'}
    assert V.reservoir_completion(result) is None


def test_reference_conversion_checks_proton_sign_and_charge():
    source = dict(id='acetate', charge='-1', smiles='CC(=O)[O-]')
    row = dict(charge=-1, canonical_input='CC(=O)[O-]', reference_smiles='CC(=O)O',
               reference_charge=0, proton_exchange=-1, n_atoms=7, n_heavy_atoms=4)
    V.reference_bookkeeping(source, row)
    with pytest.raises(AssertionError):
        V.reference_bookkeeping(source, dict(row, proton_exchange=1))


def test_explicit_proton_and_water_completions_are_exact_not_energy_fits():
    proton = dict(all_formulas_known=True, element_residual={'H': '1'}, formal_charge_residual='1')
    assert V.reservoir_completion(proton) == dict(delta_water='0', delta_proton='-1', kind='proton_only')
    water = dict(all_formulas_known=True, element_residual={'H': '-2', 'O': '-1'}, formal_charge_residual='0')
    assert V.reservoir_completion(water) == dict(delta_water='1', delta_proton='0', kind='water_and_proton_basis')
    assert V.reservoir_completion(dict(water, formal_charge_residual='1')) is None


def test_balanced_or_unknown_equations_do_not_get_a_completion():
    assert V.reservoir_completion(dict(all_formulas_known=False)) is None
    assert V.reservoir_completion(dict(all_formulas_known=True, element_residual={}, formal_charge_residual='0')) is None
