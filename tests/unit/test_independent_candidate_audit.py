"""Protect chemical overlap decisions used before experimental enrollment."""
import csv
import importlib.util
from pathlib import Path

import pytest

SOURCE = Path(__file__).resolve().parents[2]/'benchmark/audit_independent_candidates.py'
SPEC = importlib.util.spec_from_file_location('candidate_audit', SOURCE)
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


def test_partial_abstract_row_retains_concrete_terms_on_both_sides(tmp_path):
    path = tmp_path/'training.csv'
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=['description', 'reaction'])
        writer.writeheader()
        writer.writerow(dict(description='abstract partner', reaction='CN + *S = C=N + *[SH]'))
    domain = AUDIT.training_domain(path)
    assert {'CN', 'C=N'} <= set(domain['molecules'])
    assert len(domain['unresolved_terms']) == 2
    assert not domain['reactions']  # No complete equation fabricated from a subset.


def test_extent_direction_and_explicit_hydrogens_do_not_create_novelty():
    first = AUDIT.signature(AUDIT.equation_parts('CC(=O)OC + O = CC(=O)O + CO'))
    reversed_doubled = AUDIT.signature(AUDIT.equation_parts('2 CO + 2 CC(=O)[O-] = 2 O + 2 CC(=O)OC'))
    assert first == reversed_doubled
    assert AUDIT.molecular_identity('[H]C([H])([H])O[H]')[0] == 'CO'


def test_enantiomers_remain_distinct_in_the_operational_family_key():
    left = AUDIT.molecular_identity('C[C@H](N)c1ccccc1')[0]
    right = AUDIT.molecular_identity('C[C@@H](N)c1ccccc1')[0]
    assert left != right


def test_reference_equation_requires_correct_proton_reservoir():
    parts = AUDIT.equation_parts('CC(=O)OC + O = CC(=O)[O-] + CO')
    with pytest.raises(ValueError, match='not element/charge balanced'):
        AUDIT.balance(parts)
    AUDIT.balance(parts, n_proton=1)


def test_structural_screen_cannot_admit_experimental_reference(tmp_path):
    path = tmp_path/'training.csv'
    path.write_text('description,reaction\nethanol,CCO\n')
    domain = AUDIT.training_domain(path)
    candidate = dict(id='balanced_control', species=[
        dict(name='acetone', coefficient=-1, smiles='CC(=O)C'),
        dict(name='hydrogen', coefficient=-1, smiles='[H][H]'),
        dict(name='isopropanol', coefficient=1, smiles='CC(O)C')])
    report = AUDIT.audit_candidate(candidate, domain)
    assert report['chemically_balanced']
    assert report['primary_benchmark_novelty_certified'] is False
