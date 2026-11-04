"""Guard exposure and training-exclusion decisions before target errors exist."""
import csv
import importlib.util
import json
from fractions import Fraction
from pathlib import Path
from types import SimpleNamespace

import pytest

SOURCE = Path(__file__).resolve().parents[2] / 'benchmark/benchmark_feasibility.py'
SPEC = importlib.util.spec_from_file_location('benchmark_feasibility', SOURCE)
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


@pytest.fixture(autouse=True)
def clear_normalization_cache():
    AUDIT.family.cache_clear()
    yield
    AUDIT.family.cache_clear()


def table(tmp_path, rows):
    path = tmp_path / 'training.csv'
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=[*AUDIT.ALLOWED_COLUMNS, 'dG', 'error'])
        writer.writeheader()
        for row in rows:
            writer.writerow({'type': 'reaction dg', 'dG': 'SECRET_TARGET',
                             'error': 'SECRET_ERROR', **row})
    return path


@pytest.mark.parametrize('first,second', [
    ('C[n+]1ccc(C(=O)[O-])cc1', 'C[n+]1ccc(C(=O)O)cc1'),
    ('[NH3+]CC(=O)[O-]', 'NCC(=O)O'),
    ('OCO', 'C=O'),
    ('CC(O)(O)C', 'CC(=O)C'),
])
def test_reference_aliases_reconcile_charges_and_hydrates(first, second):
    assert AUDIT.family(first)['pool'] == AUDIT.family(second)['pool']


def test_pool_retains_nonanomeric_stereochemistry_broad_exposure_does_not():
    left, right = map(AUDIT.family, ['C[C@H](O)C(=O)O', 'C[C@@H](O)C(=O)O'])
    assert left['pool'] != right['pool']
    assert left['broad'] == right['broad']


def test_sugar_ring_open_aliases_preserve_family_identity():
    opened = 'O=C[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO'
    states = AUDIT.enumerate_sugar_states(opened)
    assert states.complete and len(states.states) > 2
    assert {AUDIT.family(s['smiles'])['pool'] for s in states.states} == {opened}
    a, b = states.states[:2]
    glucose = AUDIT.reaction_key([(-1, a['smiles']), (1, b['smiles'])], 'pool')
    identity = AUDIT.reaction_key([(-1, 'CCO'), (1, 'CCO')], 'pool')
    assert glucose != identity


def test_rational_extent_and_direction_are_invariant():
    first = AUDIT.equation('CC(=O)OC + O = CC(=O)O + CO')
    reverse = AUDIT.equation('1/2 CO + 0.5 CC(=O)[O-] = 1/2 O + 1/2 CC(=O)OC')
    assert AUDIT.reaction_key(first, 'pool') == AUDIT.reaction_key(reverse, 'pool')
    assert reverse[0][0] == -Fraction(1, 2)
    with pytest.raises(ValueError):
        AUDIT.equation('1/0 CO = C=O')


def test_incomplete_tautomers_remain_unknown(monkeypatch):
    monkeypatch.setattr(AUDIT, 'TAUTOMERS', SimpleNamespace(
        Enumerate=lambda _: SimpleNamespace(status='limit reached')))
    info = AUDIT.family('CCO')
    assert info['broad'] is None
    assert 'tautomer enumeration incomplete' in info['issues']


def test_incomplete_sugar_alias_cannot_certify_an_absent_family(monkeypatch):
    monkeypatch.setattr(AUDIT, 'enumerate_sugar_states', lambda _: SimpleNamespace(
        applicable=True, complete=False, reason='unresolved state'))
    info = AUDIT.family('CCO')
    assert info['pool'] is None and info['broad'] is None


def test_isotope_radical_and_disconnected_terms_are_not_currency():
    for smiles in ['[2H]O[2H]', '[H]', 'CCO.[Na+]']:
        info = AUDIT.family(smiles)
        assert not info['reservoir'] and info['broad'] is None and info['issues']
    assert AUDIT.family('O')['reservoir'] and AUDIT.family('[H+]')['reservoir']


def test_partial_development_row_keeps_concrete_partner_and_unknown(tmp_path):
    panel = tmp_path / 'panel.json'
    panel.write_text(json.dumps({'partial': {'species': {
        'known': [-1, 0, 'CCO'], 'abstract': [1, 0, '*O']}}}))
    reactions, molecules, unresolved = AUDIT.development_panel(panel)
    assert 'CCO' in molecules['broad']
    assert reactions['partial']['pool'] is None
    assert unresolved[0]['species'] == 'abstract'


def test_projection_and_candidate_report_never_export_target_values(tmp_path):
    training = table(tmp_path, [dict(description='hydrolysis', reaction='CC(=O)OC + O = CC(=O)O + CO')])
    panel = tmp_path / 'panel.json'
    panel.write_text('{}')
    report = AUDIT.audit(training, panel)
    encoded = json.dumps(report)
    assert 'SECRET_TARGET' not in encoded and 'SECRET_ERROR' not in encoded
    assert set(AUDIT.projected_table(training)[0]) == set(AUDIT.ALLOWED_COLUMNS)
    assert report['experimental_references_admitted'] == 0
    assert report['independence_certified'] is False


def test_holdout_excludes_formation_aliases_reverse_equations_and_wildcards(tmp_path):
    training = table(tmp_path, [
        dict(description='test', reaction='CC(=O)OC + O = CC(=O)O + CO'),
        dict(description='reverse duplicate', reaction='2 CO + 2 CC(=O)[O-] = 2 CC(=O)OC + 2 O'),
        dict(description='formation', reaction='CC(=O)OC', type='formation dg'),
        dict(description='abstract partner', reaction='CCO + *O = CC=O + *[OH]'),
        dict(description='retained reaction', reaction='CCO = CC=O'),
        dict(description='retained formation', reaction='CCO', type='formation dg'),
    ])
    plan = AUDIT.holdout_plan(training, [0], ['CC(=O)OC'])
    assert plan['retained_rows'] == [4, 5]
    assert [r['row'] for r in plan['excluded_rows']] == [0, 1, 2, 3]
    assert 'held-out molecular family' in plan['excluded_rows'][2]['reasons']
    assert 'unresolved training row:' in plan['excluded_rows'][3]['reasons'][0]
    assert not plan['launch_allowed'] and not plan['model_training_performed']
    assert 'SECRET_TARGET' not in json.dumps(plan)


def test_holdout_family_must_occur_in_evaluation_equation(tmp_path):
    training = table(tmp_path, [dict(description='test', reaction='CCO = CC=O')])
    with pytest.raises(ValueError, match='does not participate'):
        AUDIT.holdout_plan(training, [0], ['c1ccccc1'])


def test_formation_targets_cannot_be_enrolled_as_reaction_tests(tmp_path):
    training = table(tmp_path, [dict(description='formation', reaction='CCO', type='formation dg')])
    with pytest.raises(ValueError, match='valid reaction target'):
        AUDIT.holdout_plan(training, [0], ['CCO'])


def test_distance_report_preserves_concrete_partners_in_excluded_rows(tmp_path):
    training = table(tmp_path, [
        dict(description='test', reaction='NCCc1ccccc1 = CC=O'),
        dict(description='abstract first', reaction='*O + OCCc1ccccc1 = *N + CC=O'),
        dict(description='retained formation', reaction='CO', type='formation dg'),
    ])
    plan = AUDIT.holdout_plan(training, [0], ['NCCc1ccccc1'])
    neighbour = plan['heldout_neighbours'][0]
    assert neighbour['nearest_retained_pool'] == 'CO'
    assert neighbour['original_other_family_nearest_pool'] == 'OCCc1ccccc1'
    assert neighbour['original_other_family_nearest_similarity'] > neighbour['nearest_similarity']


def test_named_matching_preserves_phase_stereo_and_primitive_stoichiometry():
    key = AUDIT.named_equation_key
    assert key('A(aq) + 2 B(aq) = C(aq)') == key('1/2 c(aq) = b(aq) + 1/2 a(aq)')
    assert key('A(aq) = B(aq)') != key('A(g) = B(aq)')
    assert key('D-A(aq) = B(aq)') != key('L-A(aq) = B(aq)')


def test_source_link_does_not_impute_missing_conditions_or_confirm_ancestry(tmp_path):
    source = tmp_path / 'source.csv'
    with source.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=[*AUDIT.SOURCE_COLUMNS, 'K_prime'])
        writer.writeheader()
        writer.writerow({'reaction': 'alcohol(aq) = aldehyde(aq)', 'reference': 'SOURCE',
                         'temperature': '298.15', 'p_h': '7', 'K_prime': 'SECRET_EQUILIBRIUM'})
    report = {'input_sha256': {'training': 'test'}, 'rows': [
        {'row': 0, 'description': 'aldehyde(aq) = alcohol(aq)',
         'status': 'structural_candidate_only', 'matches_development': {'exact': [], 'pool': []}}]}
    result = AUDIT.source_candidates(report, source)
    assert result['counts']['with_named_source_candidates'] == 1
    assert result['counts']['with_all_four_condition_fields_present'] == 0
    assert result['rows'][0]['ancestry_confirmed'] is False
    assert result['experimental_references_admitted'] == 0
    assert 'SECRET_EQUILIBRIUM' not in json.dumps(result)
