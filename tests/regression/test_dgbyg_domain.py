import importlib.util
import sys
from pathlib import Path

import pytest

BENCHMARK = Path(__file__).resolve().parents[2] / 'benchmark'
sys.path.insert(0, str(BENCHMARK))
try:
    spec = importlib.util.spec_from_file_location('analyze_dgbyg_domain', BENCHMARK / 'analyze_dgbyg_domain.py')
    domain = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(domain)
finally:
    sys.path.pop(0)


def test_direction_and_water_do_not_change_training_reaction_match():
    train = domain.training_domain([{'reaction': 'CCO = CC=O'}])
    reaction = {'species': {'ethanol': [1, 0, 'CCO'], 'acetaldehyde': [-1, 0, 'CC=O'],
                            'water': [-1, 0, 'O']}}
    assert domain.classify(reaction, train)['tier'] == 'reaction_seen'


def test_unseen_compounds_and_unparseable_structures_are_distinguished():
    train = domain.training_domain([{'reaction': 'CCO = CC=O'}])
    novel = {'species': {'a': [-1, 0, 'CCO'], 'b': [1, 0, 'CCCCCCCC']}}
    assert domain.classify(novel, train)['n_unseen_compounds'] == 1
    assert domain.classify({'species': {'x': [1, 0, '*']}}, train)['tier'] == 'unclassifiable'


def test_reaction_change_similarity_is_direction_invariant():
    forward = domain.reaction_signature([(-1, 'CCO'), (1, 'CC=O')])
    reverse = domain.reaction_signature([(1, 'CCO'), (-1, 'CC=O')])
    training = [domain.reaction_change(forward)]
    assert domain.nearest_change(forward, training) == 1
    assert domain.nearest_change(reverse, training) == 1


def test_stereo_erasure_is_indeterminate_not_a_novel_reaction():
    train = domain.training_domain([{'reaction': 'CCO = CC=O'}])
    rx = {'species': {'a': [-1, 0, 'C[C@H](O)C(=O)O'], 'b': [1, 0, 'C[C@@H](O)C(=O)O']}}
    result = domain.classify(rx, train)
    assert result['tier'] == 'stereochemistry_collapsed'
    assert result['reaction_stereofree'] is None


def test_chiral_change_detects_epimerization_and_is_direction_invariant():
    pairs = [(-1, 'C[C@H](O)C(=O)O'), (1, 'C[C@@H](O)C(=O)O')]
    signature = domain.reaction_signature(pairs)
    assert not domain.reaction_change(signature)
    change = domain.reaction_change(signature, chiral=True)
    assert change
    assert domain.nearest_change(signature, [change], chiral=True) == 1
    reversed_signature = domain.reaction_signature([(-c, s) for c, s in pairs])
    assert domain.nearest_change(reversed_signature, [change], chiral=True) == 1


def test_comparator_conditions_and_duplicates_are_checked(tmp_path):
    f = tmp_path / 'predictions.tsv'
    header = 'method\tstatus\treaction\tT_K\tpH\tI_M\tpMg\tdG_kJ\n'
    row = 'dGbyG_eQ_bridge\tok\tr\t298.15\t7\t0\t14\t2\n'
    f.write_text(header + row)
    assert domain.read_comparator(f) == {'r': 2}
    f.write_text(header + row + row)
    with pytest.raises(ValueError, match='duplicate'):
        domain.read_comparator(f)
    f.write_text(header + row.replace('\t0\t14\t', '\t0.25\t14\t'))
    with pytest.raises(ValueError, match='conditions'):
        domain.read_comparator(f)


def test_nad_flag_uses_species_labels_only():
    assert domain.has_named_nad({'species': {'NADPH': [-1, -4, 'example']}})
    assert not domain.has_named_nad({'species': {'NADPH analog': [-1, -4, 'example']}})
