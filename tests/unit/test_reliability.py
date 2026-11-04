"""Scientific invariants of equal-coverage reliability evaluation."""
import importlib.util
import itertools
import math
from pathlib import Path

import pytest

SOURCE = Path(__file__).resolve().parents[2]/'benchmark/reliability.py'
SPEC = importlib.util.spec_from_file_location('reliability', SOURCE)
R = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(R)


def rows():
    return [dict(reaction=str(i), reference=ref, comparator=pred, metag=meta, spread=spread)
            for i, (ref, pred, meta, spread) in enumerate([
                (1., 0., 2., .1), (5., 9., 8., .2), (30., 0., 20., .2), (3., 2., 5., .3)])]


def test_boundary_ties_equal_exhaustive_exact_size_subset_average():
    data = rows()
    actual = R.curve(data, [0, 1, 1, 2], 20.)[1]
    subsets = [[data[0], *choice] for choice in itertools.combinations(data[1:3], 1)]
    expected = math.fsum(math.fsum(R.error_values(r)['comparator'] for r in sub)/2
                         for sub in subsets)/len(subsets)
    assert actual['retained'] == 2 and actual['MAE']['comparator'] == expected
    assert actual['expected_severe_count']['comparator'] == .5


def test_fully_tied_ranking_has_same_expected_risk_at_every_coverage():
    data = rows()
    curve = R.curve(data, [0]*len(data), 20.)
    errors = math.fsum(R.error_values(r)['comparator'] for r in data)/len(data)
    assert all(p['MAE']['comparator'] == errors for p in curve)
    assert all(p['severe_rate']['comparator'] == .25 for p in curve)


def test_scores_are_label_free_and_permutation_invariant():
    data = rows()
    changed = [dict(r, reference=1000*i) for i, r in enumerate(data)]
    assert R.policy_scores(data) == R.policy_scores(changed)
    packet = R.evaluate(dict(rows=data))
    reverse = R.evaluate(dict(rows=list(reversed(data))))
    assert packet['curves'] == reverse['curves']
    for policy in R.POLICIES:
        assert packet['curves'][policy][-1]['MAE'] == packet['curves']['spread'][-1]['MAE']


def test_direction_reversal_and_units_preserve_rank_and_severe_events():
    data = rows()
    old = R.evaluate(dict(rows=data, severe_threshold_kj=20.))
    scaled = [dict(r, reference=-3*r['reference'], comparator=-3*r['comparator'],
                   metag=-3*r['metag'], spread=3*r['spread']) for r in data]
    new = R.evaluate(dict(rows=scaled, severe_threshold_kj=60.))
    for policy in R.POLICIES:
        assert new['comparator_AURC'][policy] == pytest.approx(3*old['comparator_AURC'][policy])
        for a, b in zip(old['curves'][policy], new['curves'][policy]):
            assert a['coverage'] == b['coverage']
            assert a['expected_severe_count'] == b['expected_severe_count']


def test_no_invented_intervals_without_dependence_provenance():
    assert R.evaluate(dict(rows=rows()))['bootstrap']['status'] == 'not_estimated'
    data = [dict(r, cluster='shared_source') for r in rows()]
    result = R.evaluate(dict(rows=data, cluster_definition='One shared primary source'))
    assert result['bootstrap']['status'] == 'inconclusive'
    assert 'paired_AURC_minus_spread_95_percentile' not in result['bootstrap']


def test_cluster_bootstrap_is_paired_and_reproducible():
    data = [dict(reaction=str(i), reference=float(i), comparator=float(i*i+1), metag=float(i*i+i+2),
                 spread=float(i), cluster=str(i)) for i in range(10)]
    packet = dict(rows=data, cluster_definition='Synthetic independent clusters for arithmetic test')
    a = R.evaluate(packet, bootstrap_replicates=30, seed=42)
    b = R.evaluate(packet, bootstrap_replicates=30, seed=42)
    assert a['bootstrap'] == b['bootstrap']
    for bounds in a['bootstrap']['paired_AURC_minus_spread_95_percentile'].values():
        assert bounds == [0., 0.]
    for policies in a['bootstrap']['paired_MAE_minus_spread_95_percentile'].values():
        assert all(bounds == [0., 0.] for bounds in policies.values())


def test_random_subset_zero_event_probability_is_exact_finite_population():
    report = R.evaluate(dict(rows=rows()))
    assert report['full_cohort_comparator_severe_count'] == 1
    assert [p['probability_zero_severe'] for p in report['random_subset_context']] == [.75, .5, .25, 0.]


@pytest.mark.parametrize('field,value', [('spread', -1.), ('reference', float('nan')),
                                       ('comparator', float('inf')), ('metag', True)])
def test_invalid_inputs_fail_closed(field, value):
    data = rows(); data[0][field] = value
    with pytest.raises(ValueError):
        R.evaluate(dict(rows=data))


def test_duplicate_reactions_and_incomplete_uncertainty_are_rejected():
    with pytest.raises(ValueError, match='distinct'):
        R.evaluate(dict(rows=rows()+rows()[:1]))
    data = rows(); data[0].pop('spread')
    with pytest.raises(KeyError):
        R.evaluate(dict(rows=data))
