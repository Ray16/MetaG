"""Compare label-free reliability rankings at identical retained target counts.

The tool evaluates supplied predictions; it neither loads models nor certifies
experimental independence or the meaning/calibration of a supplied spread.
Boundary ties are averaged over exact-size subsets, never broken using errors.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path
import random

POLICIES = ('spread', 'disagreement', 'rank_sum', 'prediction_magnitude')
PREDICTORS = ('comparator', 'metag', 'zero')


def finite_number(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f'{name} must be a finite number')
    return float(value)


def validate_rows(rows):
    if not isinstance(rows, list) or not rows:
        raise ValueError('A nonempty reaction list is required')
    ids, clean = set(), []
    for row in rows:
        rid = row['reaction']
        if not isinstance(rid, str) or not rid or rid in ids:
            raise ValueError('Reaction IDs must be distinct nonempty strings')
        ids.add(rid)
        value = dict(reaction=rid)
        for name in ('reference', 'metag', 'comparator', 'spread'):
            value[name] = finite_number(row[name], name)
        if value['spread'] < 0:
            raise ValueError('A spread cannot be negative')
        if 'cluster' in row:
            if not isinstance(row['cluster'], str) or not row['cluster']:
                raise ValueError('A declared cluster must be a nonempty string')
            value['cluster'] = row['cluster']
        clean.append(value)
    return clean


def doubled_midranks(values):
    """Integer twice-midrank, invariant to monotone changes and boundary ordering."""
    groups = defaultdict(list)
    for i, value in enumerate(values):
        groups[value].append(i)
    ranks = [None]*len(values)
    before = 0
    for score in sorted(groups):
        indices = groups[score]
        rank = 2*before + len(indices) + 1
        for i in indices:
            ranks[i] = rank
        before += len(indices)
    return ranks


def policy_scores(rows):
    spread = [r['spread'] for r in rows]
    disagreement = [abs(r['metag']-r['comparator']) for r in rows]
    if not all(math.isfinite(value) for value in disagreement):
        raise ValueError('Nonfinite derived disagreement')
    a, b = doubled_midranks(spread), doubled_midranks(disagreement)
    return dict(spread=spread, disagreement=disagreement,
                rank_sum=[x+y for x, y in zip(a, b)],
                prediction_magnitude=[abs(r['comparator']) for r in rows])


def error_values(row):
    values = dict(comparator=abs(row['comparator']-row['reference']),
                  metag=abs(row['metag']-row['reference']), zero=abs(row['reference']))
    if not all(math.isfinite(value) for value in values.values()):
        raise ValueError('Nonfinite derived error')
    return values


def curve(rows, scores, threshold):
    """Expected metrics for uniformly selected exact-k boundary-tier subsets."""
    if len(rows) != len(scores):
        raise ValueError('Ranking and cohort lengths differ')
    groups = defaultdict(list)
    for row, score in zip(rows, scores):
        groups[score].append(row)
    before, totals, severe = 0, dict.fromkeys(PREDICTORS, 0.), dict.fromkeys(PREDICTORS, 0.)
    result = []
    for score in sorted(groups):
        tier = groups[score]
        errors = [error_values(r) for r in tier]
        sums = {name: math.fsum(e[name] for e in errors) for name in PREDICTORS}
        counts = {name: sum(e[name] > threshold for e in errors) for name in PREDICTORS}
        for take in range(1, len(tier)+1):
            retained = before+take
            fraction = take/len(tier)
            result.append(dict(retained=retained, coverage=retained/len(rows),
                boundary_tier_size=len(tier), boundary_tier_selected=take,
                MAE={name: (totals[name]+fraction*sums[name])/retained for name in PREDICTORS},
                expected_severe_count={name: severe[name]+fraction*counts[name] for name in PREDICTORS},
                severe_rate={name: (severe[name]+fraction*counts[name])/retained for name in PREDICTORS}))
        for name in PREDICTORS:
            totals[name] += sums[name]
            severe[name] += counts[name]
        before += len(tier)
    return result


def percentile(values, fraction):
    ordered = sorted(values)
    position = (len(ordered)-1)*fraction
    lo, hi = math.floor(position), math.ceil(position)
    return ordered[lo]+(ordered[hi]-ordered[lo])*(position-lo)


def cluster_intervals(rows, scores, *, threshold, replicates, seed, definition):
    """Conditional paired intervals; fixed feature scores are resampled together."""
    if not definition or any('cluster' not in row for row in rows):
        return dict(status='not_estimated', reason='Complete, justified source/family clusters were not supplied.')
    groups = defaultdict(list)
    for i, row in enumerate(rows):
        groups[row['cluster']].append(i)
    if len(groups) < 10:
        return dict(status='inconclusive', n_clusters=len(groups), reason='Fewer than ten declared clusters.')
    if not isinstance(replicates, int) or isinstance(replicates, bool) or replicates < 2:
        raise ValueError('At least two bootstrap replicates are required')
    rng, keys = random.Random(seed), sorted(groups)
    outcomes, point_outcomes = defaultdict(list), defaultdict(list)
    for _ in range(replicates):
        indices = [i for key in rng.choices(keys, k=len(keys)) for i in groups[key]]
        selected = [rows[i] for i in indices]
        risks, points = {}, {}
        for name in POLICIES:
            values = curve(selected, [scores[name][i] for i in indices], threshold)
            risks[name] = math.fsum(v['MAE']['comparator'] for v in values)/len(values)
            points[name] = {str(f): values[math.ceil(len(values)*f)-1]['MAE']['comparator']
                            for f in (.25, .5, .75, 1.)}
        for name in POLICIES[1:]:
            outcomes[name].append(risks[name]-risks['spread'])
            for f in points[name]:
                point_outcomes[(name, f)].append(points[name][f]-points['spread'][f])
    return dict(status='conditional_intervals', n_clusters=len(groups), cluster_definition=definition,
        n_replicates=replicates, seed=seed,
        paired_AURC_minus_spread_95_percentile={name: [percentile(v, .025), percentile(v, .975)]
                                               for name, v in outcomes.items()},
        paired_MAE_minus_spread_95_percentile={name: {f: [percentile(point_outcomes[(name, f)], .025),
                                                         percentile(point_outcomes[(name, f)], .975)]
                                                      for f in ('0.25', '0.5', '0.75', '1.0')}
                                               for name in POLICIES[1:]},
        qualification='Conditional on supplied clusters and fixed scores. Does not certify source independence, adequate power, total predictive uncertainty or experimental validity.')


def evaluate(packet, *, bootstrap_replicates=10000, seed=20261104):
    rows = validate_rows(packet['rows'])
    threshold = finite_number(packet.get('severe_threshold_kj', 20.), 'severe threshold')
    if threshold < 0:
        raise ValueError('Severe threshold cannot be negative')
    scores = policy_scores(rows)
    curves = {name: curve(rows, values, threshold) for name, values in scores.items()}
    n = len(rows)
    severe = sum(error_values(r)['comparator'] > threshold for r in rows)
    highlights, aurc = {}, {}
    for name, values in curves.items():
        highlights[name] = {str(f): values[math.ceil(n*f)-1] for f in (.25, .5, .75, 1.)}
        aurc[name] = math.fsum(v['MAE']['comparator'] for v in values)/n
    return dict(schema='metag-equal-coverage-reliability-v1', n=n,
        declared_cohort_exposure=packet.get('cohort_exposure', 'unspecified'),
        supplied_spread_definition=packet.get('spread_definition', 'unspecified'),
        severe_threshold_kj=threshold, full_cohort_comparator_severe_count=severe,
        ranking='Ascending feature score; rank_sum is equal-weight sum of doubled empirical midranks.',
        ties='Expected risk over all exact-k boundary-tier subsets. No use of labels to break ties.',
        curves=curves, highlights=highlights, comparator_AURC=aurc,
        comparator_AURC_minus_spread={name: value-aurc['spread'] for name, value in aurc.items()},
        random_subset_context=[dict(retained=k, probability_zero_severe=(math.comb(n-severe, k)/math.comb(n, k)
                                  if k <= n-severe else 0.)) for k in range(1, n+1)],
        bootstrap=cluster_intervals(rows, scores, threshold=threshold, replicates=bootstrap_replicates,
                                    seed=seed, definition=packet.get('cluster_definition')),
        scope='Finite supplied cohort. Generalization, total error calibration and source/experimental independence require separate evidence. Rankings are compared using the same retained counts; no new predictor is selected or fitted.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--bootstrap-replicates', type=int, default=10000)
    parser.add_argument('--seed', type=int, default=20261104)
    args = parser.parse_args()
    packet = json.loads(args.input.read_text())
    report = evaluate(packet, bootstrap_replicates=args.bootstrap_replicates, seed=args.seed)
    report['input_sha256'] = hashlib.sha256(args.input.read_bytes()).hexdigest()
    report['analysis_source_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open('x') as stream:
        json.dump(report, stream, indent=2, allow_nan=False); stream.write('\n')
    print(json.dumps(dict(n=report['n'], comparator_AURC=report['comparator_AURC'],
                         severe=report['full_cohort_comparator_severe_count'], bootstrap=report['bootstrap']), indent=2))


if __name__ == '__main__':
    main()
