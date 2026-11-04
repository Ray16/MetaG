"""Score archived dGbyG cross-validation predictions on the MetaG reference.

The input must be the original pairing table produced by compare_metag.py, not a
table already transformed to MetaG's orientation or conditions.
"""
import argparse
import csv
import hashlib
import json
import math
import statistics
from pathlib import Path

from score import load_predictions, metrics


HERE = Path(__file__).resolve().parent


def normalized_prediction(row):
    """Paper-scale I=.25 prediction -> MetaG reaction orientation at I=0.

    scale is paper_extent / MetaG_extent. The bridge is in MetaG's extent,
    G(I=0) - G(I=.25), and therefore must be added after dividing by scale.
    """
    value = float(row['dgbyg'])
    scale = float(row['scale'])
    bridge = float(row['bridge_delta_I0_minus_I025'])
    if not all(math.isfinite(x) for x in (value, scale, bridge)) or scale == 0:
        raise ValueError('Prediction, nonzero scale, and bridge must be finite')
    return value / scale + bridge


def evaluate(rows, labels, metag, *, identity=None):
    """Score common pairs, optionally admitting only audited exact-family matches.

    Without ``identity`` this is the broader exploratory comparison. When an
    identity mapping is supplied, every input row must have an audited status;
    the recorded paper row, signed extent and bridge are checked where available.
    """
    pairs, seen = [], set()
    identity_counts = {'exact_family_equation_match': 0,
                       'bridge_family_identity_or_extent_unresolved': 0}
    for row in rows:
        rid = row['reaction']
        if rid in seen:
            raise ValueError(f'duplicate reaction: {rid}')
        seen.add(rid)
        value = normalized_prediction(row)
        if identity is not None:
            if rid not in identity:
                raise ValueError(f'missing identity audit for reaction: {rid}')
            record = identity[rid]
            status = record.get('status')
            if status not in identity_counts:
                raise ValueError(f'unknown identity audit status for {rid}: {status!r}')
            identity_counts[status] += 1
            if 'author_paper_row' in record and int(row['paper_row']) != record['author_paper_row']:
                raise ValueError(f'paper-row mismatch for {rid}')
            for source_key, row_key, tolerance in (
                ('signed_author_scale', 'scale', 1e-12),
                ('bridge_I0_noMg_minus_native_kj', 'bridge_delta_I0_minus_I025', 1e-7),
            ):
                if source_key in record and not math.isclose(float(record[source_key]), float(row[row_key]),
                                                             rel_tol=0, abs_tol=tolerance):
                    raise ValueError(f'{source_key} mismatch for {rid}')
            if status != 'exact_family_equation_match':
                continue
            if (rid in labels and metag.get(rid) is not None
                    and ('signed_author_scale' not in record or 'bridge_I0_noMg_minus_native_kj' not in record)):
                raise ValueError(f'missing audited scale/bridge for scored reaction: {rid}')
        if rid in labels and metag.get(rid) is not None:
            prediction, reference = metag[rid], labels[rid]['exp']
            if any(isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x)
                   for x in (prediction, reference)):
                raise ValueError(f'nonfinite or invalid MetaG prediction/reference for {rid}')
            pairs.append((rid, prediction, value, reference))
    result = {
        'n_input_rows': len(seen),
        'n_common': len(pairs),
        'n_reversed_common': sum(float(row['scale']) < 0 for row in rows
                                 if row['reaction'] in labels and metag.get(row['reaction']) is not None
                                 and (identity is None or identity[row['reaction']]['status'] == 'exact_family_equation_match')),
        'paired_wins': {
            'MetaG': sum(abs(m - y) < abs(d - y) for _, m, d, y in pairs),
            'dGbyG': sum(abs(d - y) < abs(m - y) for _, m, d, y in pairs),
            'ties': sum(abs(m - y) == abs(d - y) for _, m, d, y in pairs),
        },
        'MetaG': metrics([(m, y) for _, m, _, y in pairs]),
        'dGbyG_oof_harmonized': metrics([(d, y) for _, _, d, y in pairs]),
        'zero': metrics([(0.0, y) for _, _, _, y in pairs]),
        'pairs': [{'reaction': rid, 'MetaG': m, 'dGbyG': d, 'reference': y}
                  for rid, m, d, y in pairs],
    }
    if identity is not None:
        result['identity_status_counts'] = identity_counts
    return result


def ec_class(reaction, reactions):
    """Return a single top-level EC class, or unresolved for mixed/missing labels."""
    ec = str(reactions.get(reaction, {}).get('EC', ''))
    top = {item.strip().split('.')[0] for item in ec.replace('&', ';').split(';') if item.strip()}
    return next(iter(top)) if len(top) == 1 and next(iter(top)) in set('1234567') else 'unresolved'


def by_ec_class(pairs, reactions):
    """Exploratory descriptive summaries; labels may contain multiple EC numbers."""
    classes = {}
    for pair in pairs:
        key = ec_class(pair['reaction'], reactions)
        classes.setdefault(key, []).append(pair)
    return {key: {
        'MetaG': metrics([(p['MetaG'], p['reference']) for p in group]),
        'dGbyG': metrics([(p['dGbyG'], p['reference']) for p in group]),
        'MetaG_wins': sum(abs(p['MetaG'] - p['reference']) < abs(p['dGbyG'] - p['reference']) for p in group),
    } for key, group in sorted(classes.items())}


def within_class_correlation(pairs, reactions, field):
    """Pearson correlation after centering prediction and reference within EC classes."""
    groups = {}
    for row in pairs:
        groups.setdefault(ec_class(row['reaction'], reactions), []).append(row)
    centered = [(row[field] - statistics.mean(x[field] for x in group),
                 row['reference'] - statistics.mean(x['reference'] for x in group))
                for group in groups.values() for row in group]
    xx = sum(x*x for x, _ in centered)
    yy = sum(y*y for _, y in centered)
    if xx <= 0 or yy <= 0:
        return None
    return sum(x*y for x, y in centered) / math.sqrt(xx*yy)


def load_identity(path):
    """Validate the retained 164-pair identity record before using its mask."""
    document = json.loads(Path(path).read_text())
    rows = document['rows']
    if not isinstance(rows, dict) or document.get('n_requested') != len(rows):
        raise ValueError('Identity record has an inconsistent reaction count')
    counts = {'exact_family_equation_match': 0, 'bridge_family_identity_or_extent_unresolved': 0}
    author_rows = set()
    for rid, record in rows.items():
        if not isinstance(record, dict):
            raise ValueError(f'invalid identity audit record for {rid}')
        status = record.get('status')
        if status not in counts:
            raise ValueError(f'unknown identity audit status for {rid}: {status!r}')
        counts[status] += 1
        paper_row = record.get('author_paper_row')
        if isinstance(paper_row, bool) or not isinstance(paper_row, int) or paper_row < 0:
            raise ValueError(f'missing or invalid author paper row for {rid}')
        if paper_row in author_rows:
            raise ValueError(f'duplicate author paper row: {paper_row}')
        author_rows.add(paper_row)
        has_scale = 'signed_author_scale' in record
        has_bridge = 'bridge_I0_noMg_minus_native_kj' in record
        if has_scale != has_bridge:
            raise ValueError(f'incomplete scale/bridge identity record for {rid}')
        if has_scale:
            scale, bridge = record['signed_author_scale'], record['bridge_I0_noMg_minus_native_kj']
            if any(isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x)
                   for x in (scale, bridge)) or scale == 0:
                raise ValueError(f'invalid scale/bridge identity record for {rid}')
    if document.get('status_counts') != counts:
        raise ValueError('Identity status totals do not match its rows')
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pairs', type=Path, required=True)
    parser.add_argument('--reference', type=Path, default=HERE / 'reference_mg_bounded.json')
    parser.add_argument('--predictions', type=Path, default=HERE / 'metag_corrected_source_predictions.json',
                        help='MetaG predictions (default: completed corrected-source panel reported in the manuscript)')
    parser.add_argument('--identity', type=Path, default=HERE / 'dgbyg_family_identity_complete.json',
                        help='audited reaction-family identity and bridge metadata')
    parser.add_argument('--reactions', type=Path, default=HERE.parent / 'src/metag/data/reactions_opentecr_std.json')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(argv)
    with args.pairs.open(newline='') as stream:
        rows = list(csv.DictReader(stream, delimiter='\t'))
    labels = json.loads(args.reference.read_text())['2.0']['labels']
    identity = load_identity(args.identity)
    predictions = load_predictions(args.predictions)
    report = evaluate(rows, labels, predictions, identity=identity)
    exploratory = evaluate(rows, labels, predictions)
    report['n_common_including_unresolved_identity'] = exploratory['n_common']
    report['n_common_excluded_unresolved_identity'] = exploratory['n_common'] - report['n_common']
    reactions = json.loads(args.reactions.read_text())
    report['exploratory_by_ec_class'] = by_ec_class(report['pairs'], reactions)
    report['within_class_correlation'] = {
        'MetaG': within_class_correlation(report['pairs'], reactions, 'MetaG'),
        'dGbyG': within_class_correlation(report['pairs'], reactions, 'dGbyG'),
        'scope': 'Descriptive EC-class-centered Pearson correlation, not an independent validation.'}
    report['input_sha256'] = {name: hashlib.sha256(path.read_bytes()).hexdigest()
                              for name, path in [('pairs', args.pairs), ('reference', args.reference),
                                                 ('predictions', args.predictions), ('reactions', args.reactions),
                                                 ('identity', args.identity)]}
    report['limitations'] = ('Development-exposed MetaG and reference; dGbyG out-of-fold only in '
                             'its original reaction split; shared dissociation bridge; archived '
                             'exact-family admission does not propagate bridge or reference uncertainty; '
                             'archived reaction mapping and ensemble require independent audit.')
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    print(json.dumps({k: v for k, v in report.items() if k != 'pairs'}, indent=2))


if __name__ == '__main__':
    main()
