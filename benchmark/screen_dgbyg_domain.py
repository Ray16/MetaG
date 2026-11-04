"""Screen reaction inputs against a specified dGbyG training table; no labels read."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

from rdkit import rdBase
from analyze_dgbyg_domain import (classify, generalized_jaccard, parse_training_equation,
                                  reaction_change, reaction_signature, training_domain)


def screen(reactions, training_rows):
    domain = training_domain(training_rows)
    indexed = []
    for i, row in enumerate(training_rows):
        if ' = ' not in row.get('reaction', ''):
            continue
        try:
            signature = reaction_signature(parse_training_equation(row['reaction']))
            change = reaction_change(signature, chiral=True)
            if change:
                indexed.append((i, row.get('description', ''), change))
        except ValueError:
            continue
    results = {}
    for rid, rx in reactions.items():
        info = classify(rx, domain)
        try:
            signature = reaction_signature([(v[0], v[2]) for v in rx['species'].values()])
            change = reaction_change(signature, chiral=True)
            if change and indexed:
                flipped = {(bit, not sign): count for (bit, sign), count in change.items()}
                scored = [(max(generalized_jaccard(change, t), generalized_jaccard(flipped, t)), i, name)
                          for i, name, t in indexed]
                best = max(x[0] for x in scored)
                info['nearest_chiral_training_rows'] = [
                    {'row_index_zero_based': i, 'description': name}
                    for score, i, name in scored if abs(score - best) < 1e-12]
            else:
                info['nearest_chiral_training_rows'] = []
        except (ValueError, KeyError, TypeError):
            info['nearest_chiral_training_rows'] = []
        results[rid] = info
    return results


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--training', type=Path, required=True)
    p.add_argument('--reactions', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args(argv)
    with args.training.open(newline='') as f:
        rows = list(csv.DictReader(f))
    report = {'rdkit_version': rdBase.rdkitVersion,
              'input_sha256': {name: hashlib.sha256(getattr(args, name).read_bytes()).hexdigest()
                               for name in ('training', 'reactions')},
              'scope': 'Heuristic structure and reaction-change screen; neither an exposure audit nor a prediction of model error.',
              'reactions': screen(json.loads(args.reactions.read_text()), rows)}
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
