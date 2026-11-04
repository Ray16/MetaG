"""Validate inputs and score predictions with the archived benchmark metrics.

The original score.py remains byte-identical to the version pinned in the
study protocol. This standard entry point rejects malformed input before
calling its metric function.
"""
import json
import math
import sys
from pathlib import Path
from score import HERE, metrics


def finite_energy(value, context):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"{context} must be a finite numeric energy")
    return value


def load_predictions(path):
    data = json.loads(Path(path).read_text())
    if isinstance(data, list):
        rows = {}
        for row in data:
            reaction = row["reaction"]
            if reaction in rows:
                raise ValueError(f"Duplicate prediction for {reaction}")
            rows[reaction] = row.get("dG")
    else:
        rows = {key: value.get("dG") for key, value in data["predictions"].items()}
    for reaction, value in rows.items():
        if value is not None:
            finite_energy(value, f"Prediction for {reaction}")
    return rows


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    pred = load_predictions(argv[0] if argv else HERE / 'metag_corrected_source_predictions.json')
    ref = json.loads((HERE / 'reference_mg_bounded.json').read_text())
    report = {}
    for threshold, block in ref.items():
        labels = {r: finite_energy(v['exp'], f'Reference for {r}') for r, v in block['labels'].items()}
        scored = [(pred[r], y) for r, y in labels.items() if pred.get(r) is not None]
        report[threshold] = {'n_labels': len(labels), 'n_without_estimate': len(labels) - len(scored),
                             'MetaG': metrics(scored),
                             'zero_predictor_same_labels': metrics([(0., y) for _, y in scored])}
    print(json.dumps(report, indent=1))
    return report


if __name__ == '__main__':
    main()
