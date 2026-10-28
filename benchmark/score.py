"""Score MetaG predictions against the Mg-bounded openTECR reference labels.

    python benchmark/score.py                                  # the archived 0.2.0.dev24 predictions
    python benchmark/score.py metag_out/metag_all.json         # output of `metag score`

Labels: ``reference_mg_bounded.json`` (threshold "2.0" = primary policy, "1.0"/"4.0" = sensitivities;
see PREREGISTRATION.md). Every metric is reported together with the number of scored labels; reactions
without a MetaG estimate are counted, never imputed.
"""
import json
import statistics as st
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def load_predictions(path):
    data = json.loads(Path(path).read_text())
    if isinstance(data, list):                                  # `metag score` output: list of results
        return {r["reaction"]: r.get("dG") for r in data}
    return {k: v.get("dG") for k, v in data["predictions"].items()}


def metrics(pairs):
    if not pairs:
        return {"n": 0}
    err = [p - y for p, y in pairs]
    ae = [abs(e) for e in err]
    return {"n": len(ae), "MAE": round(st.mean(ae), 2), "median": round(st.median(ae), 2),
            "RMSE": round(st.mean(e * e for e in err) ** 0.5, 2), "bias": round(st.mean(err), 2)}


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    pred = load_predictions(argv[0] if argv else HERE / "metag_dev24_predictions.json")
    ref = json.loads((HERE / "reference_mg_bounded.json").read_text())
    report = {}
    for threshold, block in ref.items():
        labels = {r: v["exp"] for r, v in block["labels"].items()}
        scored = [(pred[r], y) for r, y in labels.items() if pred.get(r) is not None]
        report[threshold] = {"n_labels": len(labels), "n_without_estimate": len(labels) - len(scored),
                             "MetaG": metrics(scored),
                             "zero_predictor_same_labels": metrics([(0.0, y) for _, y in scored])}
    print(json.dumps(report, indent=1))
    return report


if __name__ == "__main__":
    main()
