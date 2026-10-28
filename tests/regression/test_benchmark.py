"""The shipped benchmark labels and archived predictions reproduce the reported headline numbers."""
import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def _score():
    spec = importlib.util.spec_from_file_location("benchmark_score", ROOT / "benchmark" / "score.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.main([])


def test_primary_policy_headline(capsys):
    report = _score()["2.0"]
    assert report["n_labels"] == 274
    assert report["MetaG"]["n"] == 239
    assert report["MetaG"]["MAE"] == 8.44
    assert report["zero_predictor_same_labels"]["MAE"] == 11.68


def test_threshold_sensitivities(capsys):
    report = _score()
    assert {k: report[k]["n_labels"] for k in report} == {"2.0": 274, "1.0": 262, "4.0": 284}
    assert report["1.0"]["MetaG"]["MAE"] == 8.3
