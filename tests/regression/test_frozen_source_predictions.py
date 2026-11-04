"""Frozen output coverage and scoring must remain consistent with the reference."""
import importlib.util
import json
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[2]
spec=importlib.util.spec_from_file_location('metag_benchmark_score',ROOT/'benchmark/score.py')
score=importlib.util.module_from_spec(spec);spec.loader.exec_module(score)
import sys
sys.path.insert(0, str(ROOT/'benchmark'))
import score_checked

def test_frozen_reaction_coverage_and_primary_reference_score():
    path=ROOT/'benchmark/metag_frozen_source_predictions.json';document=json.loads(path.read_text())
    rows=document['predictions'];assert len(rows)==364
    assert sum(r['status']=='estimated' for r in rows.values())==327
    assert sum(r['status']=='workflow_no_estimate' and r['dG'] is None for r in rows.values())==37
    assert all('exp' not in r and 'err' not in r for r in rows.values())
    labels=json.loads((ROOT/'benchmark/reference_mg_bounded.json').read_text())['2.0']['labels']
    predictions=score.load_predictions(path)
    pairs=[(predictions[r],x['exp']) for r,x in labels.items() if predictions.get(r) is not None]
    assert score.metrics(pairs)['n']==245
    assert score.metrics(pairs)['MAE']==8.68


@pytest.mark.parametrize("bad", [float("nan"), float("inf"), -float("inf"), True, "1.0"])
def test_scoring_rejects_invalid_prediction_values(tmp_path, bad):
    path = tmp_path / "predictions.json"
    path.write_text(json.dumps({"predictions": {"reaction": {"dG": bad}}}))
    with pytest.raises(ValueError, match="finite numeric energy"):
        score_checked.load_predictions(path)


def test_scoring_rejects_duplicate_list_predictions(tmp_path):
    path = tmp_path / "predictions.json"
    path.write_text(json.dumps([{"reaction": "r", "dG": 1}, {"reaction": "r", "dG": 2}]))
    with pytest.raises(ValueError, match="Duplicate prediction"):
        score_checked.load_predictions(path)


def test_scoring_keeps_missing_prediction(tmp_path):
    path = tmp_path / "predictions.json"
    path.write_text(json.dumps([{"reaction": "r", "dG": None}]))
    assert score_checked.load_predictions(path) == {"r": None}


def test_checked_scorer_matches_archived_scorer_on_frozen_panel(capsys):
    path = str(ROOT / 'benchmark/metag_frozen_source_predictions.json')
    original = score.main([path])
    checked = score_checked.main([path])
    assert checked == original
    capsys.readouterr()
