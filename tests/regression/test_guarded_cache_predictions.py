"""Coverage and score of the guarded cache replay against its fixed reference."""
import importlib.util
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
spec=importlib.util.spec_from_file_location('metag_benchmark_score',ROOT/'benchmark/score.py')
score=importlib.util.module_from_spec(spec);spec.loader.exec_module(score)

def test_guarded_replay_coverage_and_frozen_common_values():
    rows=json.loads((ROOT/'benchmark/metag_guarded_cache_predictions.json').read_text())['predictions']
    old=json.loads((ROOT/'benchmark/metag_frozen_source_predictions.json').read_text())['predictions']
    assert len(rows)==364
    assert sum(r['status']=='estimated' for r in rows.values())==325
    assert sum(r['status']=='workflow_no_estimate' and r['dG'] is None for r in rows.values())==35
    assert sum(r['status']=='dependency_missing' and r['dG'] is None for r in rows.values())==4
    assert all(row['dG']==old[k]['dG'] for k,row in rows.items() if row['status']=='estimated')
    labels=json.loads((ROOT/'benchmark/reference_mg_bounded.json').read_text())['2.0']['labels']
    pairs=[(rows[k]['dG'],value['exp']) for k,value in labels.items() if rows[k]['status']=='estimated']
    assert score.metrics(pairs)['n']==243
    assert score.metrics(pairs)['MAE']==8.71
