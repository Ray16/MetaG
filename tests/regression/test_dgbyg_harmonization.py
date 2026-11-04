import importlib.util
import sys
from pathlib import Path

import pytest

BENCHMARK = Path(__file__).resolve().parents[2] / 'benchmark'
sys.path.insert(0, str(BENCHMARK))
try:
    spec = importlib.util.spec_from_file_location('score_dgbyg_oof', BENCHMARK / 'score_dgbyg_oof.py')
    score = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(score)
finally:
    sys.path.pop(0)


def test_paper_extent_and_ionic_bridge_are_applied_in_the_right_order():
    row = {'reaction': 'r', 'dgbyg': '-14', 'scale': '-2', 'bridge_delta_I0_minus_I025': '3'}
    assert score.normalized_prediction(row) == 10
    report = score.evaluate([row], {'r': {'exp': 9}}, {'r': 11})
    assert report['n_reversed_common'] == 1
    assert report['dGbyG_oof_harmonized']['MAE'] == 1
    assert report['MetaG']['MAE'] == 2


def test_duplicate_and_invalid_predictions_fail_closed():
    row = {'reaction': 'r', 'dgbyg': '1', 'scale': '1', 'bridge_delta_I0_minus_I025': '0'}
    with pytest.raises(ValueError, match='duplicate'):
        score.evaluate([row, row], {}, {})
    with pytest.raises(ValueError, match='nonzero'):
        score.normalized_prediction(dict(row, scale='0'))


def test_ec_grouping_does_not_force_cross_class_assignments():
    pairs = [{'reaction': 'a', 'MetaG': 0, 'dGbyG': 2, 'reference': 1},
             {'reaction': 'b', 'MetaG': 0, 'dGbyG': 1, 'reference': 1}]
    result = score.by_ec_class(pairs, {'a': {'EC': '2.1.1.1;3.1.1.1'}, 'b': {'EC': '3.1.1.1'}})
    assert result['unresolved']['MetaG']['n'] == 1
    assert result['3']['dGbyG']['MAE'] == 0


def test_exact_family_filter_and_recorded_bridge_checks():
    rows = [
        {'reaction':'a','paper_row':'10','dgbyg':'-14','scale':'-2','bridge_delta_I0_minus_I025':'3'},
        {'reaction':'b','paper_row':'11','dgbyg':'3','scale':'1','bridge_delta_I0_minus_I025':'0'},
    ]
    identity = {
        'a': {'status':'exact_family_equation_match','author_paper_row':10,
              'signed_author_scale':-2.,'bridge_I0_noMg_minus_native_kj':3.},
        'b': {'status':'bridge_family_identity_or_extent_unresolved','author_paper_row':11},
    }
    labels = {'a':{'exp':9},'b':{'exp':1}}
    predictions = {'a':11,'b':2}
    exact = score.evaluate(rows, labels, predictions, identity=identity)
    assert exact['n_input_rows'] == 2 and exact['n_common'] == 1
    assert exact['n_reversed_common'] == 1
    assert exact['pairs'][0]['reaction'] == 'a'
    assert exact['identity_status_counts'] == {'exact_family_equation_match':1,
                                               'bridge_family_identity_or_extent_unresolved':1}
    assert score.evaluate(rows, labels, predictions)['n_common'] == 2
    with pytest.raises(ValueError, match='paper-row mismatch'):
        score.evaluate([dict(rows[0], paper_row='12')], labels, predictions, identity=identity)
    with pytest.raises(ValueError, match='bridge.*mismatch'):
        score.evaluate([dict(rows[0], bridge_delta_I0_minus_I025='3.1')], labels, predictions, identity=identity)
    with pytest.raises(ValueError, match='missing identity'):
        score.evaluate(rows, labels, predictions, identity={'a':identity['a']})
    with pytest.raises(ValueError, match='unknown identity'):
        score.evaluate([rows[0]], labels, predictions, identity={'a':{'status':'unknown'}})
    with pytest.raises(ValueError, match='missing audited scale/bridge'):
        score.evaluate([rows[0]], labels, predictions,
                       identity={'a':{'status':'exact_family_equation_match','author_paper_row':10}})
    with pytest.raises(ValueError, match='invalid MetaG prediction/reference'):
        score.evaluate([rows[0]], {'a':{'exp':float('nan')}}, predictions, identity=identity)


def test_ec_centered_correlation_handles_class_offsets_and_zero_variance():
    rows = [{'reaction':'a','MetaG':1,'reference':1}, {'reaction':'b','MetaG':2,'reference':2},
            {'reaction':'c','MetaG':100,'reference':-100}, {'reaction':'d','MetaG':101,'reference':-99}]
    reactions = {'a':{'EC':'2.1.1.1'},'b':{'EC':'2.1.1.1'},
                 'c':{'EC':'3.1.1.1'},'d':{'EC':'3.1.1.1'}}
    assert score.within_class_correlation(rows, reactions, 'MetaG') == pytest.approx(1.)
    assert score.within_class_correlation(rows[:1], reactions, 'MetaG') is None


def test_identity_document_rejects_incomplete_or_inconsistent_metadata(tmp_path):
    import json
    path=tmp_path/'identity.json'
    rows={'a':{'status':'exact_family_equation_match','author_paper_row':1,
               'signed_author_scale':1.,'bridge_I0_noMg_minus_native_kj':0.}}
    data={'n_requested':1,'status_counts':{'exact_family_equation_match':1,
                                           'bridge_family_identity_or_extent_unresolved':0},'rows':rows}
    path.write_text(json.dumps(data))
    assert score.load_identity(path)==rows
    data['status_counts']['exact_family_equation_match']=0
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError, match='status totals'):score.load_identity(path)
    data['status_counts']['exact_family_equation_match']=1
    del rows['a']['bridge_I0_noMg_minus_native_kj']
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError, match='incomplete scale/bridge'):score.load_identity(path)
    rows['a']['bridge_I0_noMg_minus_native_kj']=0.
    rows['b']={'status':'exact_family_equation_match','author_paper_row':1}
    data['n_requested']=2;data['status_counts']['exact_family_equation_match']=2
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError, match='duplicate author paper row'):score.load_identity(path)


def test_cli_defaults_to_corrected_predictions_and_records_their_hash(tmp_path):
    import csv
    import hashlib
    import json

    identity_path = BENCHMARK / 'dgbyg_family_identity_complete.json'
    identity = json.loads(identity_path.read_text())['rows']
    labels = json.loads((BENCHMARK / 'reference_mg_bounded.json').read_text())['2.0']['labels']
    source = BENCHMARK / 'metag_corrected_source_predictions.json'
    predictions = json.loads(source.read_text())['predictions']
    rid = next(k for k, record in identity.items()
               if record['status'] == 'exact_family_equation_match'
               and 'signed_author_scale' in record and k in labels
               and predictions.get(k, {}).get('dG') is not None)
    record = identity[rid]
    pairs = tmp_path / 'pairs.tsv'
    with pairs.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, delimiter='\t', fieldnames=[
            'reaction', 'paper_row', 'scale', 'bridge_delta_I0_minus_I025', 'dgbyg'])
        writer.writeheader()
        writer.writerow({'reaction': rid, 'paper_row': record['author_paper_row'],
                         'scale': record['signed_author_scale'],
                         'bridge_delta_I0_minus_I025': record['bridge_I0_noMg_minus_native_kj'],
                         'dgbyg': 0})
    output = tmp_path / 'scored.json'
    score.main(['--pairs', str(pairs), '--output', str(output)])
    result = json.loads(output.read_text())
    assert result['input_sha256']['predictions'] == hashlib.sha256(source.read_bytes()).hexdigest()
    assert result['input_sha256']['identity'] == hashlib.sha256(identity_path.read_bytes()).hexdigest()
    assert result['n_common'] == 1
    assert result['pairs'][0]['MetaG'] == predictions[rid]['dG']


def test_complete_identity_preserves_archived_rows_and_supplies_all_bridges():
    import hashlib
    import json
    import math

    old_path = BENCHMARK / 'dgbyg_family_identity.json'
    old = json.loads(old_path.read_text())
    complete = json.loads((BENCHMARK / 'dgbyg_family_identity_complete.json').read_text())
    assert complete['previous_identity_sha256'] == hashlib.sha256(old_path.read_bytes()).hexdigest()
    assert complete['n_requested'] == old['n_requested'] == 164
    assert complete['status_counts'] == old['status_counts']
    assert set(complete['rows']) == set(old['rows'])
    filled = 0
    for rid, before in old['rows'].items():
        after = complete['rows'][rid]
        assert after['status'] == before['status']
        assert after['author_paper_row'] == before['author_paper_row']
        assert math.isfinite(after['signed_author_scale']) and after['signed_author_scale'] != 0
        assert math.isfinite(after['bridge_I0_noMg_minus_native_kj'])
        if 'signed_author_scale' in before:
            assert after['signed_author_scale'] == before['signed_author_scale']
            assert after['bridge_I0_noMg_minus_native_kj'] == before['bridge_I0_noMg_minus_native_kj']
        else:
            filled += 1
    assert filled == 16
    assert score.load_identity(BENCHMARK / 'dgbyg_family_identity_complete.json') == complete['rows']
