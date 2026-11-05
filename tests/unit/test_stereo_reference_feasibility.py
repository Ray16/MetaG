"""Reference feasibility must not confuse training overlap with independence."""
import csv
import json
from pathlib import Path
import pytest
from benchmark.stereo_reference_feasibility import audit, family_holdout, stereo_candidates, structural_rows
from benchmark.benchmark_feasibility import ALLOWED_COLUMNS

ROOT=Path(__file__).resolve().parents[2]


def row(smiles, kind='formation dg'):
    return dict(description='fixture',reaction=smiles,type=kind,T='298.15',I='0.25',pH='7',pMg='14',EC='')


def test_complete_projected_formation_review_has_no_independent_released_model_targets():
    result=audit(ROOT/'benchmark/feasibility_validation/training671_structural.csv',
                 ROOT/'benchmark/feasibility_validation/historical_inputs.json')
    counts=result['counts']
    assert counts['formation_rows']==222
    assert counts['stereo_distinct_pairs']==17 and counts['nonstereo_graph_groups']==11
    assert counts['pairs_with_unassigned_atom_stereocentres']==9
    assert counts['pair_previews_below_0_5']==1
    assert counts['pairs_with_any_pool_historical_structure_witness']==13
    assert counts['released_training_independent_pairs']==0
    assert result['experimental_references_admitted']==0 and not result['launch_allowed']
    # An absent witness is explicitly not treated as proof of historical independence.
    absent=[p for p in result['pairs'] if not any(p[k]['historical_structure_witnesses']['broad'] for k in ['left','right'])]
    assert len(absent)==1 and not absent[0]['released_model_training_independent']


def test_family_holdout_removes_both_formation_labels_reaction_labels_and_unknown_rows():
    a,b='F/C=C/F','F/C=C\\F'
    rows=[row(a),row(b),row(a+' = '+b,'reaction dg'),row('*'),row('CCO')]
    plan=family_holdout(rows,[a,b])
    assert plan['retained_rows']==[4]
    assert {r['row'] for r in plan['excluded_rows']}=={0,1,2,3}
    assert all('held-out molecular family' in r['reasons'] for r in plan['excluded_rows'][:3])
    assert not plan['launch_allowed'] and not plan['training_performed']


def test_exact_enantiomers_are_filtered_but_incomplete_encodings_are_flagged():
    pairs,_=stereo_candidates([row('C[C@H](O)CC'),row('C[C@@H](O)CC')])
    assert pairs==[]
    pairs,_=stereo_candidates([row('C[C@H](O)CC'),row('CC(O)CC')])
    assert len(pairs)==1 and pairs[0]['right']['unspecified_atom_stereocentres']
    ez,_=stereo_candidates([row('F/C=C/F'),row('F/C=C\\F')])
    assert len(ez)==1


def test_numerical_label_columns_are_not_accepted_in_the_structural_input(tmp_path):
    path=tmp_path/'training.csv'
    with path.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=[*ALLOWED_COLUMNS,'dG'])
        writer.writeheader();writer.writerow(dict(row('CCO'),dG='-123'))
    with pytest.raises(ValueError,match='outcome-free'):
        structural_rows(path)
