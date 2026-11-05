"""Curation must conserve source defects and retain explicit chemical identity."""
from fractions import Fraction
from pathlib import Path
import copy
import json

import pytest

from benchmark.modelseed_cofactor_curation import (
    column_balance, rebase_equation, review_conversion, validate_registry,
)
from benchmark.audit_modelseed_cofactor_curation import audit, cold_targets

PACKET=Path(__file__).resolve().parents[2]/'benchmark/modelseed_cofactor_curation_validation'
read=lambda name:json.loads((PACKET/name).read_text())


def test_owned_primary_structures_preserve_eight_carbon_centers_and_proton_deltas():
    registry=read('registry.json');source=read('original_compound_rows.json')
    validate_registry(registry,source)
    for cid,row in registry.items():
        target=row['target'];original=row['source']
        result=review_conversion(int(original['charge']),original['smiles'],target['charge'],target['smiles'])
        assert result['declared_carbon_stereocenters_preserved']==8
        assert result['delta_H']==(1 if cid in {'cpd00004','cpd00005'} else 0)
        assert not result['energy_equivalence_established']


def test_source_drift_and_unauthorized_cross_reference_fail():
    registry=read('registry.json');original=read('original_compound_rows.json')
    changed=copy.deepcopy(original);changed['cpd00003']['charge']='0'
    with pytest.raises(ValueError,match='source annotation changed'):
        validate_registry(registry,changed)
    altered=copy.deepcopy(registry);altered['cpd00003']['chebi_id']='CHEBI:57945'
    with pytest.raises(ValueError,match='cross-reference'):
        validate_registry(altered,original)


@pytest.mark.parametrize('source,target', [
    ('C[C@H](O)C(=O)O','C[C@@H](O)C(=O)O'),
    ('C/C=C/C','C/C=C\\C'),
    ('CCO','COC'),
    ('CCO','[13CH3]CO'),
])
def test_different_stereo_connectivity_or_isotope_is_not_a_proton_basis_change(source,target):
    with pytest.raises(ValueError):review_conversion(0,source,0,target)


def test_a_redox_change_cannot_be_relabelled_as_protonation():
    ox=read('registry.json')['cpd00003']['target'];red=read('registry.json')['cpd00004']['target']
    with pytest.raises(ValueError,match='ordinary protons'):
        review_conversion(ox['charge'],ox['smiles'],red['charge'],red['smiles'])


def test_product_proton_rebasing_preserves_balanced_and_unbalanced_columns():
    # Ester hydrolysis: methyl acetate+water -> acetate+methanol+H+.
    original={'ester':dict(smiles='COC(C)=O',charge='0'),
              'water':dict(smiles='O',charge='0'),
              'acetate':dict(smiles='CC(=O)[O-]',charge='-1'),
              'methanol':dict(smiles='CO',charge='0')}
    target=copy.deepcopy(original);target['acetate']=dict(smiles='CC(=O)O',charge='0')
    terms={('ester','0'):Fraction(-1),('water','0'):Fraction(-1),
           ('acetate','0'):Fraction(1),('methanol','0'):Fraction(1),('cpd00067','0'):Fraction(1)}
    converted,changes=rebase_equation(terms,{'acetate':{'delta_H':1}})
    assert changes=={'0':'-1'} and ('cpd00067','0') not in converted
    assert column_balance(terms,original)==column_balance(converted,target)==({},'0')
    broken=dict(terms);broken[('methanol','0')]=Fraction(2)
    bad,_=rebase_equation(broken,{'acetate':{'delta_H':1}})
    assert column_balance(broken,original)==column_balance(bad,target)
    assert column_balance(bad,target)!=({},'0')


def test_fractional_reversal_scaling_and_compartment_local_protons():
    terms={('red','a'):Fraction(-1,2),('red','b'):Fraction(1,2)}
    new,changes=rebase_equation(terms,{'red':{'delta_H':1}})
    assert changes=={'a':'1/2','b':'-1/2'}
    assert new[('cpd00067','a')]==Fraction(1,2) and new[('cpd00067','b')]==Fraction(-1,2)
    reverse={key:-2*nu for key,nu in terms.items()}
    other,_=rebase_equation(reverse,{'red':{'delta_H':1}})
    assert other=={key:-2*nu for key,nu in new.items()}


def test_independent_offline_reconstruction_retains_all_inputs_and_boundaries():
    result=audit(PACKET)
    assert result['n_requested']==43650 and result['affected_inputs']==7511
    assert result['newly_admitted_inputs']==3147 and result['lost_inputs']==0
    assert result['known_residuals_preserved']==5854 and result['curated_input_passes']==9719
    assert result['anchor_retaining_inputs']==12 and result['physical_predictions_added']==0


def test_frozen_target_plans_pass_original_gates_without_quantum_runtime_or_network():
    result=cold_targets(PACKET)
    assert result['targets']==4 and result['original_stereo_gate_used']
    assert result['quantum_imports_forbidden'] and result['network_forbidden']
