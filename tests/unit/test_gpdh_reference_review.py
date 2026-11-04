"""Reference-state normalization must not turn a held candidate into a Gibbs label."""
import importlib.util
import json
import math
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[2]
spec=importlib.util.spec_from_file_location('gpdh_reference',ROOT/'benchmark/gpdh_exchange_reference.py')
module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
PACKET=ROOT/'benchmark/gpdh_exchange_reference_review'


def test_total_substrate_normalization_and_direction():
    facts=json.loads((PACKET/'source_facts.json').read_text())
    rates=facts['rates_M_inv_s_inv']
    result=module.nominal_exchange(rates,.06,.45)
    # Independent expression uses restored total-concentration rate slopes.
    ratio=(rates['forward_G3P']/(rates['reverse_DHAP_carbonyl']*.45)) / (
        rates['forward_EtG']/(rates['reverse_GA_carbonyl']*.06))
    assert result['total_substrate_exchange_ratio']==pytest.approx(ratio)
    assert result['total_minus_RTln_ratio_kj']==pytest.approx(-.00831446261815324*298.15*math.log(ratio))
    assert result['normalization_change_kj']>4.9


def test_source_scenarios_do_not_admit_equilibrium_labels():
    result=module.audit(PACKET)
    assert result['primary_references_admitted']==result['secondary_thermodynamic_labels_admitted']==0
    assert not result['predictions_examined']
    assert not result['primary_sources_verified']
    assert result['scenario_difference_kj']==pytest.approx(.00831446261815324*298.15*math.log(.55/.45))


def test_source_candidate_chemical_identity_and_condition_request():
    from rdkit import Chem
    from metag.chem import is_balanced
    from metag.reactions import validate_reaction
    from metag.validation.calculation_conditions import calculation_condition_issues
    request=json.loads((PACKET/'candidate_request.json').read_text())
    validate_reaction(request)
    assert is_balanced(request['species'],request['n_Hplus'])
    g3p=Chem.MolFromSmiles(request['species']['R_glycerol_3_phosphate'][2])
    assert [cip for _,cip in Chem.FindMolChiralCenters(g3p)]==['R']
    assert len(calculation_condition_issues(request))==2
    assert 'exp' not in request
