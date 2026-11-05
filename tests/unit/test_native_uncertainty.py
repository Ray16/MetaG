"""Retain shared-member covariance and do not promote development evidence."""
import importlib.util
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[2]
spec=importlib.util.spec_from_file_location('native_uncertainty_audit',ROOT/'benchmark/native_uncertainty.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def test_common_formation_baseline_cancels_before_uncertainty():
    molecules={'a':{'members_kj':[1000.,2000.,3000.]},
               'b':{'members_kj':[1007.,2007.,3007.]}}
    values=module.reaction_members(molecules,[{'molecule':'a','coefficient':-1},
                                             {'molecule':'b','coefficient':1}])
    assert values==[7.,7.,7.]
    assert module.population_summary(values)==(7.,0.)
    assert module.population_summary(molecules['a']['members_kj'])[1]>800


def test_reversal_and_extent_preserve_uncertainty_units():
    molecules={'a':{'members_kj':[2.,4.,9.]},'b':{'members_kj':[5.,8.,10.]}}
    forward=module.reaction_members(molecules,[{'molecule':'a','coefficient':-1},
                                              {'molecule':'b','coefficient':1}])
    reversed_twice=module.reaction_members(molecules,[{'molecule':'a','coefficient':2},
                                                     {'molecule':'b','coefficient':-2}])
    mean,sd=module.population_summary(forward)
    opposite,scaled_sd=module.population_summary(reversed_twice)
    assert opposite==pytest.approx(-2*mean)
    assert scaled_sd==pytest.approx(2*sd)


def test_molecular_member_dimensions_must_match():
    molecules={'a':{'members_kj':[1.,2.]},'b':{'members_kj':[1.,2.,3.]}}
    with pytest.raises(ValueError,match='Member dimensions differ'):
        module.reaction_members(molecules,[{'molecule':'a','coefficient':-1},
                                          {'molecule':'b','coefficient':1}])


@pytest.mark.parametrize('values',[[],[float('nan')],[True]])
def test_undefined_member_statistics_are_refused(values):
    with pytest.raises(ValueError):module.population_summary(values)


def test_complete_released_member_packet_stays_developmental():
    report=module.audit(ROOT/'benchmark/native_uncertainty_validation')
    assert report['reaction_member_sums_checked']==6900
    assert report['risk_entries_independently_checked']==2484
    assert report['native_full_cohort_MAE_kj']==pytest.approx(2.979300287040526)
    assert report['comparator_AURC']['disagreement']<report['comparator_AURC']['spread']
    assert not report['independent_reliability_advantage_demonstrated']
