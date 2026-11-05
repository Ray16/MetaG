"""Operational coverage requires source weights, exact sampling units and no budget censoring."""
from math import comb
import pytest
from benchmark.modelseed_sampling import coverage,ec_stratum,hypergeometric_interval,select

SEED='13'*32


def test_source_strata_and_order_independent_sampling_with_census_of_small_strata():
    assert ec_stratum('1.1.1.1|2.7.1.2')=='multiple'
    assert ec_stratum('null')=='unannotated' and ec_stratum('7.1.1.1')=='7'
    frame=[dict(reaction=f'a{i}',stratum='1') for i in range(20)]+[dict(reaction='small',stratum='6')]
    a=select(frame,4,SEED);b=select(list(reversed(frame)),4,SEED)
    assert a==b and a['sample_size']==5
    assert a['strata']['1']['inclusion_probability']=='1/5'
    assert a['strata']['6']['inclusion_probability']=='1'
    assert len({r['reaction'] for r in a['selected']})==5


def test_sampling_rejects_duplicate_ids_and_outcome_columns():
    with pytest.raises(ValueError,match='Duplicate'):select([dict(reaction='a',stratum='1')]*2,4,SEED)
    with pytest.raises(ValueError,match='outcome-free'):select([dict(reaction='a',stratum='1',dG=-3)],4,SEED)


def test_exact_intervals_have_at_least_nominal_design_coverage_without_iid_assumptions():
    N,n=9,4
    for successes in range(N+1):
        covered=0.
        for observed in range(n+1):
            if not 0<=observed<=successes or not 0<=n-observed<=N-successes:continue
            probability=comb(successes,observed)*comb(N-successes,n-observed)/comb(N,n)
            interval=hypergeometric_interval(N,n,observed)
            if interval['lower_success_count']<=successes<=interval['upper_success_count']:covered+=probability
        assert covered>=.95-1e-12
    for k in range(10):
        interval=hypergeometric_interval(9,9,k)
        assert interval['lower_success_count']==interval['upper_success_count']==k


def test_stratum_weights_prevent_equal_sample_size_bias_and_retain_full_denominator():
    frame=[dict(reaction=f'large{i}',stratum='1') for i in range(900)]+[dict(reaction=f'small{i}',stratum='6') for i in range(100)]
    design=select(frame,4,SEED)
    outcomes={r['reaction']:r['stratum']=='1' for r in design['selected']}
    result=coverage(design,outcomes,10000,True)
    assert result['sample_estimated']/result['sample_requested']==.5
    assert result['eligible_frame_point']==pytest.approx(.9)
    assert result['source_population_point']==pytest.approx(.09)
    assert result['population_numeric_interval']['full_source'][1]<=.1


def test_incomplete_raw_opportunities_withhold_population_inference_and_never_zero_fill():
    design=select([dict(reaction='r',stratum='1')],4,SEED)
    result=coverage(design,{'r':False},10,False)
    assert result['status']=='inconclusive_common_budget_censoring'
    assert result['population_numeric_interval'] is None and 'source_population_point' not in result
    with pytest.raises(ValueError,match='boolean numerical-availability'):coverage(design,{},10,True)
    with pytest.raises(ValueError,match='Explicit boolean'):coverage(design,{'r':False},10,'false')
