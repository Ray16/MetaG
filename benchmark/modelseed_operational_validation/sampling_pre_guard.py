"""Outcome-free stratified sampling and finite-population operational coverage intervals."""
from collections import Counter, defaultdict
from fractions import Fraction
import hashlib
import random
import re


def ec_stratum(value):
    values=str(value or '').split('|')
    classes=set()
    for value in values:
        match=re.fullmatch(r'\s*([1-7])(?:\.[\d-]+){0,3}\s*',value)
        if match:classes.add(match[1])
    return next(iter(classes)) if len(classes)==1 else ('multiple' if classes else 'unannotated')


def select(frame,per_stratum,seed):
    """Source IDs are the units; repeated chemical equations remain source records."""
    if type(per_stratum) is not int or per_stratum<1:raise ValueError('Positive sample size required')
    if not re.fullmatch(r'[0-9a-f]{64}',seed):raise ValueError('Recorded256bit hexadecimal seed required')
    groups=defaultdict(list);seen=set()
    for row in frame:
        if set(row)!={'reaction','stratum'}:raise ValueError('Requires outcome-free source-ID/stratum projection')
        rid=row['reaction'];stratum=row['stratum']
        if not isinstance(rid,str) or not rid or rid in seen:raise ValueError('Duplicate or invalid source ID')
        if not isinstance(stratum,str) or not stratum:raise ValueError('Invalid stratum')
        groups[stratum].append(rid);seen.add(rid)
    rng=random.Random(int(seed,16));selected=[];strata={}
    for stratum,members in sorted(groups.items()):
        members=sorted(members);n=min(per_stratum,len(members));chosen=rng.sample(members,n)
        strata[stratum]=dict(population=len(members),sample_size=n,inclusion_probability=str(Fraction(n,len(members))))
        selected.extend(dict(reaction=rid,stratum=stratum,inclusion_probability=str(Fraction(n,len(members)))) for rid in chosen)
    random.Random(int(hashlib.sha256((seed+'|queue').encode()).hexdigest(),16)).shuffle(selected)
    return dict(population=len(seen),sample_size=len(selected),strata=strata,selected=selected,
        seed=seed,selection_scope='Simple random source-ID samples without replacement within EC strata; no energy/label selection.')


def hypergeometric_interval(N,n,k,alpha=.05):
    """Invert exact equal-tailed hypergeometric tests for population success count."""
    from scipy.stats import hypergeom
    if any(type(x) is not int for x in [N,n,k]) or not 0<=k<=n<=N or N<1 or n<1:
        raise ValueError('Require integer0<=k<=n<=N and nonempty sample/population')
    if not 0<alpha<1:raise ValueError('Alpha must lie in (0,1)')
    low,high=k,N-n+k
    if k==0:lower=0
    else:
        while low<high:
            mid=(low+high)//2
            if hypergeom.sf(k-1,N,mid,n)>=alpha/2:high=mid
            else:low=mid+1
        lower=low
    low,high=k,N-n+k
    if k==n:upper=N
    else:
        while low<high:
            mid=(low+high+1)//2
            if hypergeom.cdf(k,N,mid,n)>=alpha/2:low=mid
            else:high=mid-1
        upper=low
    return dict(lower_success_count=lower,upper_success_count=upper,
        lower=lower/N,upper=upper/N,confidence=1-alpha)


def coverage(design,outcomes,source_total,all_required_primitives_terminal):
    """Withhold population inference if the common budget censors any raw opportunity."""
    expected={r['reaction']:r for r in design['selected']}
    if set(outcomes)!=set(expected) or any(type(v) is not bool for v in outcomes.values()):
        raise ValueError('One boolean numerical-availability outcome required per selected source ID')
    if type(source_total) is not int or source_total<design['population']:
        raise ValueError('Full source denominator must include the eligible frame')
    if not all_required_primitives_terminal:
        return dict(status='inconclusive_common_budget_censoring',sample_requested=len(expected),
            sample_estimated=sum(outcomes.values()),population_numeric_interval=None,
            reason='Not every required raw input received its declared attempt; no population numerical-availability estimate.')
    H=len(design['strata'])
    if not H:raise ValueError('No strata')
    point=0.;lower=upper=0;rows={}
    for stratum,record in design['strata'].items():
        members=[r['reaction'] for r in design['selected'] if r['stratum']==stratum]
        N,n=record['population'],record['sample_size'];assert len(members)==n
        k=sum(outcomes[rid] for rid in members)
        interval=hypergeometric_interval(N,n,k,.05/H)
        point+=N*k/n;lower+=interval['lower_success_count'];upper+=interval['upper_success_count']
        rows[stratum]=dict(population=N,sampled=n,estimated=k,interval=interval)
    frame=design['population']
    return dict(status='complete_probability_sample',strata=rows,sample_requested=len(expected),
        sample_estimated=sum(outcomes.values()),eligible_population=frame,source_population=source_total,
        eligible_frame_point=point/frame,source_population_point=point/source_total,
        population_numeric_interval=dict(eligible_frame=[lower/frame,upper/frame],full_source=[lower/source_total,upper/source_total],
            confidence_at_least=.95,method='Bonferroni simultaneous exact finite-population hypergeometric intervals, aggregated with known stratum weights.'),
        scope='Numerical availability conditional on the fixed input/state/compute policy; not experimental accuracy, total uncertainty or a genome-scale walltime extrapolation.')
