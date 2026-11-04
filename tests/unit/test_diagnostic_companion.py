"""Analytic, adversarial and invariance controls for portable diagnostics; no QM required."""
import copy,importlib.util,math
from pathlib import Path
import pytest
H=Path(__file__).resolve().parents[2];spec=importlib.util.spec_from_file_location('diagnostic_companion',H/'benchmark/thermochemistry_diagnostics.py');D=importlib.util.module_from_spec(spec);spec.loader.exec_module(D)

def transfer():
 return {'charge':0,'standard_state':'1M_gas_to_1M_solution','rt_kj':2.,'states':[{'id':'a','electronic_kj':0.,'thermal_kj':0.,'solvation_kj':0.},{'id':'b','electronic_kj':2*math.log(3),'thermal_kj':0.,'solvation_kj':2*math.log(2)}],'reference_center_kj':3.}

def network():
 return {'conditions':{'T_K':298.15,'phase':'aqueous'},'energy_units':'kJ/mol','energies_kj':{'ab':1.,'bc':2.,'ac':4.},'stoichiometry':{'ab':{'A':-1,'B':1},'bc':{'B':-1,'C':1},'ac':{'A':-1,'C':1}}}

def test_analytic_two_state_partition_and_reference_gap():
 r=D.transfer_diagnostic(transfer());assert r['transfer_kj']==pytest.approx(-2*math.log(7/8));assert r['minimum_center_gap_under_reweighting_kj']==pytest.approx(3-2*math.log(2))

def test_common_energy_origin_cancels():
 p=transfer();a=D.transfer_diagnostic(p)
 for s in p['states']:s['electronic_kj']-=1.e6
 b=D.transfer_diagnostic(p);assert b['transfer_kj']==pytest.approx(a['transfer_kj'],abs=1e-9)

def test_owned_interval_can_overlap_despite_outside_center():
 p=transfer();p['reference_interval_kj']=[1.,4.];r=D.transfer_diagnostic(p);assert r['center_excluded_on_retained_basis'] and not r['supplied_interval_excluded']

@pytest.mark.parametrize('change',[lambda p:p.update(charge=-1),lambda p:p.update(rt_kj=0.),lambda p:p['states'].append(copy.deepcopy(p['states'][0])),lambda p:p['states'][0].update(solvation_kj=float('nan')),lambda p:p.update(standard_state='unspecified')])
def test_transfer_rejects_unsupported_or_incomplete_definitions(change):
 p=transfer();change(p)
 with pytest.raises(ValueError):D.transfer_diagnostic(p)

def test_three_reaction_radius_has_independent_analytic_oracle():
 r=D.cycle_diagnostic(network());assert r['nullity']==1 and r['radius_kj']==pytest.approx(1/3);assert r['dual_lower_bound_kj']==pytest.approx(1/3);assert r['primal_dual_gap_kj']==pytest.approx(0.,abs=1e-9)

def test_complete_nullspace_detects_hidden_relation():
 p=network();p['stoichiometry']['ac_copy']=dict(p['stoichiometry']['ac']);p['energies_kj']['ac_copy']=8.
 r=D.cycle_diagnostic(p);assert r['nullity']==2;assert r['radius_kj']>=2.-1e-8

def test_coherent_values_and_reaction_reversal():
 p=network();p['energies_kj']['ac']=3.;assert D.cycle_diagnostic(p)['radius_kj']<1e-9
 p=network();p['stoichiometry']['ab']={s:-v for s,v in p['stoichiometry']['ab'].items()};p['energies_kj']['ab']=-1.;assert D.cycle_diagnostic(p)['radius_kj']==pytest.approx(1/3)

def test_cycle_basis_change_does_not_change_radius():
 c=[[1.,1.,-1.,0.],[0.,0.,-1.,1.]];g=[1.,2.,4.,8.]
 a=D.radius_from_cycles(c,g);b=D.radius_from_cycles([[v*3 for v in c[0]],[v+2*w for v,w in zip(c[1],c[0])]],g)
 assert a['radius_kj']==pytest.approx(b['radius_kj'],abs=1e-9)

def test_no_cycles_is_uninformative_not_a_successful_accuracy_test():
 p=network();p['stoichiometry'].pop('ac');p['energies_kj'].pop('ac');r=D.cycle_diagnostic(p);assert r['nullity']==0 and r['diagnosis']=='no_testable_cycles'

@pytest.mark.parametrize('change',[lambda p:p['energies_kj'].update(ac=None),lambda p:p['energies_kj'].update(ac=float('inf')),lambda p:p['energies_kj'].pop('ac'),lambda p:p.update(conditions={}),lambda p:p.update(energy_units='kcal/mol')])
def test_network_missing_values_and_mixed_conventions_are_not_zero_imputed(change):
 p=network();change(p)
 with pytest.raises((ValueError,TypeError)):D.cycle_diagnostic(p)


def test_ill_conditioned_integer_basis_has_the_same_analytic_radius():
    n=36;cycles=[[int(i==j)-int(j==n-1) for j in range(n)] for i in range(n-1)]
    mixed=[list(r) for r in cycles]
    for i in range(1,len(mixed)):mixed[i]=[a+2*b for a,b in zip(cycles[i],cycles[i-1])]
    values=list(range(n));a=D.radius_from_cycles(cycles,values);b=D.radius_from_cycles(mixed,values)
    assert a['radius_kj']==pytest.approx(17.5,abs=1e-9)
    assert b['radius_kj']==pytest.approx(17.5,abs=1e-9)


def test_acyclic_branch_does_not_validate_or_influence_cycle_energy():
    p=network();p['stoichiometry']['de']={'D':-1,'E':1};p['energies_kj']['de']=999.
    r=D.cycle_diagnostic(p);assert r['n_reactions']==4 and r['n_cycle_participating_reactions']==3
    assert 'de' not in r['cycle_participating_reaction_ids']
    assert r['radius_kj']==pytest.approx(1/3)
