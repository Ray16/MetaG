"""An independent audit must catch coherent offsets and resealed bad records."""
import copy
import importlib.util
import json
import math
from pathlib import Path

import pytest

from metag.constants import CONDITIONS,G_HPLUS,STD_STATE_KJ
from metag.settings import Settings
from metag.state_function import StateFunctionEstimator

PATH = Path(__file__).parents[2]/'benchmark/verify_state_function_ledger.py'
spec = importlib.util.spec_from_file_location('verify_state_function_ledger',PATH)
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


class Backend:
    model_name = 'ledger-controls'
    producer = {'cacheable':False,'kind':'synthetic'}

    def __init__(self):
        self.settings = Settings()

    def species(self,q,smiles,label,warnings,diagnostics):
        return -1000.+sum(map(ord,smiles))%31,.2

    def liquid_water(self):
        return -900.


def payload(panel=None):
    forward = {'species':{'a':[-1,0,'CCO'],'b':[1,0,'COC']}}
    reverse = {'species':{'a':[1,0,'CCO'],'b':[-1,0,'COC']}}
    panel = panel or {'forward':forward,'reverse':reverse}
    estimator = StateFunctionEstimator(Backend())
    results = {}
    for rid,rx in panel.items():
        row = estimator.score_reaction(rx)
        assert row['dG'] is not None,row['suspect']
        results[rid] = {'status':'assembled','dG':row['dG'],'dG_unrounded':row['dG_unrounded'],
            'sampling_sigma_kj':row['U_samp'],
            'assignments':{n:r['state_id'] for n,r in row['state_records'].items()}}
    return {'schema':'metag-owned-state-function-ledger-v1','arm':'full','conditions':dict(CONDITIONS),
            'proton_energy_kj':G_HPLUS,'std_state_kj':STD_STATE_KJ,'KT_kj':2.4789,
            'inputs':panel,'results':results,'ledger':estimator.ledger()}


def reseal(data,key):
    row = data['ledger']['states'][key]
    row.pop('record_sha256')
    row['record_sha256'] = audit.digest(row)
    data['ledger'].pop('sha256')
    data['ledger']['sha256'] = audit.digest(data['ledger'])


def test_independent_reassembly_and_closed_cycle_pass():
    result = audit.verify(payload(),cycles={'independent_cycles':[{'forward':1,'reverse':1}]})
    assert result['n_reactions_checked']==2 and result['n_evaluable_cycles']==1
    assert result['max_cycle_residual_kj']==0.
    assert result['n_raw_records_authenticated']==0


def test_opposite_reaction_offsets_close_the_cycle_but_fail_molecular_reassembly():
    data = payload()
    for rid,offset in [('forward',1.),('reverse',-1.)]:
        data['results'][rid]['dG_unrounded'] += offset
        data['results'][rid]['dG'] = round(data['results'][rid]['dG_unrounded'],1)
    assert data['results']['forward']['dG_unrounded']+data['results']['reverse']['dG_unrounded']==0.
    with pytest.raises(ValueError,match='reaction assembly'):
        audit.verify(data)


def test_resealed_state_offset_still_fails_component_arithmetic():
    data = payload()
    key = next(iter(data['ledger']['states']))
    data['ledger']['states'][key]['potential_kj'] += 1.
    reseal(data,key)
    with pytest.raises(ValueError,match='molecular state potential'):
        audit.verify(data)


def test_resealed_sampling_spread_must_match_partition_derivatives():
    data = payload()
    key = next(iter(data['ledger']['states']))
    row = data['ledger']['states'][key]
    sample = next(iter(row['sampling_inputs'].values()))
    sample['sensitivity'] += .2
    row['sampling_sigma_kj'] = abs(sample['sensitivity']*sample['sigma_kj'])
    reseal(data,key)
    with pytest.raises(ValueError,match='sampling partition derivative'):
        audit.verify(data)


def test_incomplete_reaction_assignment_is_rejected():
    data = payload()
    data['results']['forward']['assignments'].pop('a')
    with pytest.raises(ValueError,match='missing molecular assignments'):
        audit.verify(data)


def test_unbalanced_reaction_cannot_pass_with_a_matching_energy_sum():
    data = payload()
    data['inputs']['forward']['species']['a'][0] = -2
    row = data['results']['forward']
    state = data['ledger']['states'][row['assignments']['a']]
    row['dG_unrounded'] -= state['potential_kj']
    row['dG'] = round(row['dG_unrounded'],1)
    with pytest.raises(ValueError,match='element/charge balance'):
        audit.verify(data)


def test_absence_remains_an_unavailable_result():
    data = payload()
    data['results']['forward'] = {'status':'archived_molecular_record_unavailable','dG':None}
    assert audit.verify(data)['n_reactions_checked']==1


def test_unlisted_intervention_and_changed_conditions_are_rejected():
    data = payload()
    with pytest.raises(ValueError,match='absent from declared protocol'):
        audit.verify(data,protocol={'arms':['another']})
    data['conditions']['pH'] = 8.
    with pytest.raises(ValueError,match='Unsupported conditions'):
        audit.verify(data)


def pool_payload():
    glc = 'OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O'
    fru = 'OC[C@H]1OC(O)(CO)[C@@H](O)[C@@H]1O'
    reaction = {'species':{'a':[-1,0,glc],'b':[1,0,fru]},
                'species_state_pools':{n:{'kind':'reducing_sugar_equilibrium','source':'Synthetic pool'}
                                       for n in ('a','b')}}
    return payload({'isomerization':reaction})


def test_missing_pool_constituent_is_rejected_even_after_resealing():
    data = pool_payload()
    key = next(k for k,r in data['ledger']['states'].items() if r['plan'].get('kind'))
    data['ledger']['states'][key]['states'].pop()
    reseal(data,key)
    with pytest.raises(ValueError,match='Incomplete declared pool state basis'):
        audit.verify(data)


def test_declared_lowest_state_control_retains_all_states_but_only_selected_derivatives():
    data = pool_payload()
    data['arm'] = 'lowest_sugar_state'
    key = next(k for k,r in data['ledger']['states'].items() if r['plan'].get('kind'))
    row = data['ledger']['states'][key]
    winner = min(row['states'],key=lambda r:(r['potential_kj'],r['smiles']))
    child = data['ledger']['states'][winner['state_id']]
    row.update(potential_kj=child['potential_kj'],sampling_inputs=copy.deepcopy(child['sampling_inputs']),
               sampling_sigma_kj=child['sampling_sigma_kj'],pool_intervention='lowest_sugar_state',
               selected_constituent=winner['state_id'])
    for member in row['states']:
        member['conditional_population'] = float(member['state_id']==winner['state_id'])
    reseal(data,key)
    # Isolate state-partition verification from a separate reaction result.
    data['results']['isomerization'] = {'status':'not_scored_in_this_control','dG':None}
    assert audit.verify(data)['n_pools_checked']==2
