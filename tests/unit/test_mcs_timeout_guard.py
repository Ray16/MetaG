"""An incomplete MCS must never supply an atom map to the production router."""
from types import SimpleNamespace
from unittest.mock import patch
import pytest
from rdkit import Chem
from metag.routing import truncate, planner

def test_cancelled_search_never_yields_atom_map():
    incomplete=SimpleNamespace(canceled=True,numAtoms=1,smartsString='[#6]')
    with patch.object(truncate,'_mcs',return_value=incomplete):
        with pytest.raises(RuntimeError,match='partial atom map rejected'):
            truncate.mcs_atom_map(Chem.MolFromSmiles('CCO'),Chem.MolFromSmiles('COC'))

@pytest.fixture
def fresh_pair_cache():
    truncate._pair_maps.cache_clear()
    yield
    truncate._pair_maps.cache_clear()


def test_cancelled_search_fails_closed_through_planner(fresh_pair_cache):
    incomplete=SimpleNamespace(canceled=True,numAtoms=1,smartsString='[#6]')
    rx={'species':{'a':[-1,0,'CCO'],'b':[1,0,'COC']},'n_Hplus':0,'note':''}
    routes={'errors':[],'warnings':[],'truncated':False,'trunc_radius':None,'reverted':[]}
    with patch('metag.routing.truncation_gate.prefer_full',return_value=False), \
         patch.object(truncate,'truncation_radius',return_value=2), \
         patch.object(truncate,'_mcs',return_value=incomplete):
        result=planner._route_truncation(rx,routes,True,lambda *_:None)
    assert result is rx                                  # full molecules: no partial atom map used
    assert not routes['truncated']
    assert not routes['errors']                          # an exhausted search is not a refusal
    assert any('partial atom map rejected' in w for w in routes['warnings'])


def test_budget_is_deterministic_and_raises_budget_error():
    a = Chem.MolFromSmiles('Nc1nc(=O)c2c([nH]1)NC[C@@H]1CN(c3ccc(C(=O)NC)cc3)CN21')
    b = Chem.MolFromSmiles('Nc1nc(=O)c2c([nH]1)NC[C@@H]1C[N+](c3ccc(C(=O)NC)cc3)=CN21')
    for _ in range(2):                                   # same outcome on every call, independent of load
        with pytest.raises(truncate.MCSBudgetExceeded, match='deterministic budget'):
            truncate._mcs(a, b, budget=50)
    assert truncate._mcs(Chem.MolFromSmiles('CCO'), Chem.MolFromSmiles('CCN')).numAtoms == 2


def test_budget_exhaustion_routes_full_molecules(fresh_pair_cache):
    rx = {'species': {'a': [-1, 0, 'CCO'], 'b': [1, 0, 'COC']}, 'n_Hplus': 0, 'note': ''}
    routes = {'errors': [], 'warnings': [], 'truncated': False, 'trunc_radius': None, 'reverted': []}
    with patch('metag.routing.truncation_gate.prefer_full',
               side_effect=truncate.MCSBudgetExceeded('budget')):
        result = planner._route_truncation(rx, routes, True, lambda *_: None)
    assert result is rx and routes['prefer_full'] and not routes['errors']
