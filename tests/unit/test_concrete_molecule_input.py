"""Undefined R groups must fail before model or molecular work is attempted."""
import json
from pathlib import Path

import pytest
from rdkit import Chem

from metag.pipeline import score_reaction
from metag.reactions import validate_reaction
from metag.settings import Settings
from metag.state_function import StateFunctionEstimator


class Backend:
    model_name = 'input-boundary-test'
    producer = {'cacheable':False,'kind':'synthetic'}

    def __init__(self):
        self.settings = Settings()
        self.calls = []

    def species(self,*args):
        self.calls.append(args)
        raise AssertionError('Undefined molecular input reached the backend')


@pytest.mark.parametrize('smiles',['*','CC*','CC(=O)O*'])
def test_wildcards_are_rejected_before_legacy_or_state_function_molecular_work(smiles):
    rx = {'species':{'left':[-1,0,smiles],'right':[1,0,smiles]}}
    backend = Backend()
    with pytest.raises(ValueError,match='wildcard/dummy'):
        score_reaction(rx,backend,log=lambda *_:None)
    with pytest.raises(ValueError,match='wildcard/dummy'):
        StateFunctionEstimator(backend).score_reaction(rx)
    assert not backend.calls


def test_guard_does_not_change_the_evaluated_364_input_domain():
    path = Path(__file__).parents[2]/'src/metag/data/reactions_opentecr_std.json'
    panel = json.loads(path.read_text())
    assert len(panel)==364
    for reaction in panel.values():
        validate_reaction(reaction)
        assert all(atom.GetAtomicNum()>0 for _,_,smiles in reaction['species'].values()
                   for atom in Chem.MolFromSmiles(smiles).GetAtoms())
