"""Mapped input must never silently switch the water standard-state convention."""
import copy
import pytest
from metag.pipeline import score_reaction
from metag.reactions import ReactionInput


class NoMolecularWork:
    pu = None

    def __getattr__(self, name):
        raise AssertionError('Input should be rejected before backend access: ' + name)


@pytest.mark.parametrize('smiles', ['[OH2:1]', '[OH2:0]', '[H:2]O[H:3]',
                                    '[H:2][O:1][H:3]', '[CH3:1][OH:2]'])
def test_atom_maps_rejected_before_backend_or_routing(smiles):
    reaction = {'note': 'mapped input control', 'n_Hplus': 0,
                'species': {'mapped': [-1, 0, smiles], 'copy': [1, 0, smiles]}}
    original = copy.deepcopy(reaction)
    with pytest.raises(ValueError, match='atom-mapped input SMILES are unsupported'):
        score_reaction(reaction, NoMolecularWork(), log=lambda *_: None)
    with pytest.raises(ValueError, match='atom-mapped input SMILES are unsupported'):
        ReactionInput.from_mapping(reaction)
    assert reaction == original
