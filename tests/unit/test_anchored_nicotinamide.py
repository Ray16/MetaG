import math
from pathlib import Path

import pytest
from rdkit import Chem

from benchmark.anchored_nicotinamide import (
    AnchoredNicotinamideEstimator, OX_CORE, RED_CORE, anchor_plan,
)
from metag.chem import is_balanced
from metag.energetics.provenance import content_hash
from metag.settings import Settings
from metag.state_function import _combine_sampling


RED = 'NC(=O)C1=CN([C@@H]2O[C@H](CO)[C@@H](O)[C@H]2O)C=CC1'
OX = 'NC(=O)c1ccc[n+]([C@@H]2O[C@H](CO)[C@@H](O)[C@H]2O)c1'
canonical = lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(s), isomericSmiles=True)


def test_exact_spectator_stereochemistry_and_balanced_transfer():
    plan = anchor_plan(0, RED)
    assert plan['oxidized_sibling_smiles'] == canonical(OX)
    terms = {'target': [-1, 0, RED]}
    terms.update({str(i): [r['coefficient'], r['charge'], r['smiles']]
                  for i, r in enumerate(plan['raw_terms'])})
    assert is_balanced(terms, 0)
    epimer = RED.replace('[C@@H](O)', '[C@H](O)', 1)
    assert anchor_plan(0, epimer)['oxidized_sibling_smiles'] != plan['oxidized_sibling_smiles']
    assert anchor_plan(0, RED.replace('CO)', 'COP(=O)(O)O)'))['oxidized_sibling_smiles'] != plan['oxidized_sibling_smiles']
    mol = Chem.MolFromSmiles(RED)
    reordered = Chem.RenumberAtoms(mol, list(reversed(range(mol.GetNumAtoms()))))
    assert anchor_plan(0, Chem.MolToSmiles(reordered, canonical=False)) == plan


def test_scope_and_charge_are_explicit():
    assert anchor_plan(0, RED_CORE) is None
    assert anchor_plan(1, OX_CORE) is None
    assert anchor_plan(1, OX) is None
    assert anchor_plan(0, 'CCO') is None
    with pytest.raises(ValueError, match='charge'):
        anchor_plan(1, RED)


class Backend:
    settings = Settings()
    model_name = 'test-constant-primitives'
    producer = {'cacheable': False, 'kind': 'test'}
    def species(self, q, smi, label, warnings, diagnostics):
        return {canonical(OX): (100., 7.), RED_CORE: (20., 2.), OX_CORE: (10., 3.),
                canonical(RED): (900., 99.)}[smi]
    def liquid_water(self):
        return 0.


def test_shared_anchor_sampling_cancels_in_redox_difference():
    estimator = AnchoredNicotinamideEstimator(Backend())
    reduced = estimator.state(0, RED)
    oxidized = estimator.state(1, OX)
    assert reduced['solution_reference_kj'] - oxidized['solution_reference_kj'] == pytest.approx(10.)
    assert reduced['sampling_sigma_kj'] == pytest.approx(math.sqrt(62.))
    sigma, inputs = _combine_sampling([(1., reduced), (-1., oxidized)])
    assert sigma == pytest.approx(math.sqrt(13.))
    assert inputs[content_hash({'charge': 1, 'smiles': canonical(OX)})]['sensitivity'] == 0.
    assert canonical(RED) not in [r['smiles'] for r in inputs.values()]


def test_same_state_survives_other_context_and_returned_record_mutation():
    estimator = AnchoredNicotinamideEstimator(Backend())
    before = estimator.state(0, RED)
    returned = estimator.state(0, RED)
    returned['potential_kj'] = -999999.
    estimator.state(1, OX)
    assert estimator.state(0, RED) == before
    reverse_order = AnchoredNicotinamideEstimator(Backend())
    reverse_order.state(1, OX)
    assert reverse_order.state(0, RED) == before


def test_owned_record_panel_and_shared_sampling_reconstruct_independently():
    from benchmark.audit_anchored_nicotinamide import audit
    root = Path(__file__).resolve().parents[2] / 'benchmark'
    result = audit(root / 'anchored_nicotinamide_validation', root / 'global_amino_reference_validation')
    assert result['candidate_estimates'] == 160
    assert result['paired_labels'] == 112
    assert result['paired_MAE_kj']['anchored'] == pytest.approx(7.9362908034874895)
    assert result['reaction_invariance_checks'] == 640 and result['cycles'] == 18
    assert result['anchor_retaining_experimental_labels'] == 0
    assert not result['production_policy_adopted']


def test_frozen_both_arm_replay_without_quantum_runtime_or_network():
    from benchmark.audit_anchored_nicotinamide import replay
    root = Path(__file__).resolve().parents[2] / 'benchmark'
    result = replay(root / 'anchored_nicotinamide_validation', root / 'global_amino_reference_validation')
    assert result['requests'] == 364
    assert result['all_outcomes_and_invariances_reproduced']
