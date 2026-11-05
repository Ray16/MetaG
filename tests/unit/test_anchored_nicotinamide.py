import math
import copy
import json
import zipfile
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


def test_portable_comparison_separates_only_host_software_and_verifies_ledger():
    from benchmark.audit_anchored_nicotinamide import portable_outcome
    root = Path(__file__).resolve().parents[2] / 'benchmark/anchored_nicotinamide_validation'
    with zipfile.ZipFile(root / 'outcomes.zip') as archive:
        original = json.loads(archive.read('full.json'))
    expected, software = portable_outcome(original)
    changed = copy.deepcopy(original)
    inventory = json.dumps(dict(software, torch='cpu-test-host', xtb=None), sort_keys=True)
    changed['ledger']['config']['software'] = inventory
    for row in changed['results'].values():
        if row['result'] is not None:
            row['result']['config']['software'] = inventory
    changed['ledger'].pop('sha256')
    changed['ledger']['sha256'] = content_hash(changed['ledger'])
    actual, recorded = portable_outcome(changed)
    assert actual == expected and recorded['torch'] == 'cpu-test-host'
    altered = copy.deepcopy(changed)
    altered['ledger']['raw_molecules'][0]['record'][0] += 1.
    with pytest.raises(AssertionError):
        portable_outcome(altered)
    # A new valid ledger hash cannot hide a changed scientific value either.
    altered['ledger'].pop('sha256')
    altered['ledger']['sha256'] = content_hash(altered['ledger'])
    assert portable_outcome(altered)[0] != expected


def test_one_ulp_sampling_drift_uses_frozen_tolerance_without_hiding_changes():
    from benchmark.audit_anchored_nicotinamide import compare_outcomes
    limits = dict(invariance_kj=1e-7, sampling_invariance_kj=1e-8)
    original = {'result': {'U_samp': 0.7804640940598064, 'dG': -10.}}
    changed = copy.deepcopy(original)
    changed['result']['U_samp'] = 0.7804640940598063
    receipt = compare_outcomes(changed, original, limits)
    assert receipt['n_numerical_differences'] == 1
    assert receipt['maximum_absolute_difference'] < 2e-16
    changed['result']['U_samp'] += 1e-6
    with pytest.raises(AssertionError):
        compare_outcomes(changed, original, limits)
    with pytest.raises(AssertionError):
        compare_outcomes({'config': {'fixed_model': 1. + 1e-10}},
                         {'config': {'fixed_model': 1.}}, limits)
