import importlib.util
import json
import sys
from pathlib import Path

BENCHMARK = Path(__file__).resolve().parents[2] / 'benchmark'
sys.path.insert(0, str(BENCHMARK))
try:
    spec = importlib.util.spec_from_file_location('cycle_audit', BENCHMARK / 'check_cycle.py')
    audit = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(audit)
finally:
    sys.path.pop(0)


def test_archived_sucrose_cycle_closes_chemically_but_not_in_energy():
    cycle = json.loads((BENCHMARK / 'reaction_cycles.json').read_text())['sucrose_adp_glucose']
    reactions = json.loads((BENCHMARK.parent / 'src/metag/data/reactions_opentecr_std.json').read_text())
    predictions = audit.load_predictions(BENCHMARK / 'metag_0.3.0_predictions.json')
    result = audit.check_cycle(cycle, reactions, predictions)
    assert result['structurally_closed'] and result['state_pool_consistent']
    assert abs(result['energy_sum_kj_mol'] - 4.4) < 1e-9
    assert result['rounding_bound_kj_mol'] < 0.16


def test_cycle_audit_detects_missing_reaction_and_pool_disagreement():
    reaction = {'n_Hplus': 0, 'species': {'a': [-1, 0, 'CCO'], 'b': [1, 0, 'CC=O']}}
    reverse = {'n_Hplus': 0, 'species': {'a': [1, 0, 'CCO'], 'b': [-1, 0, 'CC=O']},
               'species_state_pools': {'a': {'kind': 'reducing_sugar_equilibrium', 'source': 'example'}}}
    incomplete = audit.check_cycle({'one': 1}, {'one': reaction}, {'one': 3})
    assert not incomplete['structurally_closed']
    conflict = audit.check_cycle({'one': 1, 'two': 1}, {'one': reaction, 'two': reverse}, {'one': 3, 'two': -3})
    assert conflict['structurally_closed'] and not conflict['state_pool_consistent']


def test_shared_replay_without_truncation_closes_sucrose_cycle():
    sys.path.insert(0, str(BENCHMARK))
    try:
        spec = importlib.util.spec_from_file_location('cycle_replay', BENCHMARK / 'replay_cycles.py')
        replay = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(replay)
    finally:
        sys.path.pop(0)
    cycles = json.loads((BENCHMARK / 'reaction_cycles.json').read_text())
    reactions = json.loads((BENCHMARK.parent / 'src/metag/data/reactions_opentecr_std.json').read_text())
    result = replay.replay({'sucrose': cycles['sucrose_adp_glucose']}, reactions)['sucrose']['without_truncation']
    assert result['structurally_closed'] and result['state_pool_consistent']
    assert abs(result['artificial_energy_sum_kj_mol']) <= 0.15


def test_forced_neutral_base_replay_closes_the_amino_acid_cycle():
    sys.path.insert(0, str(BENCHMARK))
    try:
        spec = importlib.util.spec_from_file_location('cycle_replay_forced', BENCHMARK / 'replay_cycles.py')
        replay = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(replay)
    finally:
        sys.path.pop(0)
    cycles = json.loads((BENCHMARK / 'reaction_cycles.json').read_text())
    reactions = json.loads((BENCHMARK.parent / 'src/metag/data/reactions_opentecr_std.json').read_text())
    result = replay.replay({'amino': cycles['alanine_glutamate_redox']}, reactions, force_base=True)
    assert abs(result['amino']['without_truncation']['artificial_energy_sum_kj_mol']) <= 0.15


def test_protonation_impact_audit_restores_default_router():
    from metag.routing import pka_transform
    sys.path.insert(0, str(BENCHMARK))
    try:
        spec = importlib.util.spec_from_file_location('compare_routes', BENCHMARK / 'compare_protonation_routes.py')
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    finally:
        sys.path.pop(0)
    reactions = json.loads((BENCHMARK.parent / 'src/metag/data/reactions_opentecr_std.json').read_text())
    original = pka_transform.build_ph0_reaction
    report = module.compare({'rxn00191': reactions['rxn00191']})
    assert pka_transform.build_ph0_reaction is original
    assert report['n_changed'] == 1
