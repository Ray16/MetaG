"""Evaluate the narrow experimental transamination route on CPU, with no experimental labels."""
import argparse
import hashlib
import json
from pathlib import Path
from unittest.mock import patch

from metag.routing import pka_transform
from compare_protonation_routes import route_all
from experimental_transamination import wrap
from replay_cycles import replay

HERE = Path(__file__).resolve().parent


def evaluate(reactions, cycles):
    baseline = route_all(reactions)
    original = pka_transform.build_ph0_reaction
    with patch.object(pka_transform, 'build_ph0_reaction', wrap(original)):
        alternative = route_all(reactions)
        cycle_results = replay(cycles, reactions)
    changed, new_errors = [], []
    for rid in reactions:
        before, bdiag, _ = baseline[rid]
        after, adiag, _ = alternative[rid]
        if any(before.get(k) != after.get(k) for k in ('species', 'pka_sites', 'n_Hplus')):
            changed.append(rid)
        if adiag.get('errors') and not bdiag.get('errors'):
            new_errors.append(rid)
    return {'n_reactions': len(reactions), 'changed_ids': changed, 'new_error_ids': new_errors,
            'replay_cycles': cycle_results}


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reactions', type=Path, default=HERE.parent / 'src/metag/data/reactions_opentecr_std.json')
    p.add_argument('--cycles', type=Path, default=HERE / 'reaction_cycles.json')
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args(argv)
    result = evaluate(json.loads(args.reactions.read_text()), json.loads(args.cycles.read_text()))
    result['scope'] = __doc__
    result['input_sha256'] = {k: hashlib.sha256(getattr(args, k).read_bytes()).hexdigest()
                              for k in ('reactions', 'cycles')}
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps({'n_reactions': result['n_reactions'], 'changed_ids': result['changed_ids'],
                      'new_error_ids': result['new_error_ids'],
                      'cycle_sums': {k: v['without_truncation']['artificial_energy_sum_kj_mol']
                                     for k, v in result['replay_cycles'].items()}}, indent=2))


if __name__ == '__main__':
    main()
