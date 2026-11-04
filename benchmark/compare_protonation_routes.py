"""Compare current routing with an experimental unconditional neutral-base route.

CPU only, no species energies or labels are read. This is an impact audit, not
a recommendation to use the unrestricted alternative for production chemistry.
"""
import argparse
import hashlib
import json
from pathlib import Path
from unittest.mock import patch

from metag.routing import pka_transform
from metag.routing.planner import route_reaction

HERE = Path(__file__).resolve().parent


def route_all(reactions):
    return {rid: route_reaction(rx, allow_truncate=False, log=lambda *_: None) for rid, rx in reactions.items()}


def compare(reactions):
    baseline = route_all(reactions)
    original = pka_transform.build_ph0_reaction

    def forced(*args, **kwargs):
        kwargs['force_base'] = True
        return original(*args, **kwargs)

    with patch.object(pka_transform, 'build_ph0_reaction', forced):
        alternative = route_all(reactions)
    result = {}
    for rid in reactions:
        a, ar, _ = baseline[rid]
        b, br, _ = alternative[rid]
        fields = ('species', 'pka_sites', 'n_Hplus')
        changed = any(a.get(field) != b.get(field) for field in fields)
        if changed or ar.get('errors') or br.get('errors'):
            result[rid] = {'changed': changed,
                           'baseline_errors': ar.get('errors', []), 'alternative_errors': br.get('errors', []),
                           'baseline_ph0': ar.get('ph0'), 'alternative_ph0': br.get('ph0'),
                           'baseline_species': a.get('species'), 'alternative_species': b.get('species'),
                           'baseline_pka_sites': a.get('pka_sites'), 'alternative_pka_sites': b.get('pka_sites'),
                           'baseline_n_Hplus': a.get('n_Hplus'), 'alternative_n_Hplus': b.get('n_Hplus')}
    return {'n_reactions': len(reactions), 'n_changed': sum(v['changed'] for v in result.values()),
            'n_baseline_errors': sum(bool(v['baseline_errors']) for v in result.values()),
            'n_alternative_errors': sum(bool(v['alternative_errors']) for v in result.values()),
            'new_error_ids': [rid for rid, v in result.items() if v['alternative_errors'] and not v['baseline_errors']],
            'affected': result}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reactions', type=Path, default=HERE.parent / 'src/metag/data/reactions_opentecr_std.json')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(argv)
    result = compare(json.loads(args.reactions.read_text()))
    result['input_sha256'] = hashlib.sha256(args.reactions.read_bytes()).hexdigest()
    result['scope'] = __doc__
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'affected'}, indent=2))


if __name__ == '__main__':
    main()
