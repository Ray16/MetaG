"""Diagnose routing consistency with deterministic artificial species energies.

Values are deliberately nonphysical; differences locate route dependence and
must never be interpreted as prediction errors or experimental free energies.
"""
import argparse
import hashlib
import json
from pathlib import Path
from unittest.mock import patch

from rdkit import Chem
from metag.pipeline import score_reaction
from metag.settings import Settings

from check_cycle import check_cycle

HERE = Path(__file__).resolve().parent


class ReplayBackend:
    model_name = 'cycle-replay-nonphysical'
    producer = {'cacheable': False, 'kind': 'deterministic_replay'}
    settings = Settings()

    def species(self, q, smiles, label, warnings, diagnostics):
        canonical = Chem.MolToSmiles(Chem.MolFromSmiles(smiles), isomericSmiles=True)
        key = hashlib.sha256(f'{q}|{canonical}'.encode()).hexdigest()
        return -1000.0 + (int(key[:8], 16) % 100000) / 100, 0.0

    def liquid_water(self):
        return -1000.0


def replay(cycles, reactions, *, force_base=False):
    backend = ReplayBackend()
    out = {}
    from metag.routing import pka_transform
    original = pka_transform.build_ph0_reaction

    def forced(*args, **kwargs):
        kwargs['force_base'] = True
        return original(*args, **kwargs)

    # Diagnostic override only. Unrestricted base routing can be wrong for other chemistry.
    with patch.object(pka_transform, 'build_ph0_reaction', forced if force_base else original):
        for name, cycle in cycles.items():
            records = {}
            for allow_truncate in (True, False):
                predictions, routes = {}, {}
                for rid in cycle:
                    result = score_reaction(reactions[rid], backend, key=rid, allow_truncate=allow_truncate,
                                            log=lambda *_: None)
                    predictions[rid] = result.get('dG_raw')
                    routes[rid] = {'truncated': result['routes'].get('truncated'),
                                   'cofactor_ring': result['routes'].get('cofactor_ring'),
                                   'ph0': result['routes'].get('ph0'),
                                   'species_scored': result.get('species_scored'),
                                   'suspect': result.get('suspect')}
                report = check_cycle(cycle, reactions, predictions)
                records['default_truncation' if allow_truncate else 'without_truncation'] = {
                    'structurally_closed': report['structurally_closed'],
                    'state_pool_consistent': report['state_pool_consistent'],
                    'artificial_energy_sum_kj_mol': report['energy_sum_kj_mol'],
                    'routes': routes}
            out[name] = records
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cycles', type=Path, default=HERE / 'reaction_cycles.json')
    p.add_argument('--reactions', type=Path, default=HERE.parent / 'src/metag/data/reactions_opentecr_std.json')
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--force-base', action='store_true', help='experimental override; not a validated production route')
    args = p.parse_args(argv)
    cycles, reactions = json.loads(args.cycles.read_text()), json.loads(args.reactions.read_text())
    output = {'scope': __doc__, 'force_base': args.force_base,
              'input_sha256': {k: hashlib.sha256(getattr(args, k).read_bytes()).hexdigest()
                for k in ('cycles', 'reactions')}, 'cycles': replay(cycles, reactions, force_base=args.force_base)}
    args.output.write_text(json.dumps(output, indent=2, allow_nan=False) + '\n')
    print(json.dumps({name: {route: data['artificial_energy_sum_kj_mol'] for route, data in result.items()}
                      for name, result in output['cycles'].items()}, indent=2))


if __name__ == '__main__':
    main()
