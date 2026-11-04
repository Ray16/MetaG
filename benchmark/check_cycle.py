"""Check specified reaction cycles against original structures and predictions.

The structural check is exact up to numerical tolerance for canonical isomeric
SMILES, charge, stoichiometry, explicit proton reservoir and declared sugar pools.
It does not assert that routed species, pKa approximations or sampled energies
are consistent between reactions.
"""
import argparse
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

from rdkit import Chem

from score import load_predictions

HERE = Path(__file__).resolve().parent


def check_cycle(cycle, reactions, predictions):
    stoich = defaultdict(float)
    declarations = defaultdict(set)
    values = []
    for rid, weight in cycle.items():
        if not isinstance(weight, (int, float)) or not math.isfinite(weight) or weight == 0:
            raise ValueError(f'{rid}: finite nonzero weight required')
        rx = reactions[rid]
        if rx.get('pka_sites') or rx.get('species_protonation_models'):
            raise ValueError(f'{rid}: caller-owned protonation requires a separate cycle audit')
        prediction = predictions.get(rid)
        if prediction is None or not math.isfinite(prediction):
            raise ValueError(f'{rid}: missing or nonfinite prediction')
        for name, (coefficient, charge, smiles) in rx['species'].items():
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError(f'{rid}/{name}: invalid SMILES')
            key = (int(charge), Chem.MolToSmiles(mol, isomericSmiles=True))
            stoich[key] += weight * coefficient
            pool = rx.get('species_state_pools', {}).get(name)
            declarations[key].add(json.dumps(pool, sort_keys=True))
        stoich[('proton_reservoir', 'H+')] += weight * rx.get('n_Hplus', 0)
        values.append({'reaction': rid, 'weight': weight, 'dG': prediction})
    residual = [{'charge': q, 'smiles': smi, 'coefficient': amount}
                for (q, smi), amount in sorted(stoich.items(), key=lambda x: str(x[0]))
                if abs(amount) > 1e-9]
    conflicts = [{'charge': q, 'smiles': smi, 'assignments': sorted(assignments)}
                 for (q, smi), assignments in declarations.items() if len(assignments) > 1]
    return {'structurally_closed': not residual, 'state_pool_consistent': not conflicts,
            'residual_stoichiometry': residual, 'state_pool_conflicts': conflicts,
            'energy_sum_kj_mol': math.fsum(v['weight'] * v['dG'] for v in values),
            'rounding_bound_kj_mol': 0.05 * math.fsum(abs(v['weight']) for v in values),
            'predictions': values}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cycles', type=Path, required=True)
    parser.add_argument('--reactions', type=Path, default=HERE.parent / 'src/metag/data/reactions_opentecr_std.json')
    parser.add_argument('--predictions', type=Path, default=HERE / 'metag_0.3.0_predictions.json')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(argv)
    cycles = json.loads(args.cycles.read_text())
    reactions, predictions = json.loads(args.reactions.read_text()), load_predictions(args.predictions)
    reports = {k: check_cycle(v, reactions, predictions) for k, v in cycles.items()}
    output = {'input_sha256': {name: hashlib.sha256(getattr(args, name).read_bytes()).hexdigest()
                               for name in ('cycles', 'reactions', 'predictions')}, 'cycles': reports,
              'scope': 'Rounded archived reaction predictions; structural closure does not establish error attribution.'}
    args.output.write_text(json.dumps(output, indent=2, allow_nan=False) + '\n')
    print(json.dumps(output, indent=2))
    if not all(r['structurally_closed'] and r['state_pool_consistent'] for r in reports.values()):
        raise SystemExit('Cycle inputs or state-pool definitions are inconsistent')


if __name__ == '__main__':
    main()
