"""Frozen research reassembly; compare labels only after assembly/invariance checks."""
from pathlib import Path
from collections import Counter
from fractions import Fraction
import copy
import hashlib
import json
import math
import os
import sys
import zipfile

HERE = Path(__file__).resolve().parent
read = lambda p: json.loads(Path(p).read_text())
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()


def write(name, value):
    with (HERE / name).open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')


def main():
    freeze = read(HERE / 'freeze.json')
    for name, digest in freeze['local_files_sha256'].items():
        assert sha(HERE / name) == digest, name
    packet = Path(sys.argv[1]).resolve()
    for name, digest in freeze['dependency_files_sha256'].items():
        assert sha(packet / name) == digest, name
    # No raw-source optimization or remote service is invoked by this study.
    from global_amino_reference import audit, extract
    write('raw_record_audit.json', audit(packet))
    source = HERE / 'source'
    source.mkdir(exist_ok=False)
    extract(packet / 'assembly_source.zip', source)
    sys.path[:0] = [str(HERE), str(source / 'assembly_source')]
    from rdkit import Chem
    import sympy as sp
    from metag.settings import Settings
    from metag.state_function import StateFunctionEstimator
    from metag.routing.sugar_states import prepare_declared_sugar_inputs
    from anchored_nicotinamide import AnchoredNicotinamideEstimator

    catalog = read(packet / 'catalog596.json')
    class MissingRecord(RuntimeError):
        pass
    class Backend:
        def __init__(self):
            self.settings = Settings(**read(packet / 'raw_settings.json'))
            self.model_name = self.settings.model
            self.producer = {'cacheable': False, 'kind': 'owned-explicit-profile-reassembly',
                             'catalog_sha256': sha(packet / 'catalog596.json'),
                             'source_profiles': catalog['source_profiles']}
            self.rows = {(r['charge'], r['smiles']): r for r in catalog['records']}
            self.used = {}
        def species(self, charge, smiles, label, warnings, diagnostics):
            key = (charge, Chem.MolToSmiles(Chem.MolFromSmiles(smiles), isomericSmiles=True))
            if key not in self.rows:
                raise MissingRecord(str(charge) + '|' + key[1])
            row = self.rows[key]
            self.used[row['name']] = row['sha256']
            diagnostics.update(source_record=row['name'], source_sha256=row['sha256'],
                               raw_profile=row['raw_profile'], n_minima=row['n_minima'])
            return row['G_kj'], row['sigma_kj']
        def liquid_water(self):
            return catalog['water_liquid_1atm_kj']

    inputs = read(packet / 'inputs.json')
    allowed = {'species', 'n_Hplus', 'species_state_pools', 'species_protonation_models', 'pka_sites'}
    requests = {rid: {k: copy.deepcopy(v) for k, v in rx.items() if k in allowed}
                for rid, rx in inputs.items()}
    with zipfile.ZipFile(packet / 'outcomes.zip') as archive:
        prior = json.loads(archive.read('reassembly596/full_extended596.json'))['results']
    all_results, estimators = {}, {}
    for arm, cls in [('full', StateFunctionEstimator), ('anchored', AnchoredNicotinamideEstimator)]:
        backend = Backend()
        estimator = cls(backend)
        rows = {}
        for rid, request in requests.items():
            try:
                result = estimator.score_reaction(request, key=rid)
                rows[rid] = dict(status='estimated' if result['dG'] is not None else 'workflow_refusal', result=result)
            except MissingRecord as exc:
                rows[rid] = dict(status='missing_archived_molecule', missing=str(exc), result=None)
            if arm == 'full':
                assert rows[rid]['status'] == prior[rid]['status'], rid
                if rows[rid]['status'] == 'estimated':
                    assert abs(result['dG_unrounded'] - prior[rid]['result']['dG_unrounded']) < 1e-7, rid
                    assert abs(result['U_samp'] - prior[rid]['result']['U_samp']) < 1e-8, rid
        all_results[arm] = rows
        estimators[arm] = estimator
        write(arm + '.json', dict(results=rows, ledger=estimator.ledger(), used_raw_records=backend.used))
        print(arm, dict(Counter(r['status'] for r in rows.values())), flush=True)

    estimator = estimators['anchored']
    checks = []
    state_signatures = {}
    for rid, request in requests.items():
        saved = all_results['anchored'][rid]
        if saved['status'] != 'estimated':
            continue
        base = saved['result']
        for record in base['state_records'].values():
            key, digest = record['state_id'], record['record_sha256']
            assert state_signatures.setdefault(key, digest) == digest
        for variant in ['reverse', 'double', 'atom_order', 'species_order']:
            other = copy.deepcopy(request)
            factor = -1 if variant == 'reverse' else 2 if variant == 'double' else 1
            if variant in ['reverse', 'double']:
                for values in other['species'].values():
                    values[0] *= factor
                other['n_Hplus'] = factor * other.get('n_Hplus', 0)
            if variant == 'species_order':
                other['species'] = dict(reversed(list(other['species'].items())))
            if variant == 'atom_order':
                for values in other['species'].values():
                    mol = Chem.MolFromSmiles(values[2])
                    mol = Chem.RenumberAtoms(mol, list(reversed(range(mol.GetNumAtoms()))))
                    values[2] = Chem.MolToSmiles(mol, canonical=False, isomericSmiles=True)
            result = estimator.score_reaction(other, key=rid + '_' + variant)
            assert result['dG'] is not None, (rid, variant, result['suspect'])
            residual = result['dG_unrounded'] - factor * base['dG_unrounded']
            assert abs(residual) < 1e-7, (rid, variant, residual)
            assert abs(result['U_samp'] - abs(factor) * base['U_samp']) < 1e-8
            checks.append(dict(reaction=rid, variant=variant, residual_kj=residual))
    # A separate ledger populated in the opposite reaction order must agree.
    reverse_estimator = AnchoredNicotinamideEstimator(Backend())
    for rid, request in reversed(list(requests.items())):
        if all_results['anchored'][rid]['status'] == 'estimated':
            result = reverse_estimator.score_reaction(request, key=rid)
            assert result['dG_unrounded'] == all_results['anchored'][rid]['result']['dG_unrounded']
            for rec in result['state_records'].values():
                assert rec['record_sha256'] == state_signatures[rec['state_id']]
    ids = sorted(rid for rid, row in all_results['anchored'].items() if row['status'] == 'estimated')
    columns, keys = [], set()
    for rid in ids:
        rx, column = requests[rid], {}
        for name, (nu, q, smi) in rx['species'].items():
            if name in rx.get('species_state_pools', {}):
                _, pools = prepare_declared_sugar_inputs({'species': {'state': [1, q, smi]},
                    'species_state_pools': {'state': rx['species_state_pools'][name]}})
                key = json.dumps(['pool', q, sorted(s['smiles'] for s in pools['state']['source_states'])])
            else:
                key = json.dumps(['molecule', q, Chem.MolToSmiles(Chem.MolFromSmiles(smi), isomericSmiles=True)])
            column[key] = column.get(key, Fraction(0)) + Fraction(str(nu))
            keys.add(key)
        if rx.get('n_Hplus', 0):
            column['explicit_proton_reservoir'] = Fraction(str(rx['n_Hplus']))
            keys.add('explicit_proton_reservoir')
        columns.append(column)
    matrix = sp.Matrix([[sp.Rational(col.get(k, 0)) for col in columns] for k in sorted(keys)])
    cycles = []
    for vector in matrix.nullspace():
        denominator = math.lcm(*(int(v.q) for v in vector))
        integers = [int(v * denominator) for v in vector]
        divisor = math.gcd(*(abs(v) for v in integers))
        integers = [v // divisor for v in integers]
        residual = math.fsum(n * all_results['anchored'][rid]['result']['dG_unrounded'] for rid, n in zip(ids, integers))
        norm = sum(map(abs, integers))
        assert abs(residual) / norm < 1e-7
        cycles.append(dict(coefficients={rid: n for rid, n in zip(ids, integers) if n},
                           residual_kj=residual, normalized_residual_kj=residual / norm))
    write('invariances.json', dict(reaction_checks=checks, shared_states=len(state_signatures),
        opposite_population_order_checked=True, available_cycle_basis=cycles,
        scope='Literal molecular/pool stoichiometry and explicit proton reservoir; not accuracy evidence.'))
    write('assembly_complete.json', dict(assembly_complete=True, labels_read=False,
        n_requested=len(requests), n_reaction_checks=len(checks), n_cycles=len(cycles)))
    # Only now load the development-exposed experimental label table.
    labels = read(packet / 'reference.json')['2.0']['labels']
    values = {arm: {rid: row['result']['dG_unrounded'] for rid, row in rows.items()
                    if row['status'] == 'estimated'} for arm, rows in all_results.items()}
    def metrics(errors):
        if not errors:
            return {'n': 0}
        return dict(n=len(errors), MAE=math.fsum(map(abs, errors)) / len(errors),
            RMSE=math.sqrt(math.fsum(e * e for e in errors) / len(errors)),
            bias=math.fsum(errors) / len(errors), n_abs_error_gt20=sum(abs(e) > 20 for e in errors))
    def paired(targets):
        return {arm: metrics([values[arm][rid] - labels[rid]['exp'] for rid in targets])
                for arm in values} | {'zero': metrics([-labels[rid]['exp'] for rid in targets])}
    common = sorted(set(values['full']) & set(values['anchored']) & set(labels))
    affected = [rid for rid in common if any(r['molecular_diagnostics'].get('anchored_transfer')
                for r in all_results['anchored'][rid]['result']['state_records'].values()
                if 'molecular_diagnostics' in r)]
    unaffected = sorted(set(common) - set(affected))
    assert all(abs(values['full'][rid] - values['anchored'][rid]) < 1e-7 for rid in unaffected)
    report = dict(scope='Development-exposed cached reassembly; not independent validation or production adoption.',
        n_requested=len(requests), freeze_sha256=sha(HERE / 'freeze.json'),
        arms={arm: dict(status_counts=dict(Counter(row['status'] for row in all_results[arm].values())),
              conditional_error=metrics([v - labels[rid]['exp'] for rid, v in vals.items() if rid in labels]),
              zero_same_targets=metrics([-labels[rid]['exp'] for rid in vals if rid in labels]))
              for arm, vals in values.items()},
        paired_common=paired(common), paired_affected=paired(affected), paired_unaffected=paired(unaffected),
        paired_by_EC={ec: paired([rid for rid in common if str(inputs[rid]['EC']).split('.')[0] == ec])
                      for ec in sorted({str(inputs[rid]['EC']).split('.')[0] for rid in common})},
        common_ids=common, affected_ids=affected,
        gained_requests=sorted(set(values['anchored']) - set(values['full'])),
        lost_requests=sorted(set(values['full']) - set(values['anchored'])),
        n_cycles=len(cycles), n_reaction_checks=len(checks),
        maximum_cycle_residual_kj=max([abs(c['normalized_residual_kj']) for c in cycles] or [0]),
        production_policy_adopted=False, independent_primary_labels_admitted=0,
        limitation='The full oxidized anchor cancels within a family but is unvalidated in scaffold synthesis/exchange; omitted spectator interactions are not zero by construction.')
    write('report.json', report)
    print(json.dumps({k: v for k, v in report.items() if k not in ['common_ids', 'affected_ids', 'gained_requests', 'lost_requests']}, indent=2), flush=True)


if __name__ == '__main__':
    main()
