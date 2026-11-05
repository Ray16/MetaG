"""Audit/replay the experimental oxidized-anchor approximation; no quantum runtime."""
from decimal import Decimal, localcontext
from pathlib import Path
import argparse
import copy
import hashlib
import importlib.util
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
import zipfile

read = lambda p: json.loads(Path(p).read_text())
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()


def load_helper(packet):
    spec = importlib.util.spec_from_file_location('anchored_dependency_audit', Path(packet) / 'global_amino_reference.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def audit(packet, dependency):
    packet, dependency = Path(packet), Path(dependency)
    helper = load_helper(packet)
    helper.verify_manifest(packet)
    frozen = read(packet / 'freeze.json')
    for name, digest in frozen['local_files_sha256'].items():
        assert sha(packet / name) == digest, name
    for name, digest in frozen['dependency_files_sha256'].items():
        assert sha(dependency / name) == digest, name
    raw_audit = helper.audit(dependency)
    catalog = read(dependency / 'catalog596.json')
    raw = {(r['charge'], r['smiles']): r for r in catalog['records']}
    with zipfile.ZipFile(packet / 'outcomes.zip') as archive:
        full = json.loads(archive.read('full.json'))
        candidate = json.loads(archive.read('anchored.json'))
    inputs = read(dependency / 'inputs.json')
    labels = read(dependency / 'reference.json')['2.0']['labels']
    report = read(packet / 'report.json')
    invariances = read(packet / 'invariances.json')
    assert len(full['results']) == len(candidate['results']) == len(inputs) == 364
    with zipfile.ZipFile(dependency / 'outcomes.zip') as archive:
        original = json.loads(archive.read('reassembly596/full_extended596.json'))
    for rid, row in full['results'].items():
        old = original['results'][rid]
        assert row['status'] == old['status']
        if row['status'] == 'estimated':
            assert row['result']['dG_unrounded'] == old['result']['dG_unrounded']
    common = sorted(rid for rid in labels if all(doc['results'][rid]['status'] == 'estimated'
                     for doc in [full, candidate]))
    assert common == report['common_ids']
    for arm, doc in [('full', full), ('anchored', candidate)]:
        errors = [doc['results'][rid]['result']['dG_unrounded'] - labels[rid]['exp'] for rid in common]
        assert abs(math.fsum(map(abs, errors)) / len(errors) - report['paired_common'][arm]['MAE']) < 1e-10
        assert sum(abs(e) > 20 for e in errors) == report['paired_common'][arm]['n_abs_error_gt20']
    maximum = 0.
    mapped_ids = set()
    shared = {}
    with localcontext() as context:
        context.prec = 50
        D = lambda v: Decimal(str(v))
        standard = D('.008314') * D('298.15') * (D('.082057') * D('298.15')).ln()
        proton = D('-26.3') - D('1104.5') - D('2.303') * D('.008314') * D('298.15') * 7
        for rid, row in candidate['results'].items():
            if row['status'] != 'estimated':
                continue
            result, rx = row['result'], inputs[rid]
            terms = [D(rx.get('n_Hplus', 0)) * proton]
            sensitivities = {}
            for name, record in result['state_records'].items():
                assert shared.setdefault(record['state_id'], record['record_sha256']) == record['record_sha256']
                terms.append(D(rx['species'][name][0]) * D(record['potential_kj']))
                for identity, primitive in record['sampling_inputs'].items():
                    original_raw = raw[(primitive['charge'], primitive['smiles'])]
                    assert primitive['sigma_kj'] == original_raw['sigma_kj']
                    total = sensitivities.setdefault(identity, [D(0), D(primitive['sigma_kj'])])
                    total[0] += D(rx['species'][name][0]) * D(primitive['sensitivity'])
                transfer = record.get('molecular_diagnostics', {}).get('anchored_transfer')
                if transfer:
                    mapped_ids.add(record['state_id'])
                    molecular = sum(D(term['coefficient']) * D(raw[(term['charge'], term['smiles'])]['G_kj'])
                                    for term in transfer['raw_terms'])
                    assert abs(molecular + standard - D(record['solution_reference_kj'])) < D('1e-7')
                    assert abs(D(record['solution_reference_kj']) + D(record['protonation_kj'])
                               + D(record['reference_proton_kj']) - D(record['potential_kj'])) < D('1e-7')
            residual = abs(sum(terms) - D(result['dG_unrounded']))
            maximum = max(maximum, float(residual))
            assert residual < D('1e-7'), rid
            sigma = sum((v[0] * v[1]) ** 2 for v in sensitivities.values()).sqrt()
            assert abs(sigma - D(result['U_samp'])) < D('1e-8'), rid
    assert all(abs(r['residual_kj']) < 1e-7 for r in invariances['reaction_checks'])
    for cycle in invariances['available_cycle_basis']:
        residual = math.fsum(n * candidate['results'][rid]['result']['dG_unrounded']
                             for rid, n in cycle['coefficients'].items())
        norm = sum(map(abs, cycle['coefficients'].values()))
        assert abs(residual) / norm < 1e-7
    diagnosis = read(packet / 'applicability_diagnosis.json')
    anchors = {(p['oxidized_sibling_charge'], p['oxidized_sibling_smiles'])
               for p in candidate['ledger']['anchored_reference_plans']}
    for row in diagnosis['rows']:
        result = candidate['results'][row['reaction']]['result']
        cancelled = not any(abs(p['sensitivity']) > 1e-12 and (p['charge'], p['smiles']) in anchors
                            for p in result['sampling_inputs'].values())
        assert cancelled == row['full_anchor_cancels']
    assert not report['production_policy_adopted'] and report['independent_primary_labels_admitted'] == 0
    return dict(status='passed',raw_records_verified=raw_audit['authenticated_raw_records'],
        requests=364, candidate_estimates=sum(r['status']=='estimated' for r in candidate['results'].values()),
        paired_labels=len(common), paired_MAE_kj={a: report['paired_common'][a]['MAE'] for a in ['full','anchored','zero']},
        exact_shared_mapped_states=len(mapped_ids), maximum_decimal_assembly_residual_kj=maximum,
        reaction_invariance_checks=len(invariances['reaction_checks']), cycles=len(invariances['available_cycle_basis']),
        anchor_retaining_experimental_labels=diagnosis['anchor_retaining_label_count'],
        production_policy_adopted=False, independent_primary_labels_admitted=0)


def portable_outcome(document):
    """Verify self-hash and separate only replay-host software inventory.

    The original raw producers, frozen source identity, settings, states,
    diagnostics and all numerical outputs remain part of comparison.
    A data-only replay need not have the original quantum software installed.
    """
    doc = copy.deepcopy(document)
    ledger = doc['ledger']
    digest = ledger.pop('sha256')
    assert hashlib.sha256(json.dumps(ledger, sort_keys=True, allow_nan=False).encode()).hexdigest() == digest
    software = ledger['config'].pop('software')
    assert isinstance(json.loads(software), dict)
    for row in doc['results'].values():
        if row['result'] is not None:
            assert row['result']['config'].pop('software') == software
    # Derived state hashes can change when a floating-point diagnostic differs
    # by one ulp. Verify every state and replace its *references* by the exact
    # molecular state identity; compare the complete payload separately below.
    state_hashes = {}
    for identity, record in ledger['states'].items():
        assert record['state_id'] == identity
        body = dict(record)
        state_digest = body.pop('record_sha256')
        assert hashlib.sha256(json.dumps(body, sort_keys=True, allow_nan=False).encode()).hexdigest() == state_digest
        state_hashes[state_digest] = identity
    def normalize(value):
        if isinstance(value, dict):
            if 'record_sha256' in value and 'plan' in value:
                body = dict(value)
                state_digest = body.pop('record_sha256')
                assert hashlib.sha256(json.dumps(body, sort_keys=True, allow_nan=False).encode()).hexdigest() == state_digest
                assert state_hashes[state_digest] == value['state_id']
            result = {}
            for key, item in value.items():
                if key == 'record_sha256':
                    assert item in state_hashes
                    result[key] = 'verified-state:' + state_hashes[item]
                else:
                    result[key] = normalize(item)
            return result
        if isinstance(value, list):
            return [normalize(item) for item in value]
        return value
    return normalize(doc), json.loads(software)


def first_difference(left, right, path='root'):
    """Bounded diagnostic for a portability failure without dumping the ledger."""
    if type(left) is not type(right):
        return dict(path=path, left_type=type(left).__name__, right_type=type(right).__name__)
    if isinstance(left, dict):
        if left.keys() != right.keys():
            return dict(path=path, missing=sorted(right.keys()-left.keys()), extra=sorted(left.keys()-right.keys()))
        for key in left:
            result = first_difference(left[key], right[key], path + '.' + str(key))
            if result:
                return result
    elif isinstance(left, list):
        if len(left) != len(right):
            return dict(path=path, left_length=len(left), right_length=len(right))
        for i, (a, b) in enumerate(zip(left, right)):
            result = first_difference(a, b, path + '[' + str(i) + ']')
            if result:
                return result
    elif left != right:
        return dict(path=path, left=left, right=right)
    return None


def compare_outcomes(left, right, thresholds):
    """Check every field using the original protocol's absolute tolerances.

    Only derived floating-point quantities use a tolerance. Discrete identity,
    all configuration (including raw-producer inventories), keys and structure
    compare exactly. Ledger/state self-hashes are verified before this call.
    """
    receipt = dict(n_numerical_differences=0, maximum_absolute_difference=0., examples=[])
    def check(a, b, path):
        assert type(a) is type(b), first_difference(a, b, path)
        if isinstance(a, dict):
            assert a.keys() == b.keys(), first_difference(a, b, path)
            if path.endswith('.config'):
                assert a == b, first_difference(a, b, path)
                return
            for key in a:
                check(a[key], b[key], path + '.' + str(key))
        elif isinstance(a, list):
            assert len(a) == len(b), first_difference(a, b, path)
            for i, (aa, bb) in enumerate(zip(a, b)):
                check(aa, bb, path + '[' + str(i) + ']')
        elif isinstance(a, float):
            assert math.isfinite(a) and math.isfinite(b), path
            delta = abs(a - b)
            sampling = any(token in path for token in ['U_samp', 'sigma', 'sampling', 'sensitivity', 'population'])
            tolerance = thresholds['sampling_invariance_kj' if sampling else 'invariance_kj']
            assert delta <= tolerance, dict(path=path, left=a, right=b, difference=delta, tolerance=tolerance)
            if delta:
                receipt['n_numerical_differences'] += 1
                receipt['maximum_absolute_difference'] = max(receipt['maximum_absolute_difference'], delta)
                if len(receipt['examples']) < 10:
                    receipt['examples'].append(dict(path=path, difference=delta, tolerance=tolerance))
        else:
            assert a == b, first_difference(a, b, path)
    check(left, right, 'root')
    return receipt


def replay(packet, dependency):
    packet, dependency = Path(packet).resolve(), Path(dependency).resolve()
    with tempfile.TemporaryDirectory(prefix='metag-anchor-replay-') as name:
        root = Path(name)
        for filename in ['protocol.json', 'freeze.json', 'run_study.py',
                         'anchored_nicotinamide.py', 'global_amino_reference.py']:
            shutil.copyfile(packet / filename, root / filename)
        (root / 'sitecustomize.py').write_text('''import importlib.abc,sys,socket
class NoQuantum(importlib.abc.MetaPathFinder):
 def find_spec(self,fullname,path=None,target=None):
  if fullname.split('.')[0] in {'torch','fairchem','pyscf','gpu4pyscf','cupy'}:
   raise RuntimeError('Quantum runtime forbidden in owned-record replay: '+fullname)
sys.meta_path.insert(0,NoQuantum())
def deny_network(*args,**kwargs): raise RuntimeError('Network forbidden in owned-record replay')
socket.socket.connect=deny_network
''')
        env = dict(os.environ, PYTHONPATH=str(root), PYTHONDONTWRITEBYTECODE='1',
                   OMP_NUM_THREADS='1', OMP_THREAD_LIMIT='1', MKL_NUM_THREADS='1',
                   OPENBLAS_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1')
        result = subprocess.run([sys.executable, str(root / 'run_study.py'), str(dependency)],
                                cwd=root, env=env, capture_output=True, text=True, timeout=180)
        if result.returncode:
            raise RuntimeError(result.stdout + result.stderr)
        thresholds = read(packet / 'protocol.json')['thresholds']
        comparisons = {}
        for filename in ['report.json', 'invariances.json', 'assembly_complete.json']:
            comparisons[filename] = compare_outcomes(read(root / filename), read(packet / filename), thresholds)
        with zipfile.ZipFile(packet / 'outcomes.zip') as archive:
            inventories = {}
            for filename in ['full.json', 'anchored.json']:
                new, replay_software = portable_outcome(read(root / filename))
                old, original_software = portable_outcome(json.loads(archive.read(filename)))
                comparisons[filename] = compare_outcomes(new, old, thresholds)
                inventories[filename] = dict(original=original_software, replay=replay_software)
        return dict(status='passed',requests=364, all_outcomes_and_invariances_reproduced=True,
                    quantum_imports_forbidden=True,network_forbidden=True,
                    software_inventory_comparison=inventories,
                    numerical_comparisons=comparisons,
                    comparison_scope='All scientific fields checked under frozen energy/sampling tolerances; identities/configuration/raw provenance exact. Verified derived hashes resolve to their complete state payloads. Replay-host software inventories reported separately.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--packet', type=Path, default=Path(__file__).with_name('anchored_nicotinamide_validation'))
    parser.add_argument('--dependency-packet', type=Path, default=Path(__file__).with_name('global_amino_reference_validation'))
    parser.add_argument('--cold-replay', action='store_true')
    args = parser.parse_args()
    result = audit(args.packet, args.dependency_packet)
    if args.cold_replay:
        result['cold_replay'] = replay(args.packet, args.dependency_packet)
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
