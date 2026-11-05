"""Authenticate the secondary fresh-water reference, without inference or network."""
from decimal import Decimal, localcontext
import copy
import hashlib
import json
from pathlib import Path
import zipfile

from benchmark.operational_records import authenticate


def main():
    project = Path(__file__).resolve().parents[1]
    packet = Path(__file__).with_name('modelseed_operational_water_reference_validation')
    read = lambda p: json.loads(p.read_text())
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    manifest = read(packet / 'manifest.json')
    assert sha(project / 'benchmark/operational_records.py') == manifest['raw_auditor_source_sha256']
    assert {str(p.relative_to(packet)) for p in packet.rglob('*') if p.is_file()} == set(manifest['files']) | {'manifest.json'}
    for name, row in manifest['files'].items():
        path = packet / name
        assert path.stat().st_size == row['bytes'] and sha(path) == row['sha256'], name
    source = project / manifest['precompute_source']['path']
    assert sha(source) == manifest['precompute_source']['sha256']
    with zipfile.ZipFile(source) as archive:
        prefix = 'water_repair_v100/'
        freeze_bytes = archive.read(prefix + 'freeze.json')
        frozen = json.loads(freeze_bytes)
        for name, digest in frozen['files_sha256'].items():
            assert hashlib.sha256(archive.read(prefix + name)).hexdigest() == digest, name
        helper = archive.read(prefix + 'owned_source/benchmark/operational_water.py')
        assert hashlib.sha256(helper).hexdigest() == manifest['water_helper_source_sha256']
        protocol_sha256 = hashlib.sha256(archive.read(prefix + 'protocol.json')).hexdigest()
    phase = packet / 'phase'
    startup, parent = read(phase / 'startup.json'), read(packet / 'parent_H100_startup.json')
    runtime = copy.deepcopy(startup['runtime'])
    assert startup['status'] == 'passed' and max(map(abs, startup['point_energy_differences_kj'])) <= .001
    assert runtime['hardware_capability'] == [7, 0] and runtime['nitro_eligible_here'] is False
    hashed = lambda x: hashlib.sha256(json.dumps(x, sort_keys=True, allow_nan=False).encode()).hexdigest()
    assert hashed(runtime) == startup['producer']['model']['inference_sha256']
    runtime['hardware_capability'], runtime['nitro_eligible_here'] = [9, 0], True
    assert hashed(runtime) == parent['producer']['model']['inference_sha256']
    normalized = copy.deepcopy(startup['producer'])
    normalized['model']['inference_sha256'] = parent['producer']['model']['inference_sha256']
    assert normalized == parent['producer']
    task, = read(phase / 'tasks.json')
    result = read(phase / 'results' / (task['slot'] + '.json'))
    assert result['freeze_sha256'] == hashlib.sha256(freeze_bytes).hexdigest()
    assert result['protocol_sha256'] == protocol_sha256
    raw, catalog = authenticate(task, result, phase, {}, startup['producer']['model'])
    assert raw['settings'] == read(phase / 'raw_settings.json')
    assert result['producer'] == startup['producer']
    with localcontext() as context:
        context.prec = 50
        total = sum(Decimal(str(raw[k])) for k in [
            'electronic_kj', 'experimental_water_hydration_kj', 'thermal_kj', 'concentration_kj'])
        error = abs(total - Decimal(str(raw['G'])))
        assert error < Decimal('1e-7')
    proof, launch = read(phase / 'terminal_process_audit.json'), read(phase / 'launch_manifest.json')
    assert launch['freeze_sha256'] == hashlib.sha256(freeze_bytes).hexdigest()
    assert proof['terminal'] and proof['run_id'] == launch['run_id']
    assert proof['launch_manifest_sha256'] == sha(phase / 'launch_manifest.json')
    assert len(proof['ownership_checks']) == 3 and not any(r['same_live_process'] for r in proof['ownership_checks'])
    assert not proof['remaining_marked_predictors']
    assert read(packet / 'original_water_import_failure.json')['status'] == 'execution_failure'
    assert read(packet / 'H100_water_resource_interruption.json')['status'] == 'external_resource_interruption'
    print(json.dumps(dict(status='passed', records_authenticated=1,
                          decimal_component_error_kj=str(error), n_imag=raw['n_imag'],
                          symmetrynumber=raw['symmetrynumber'], original_failures_preserved=True,
                          experimental_accuracy_tested=False), indent=2))


if __name__ == '__main__':
    main()
