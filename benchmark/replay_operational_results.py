"""Reconstruct the two terminal operational contexts from archived raw records."""
from pathlib import Path
import argparse
import contextlib
import hashlib
import importlib.abc
import io
import json
import os
import runpy
import sys
import tempfile
import zipfile


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class NoInference(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split('.')[0] in {'torch', 'fairchem', 'cupy', 'pyscf'}:
            raise RuntimeError('Inference import forbidden in retained-record replay')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--packet', type=Path,
                        default=Path(__file__).with_name('modelseed_operational_result_validation'))
    args = parser.parse_args()
    packet = args.packet.resolve()
    project = Path(__file__).resolve().parents[1]
    manifest = json.loads((packet / 'manifest.json').read_text())
    assert {p.name for p in packet.iterdir() if p.is_file()} == set(manifest['files']) | {'manifest.json'}
    for name, row in manifest['files'].items():
        path = packet / name
        assert path.stat().st_size == row['bytes'] and digest(path) == row['sha256'], name
    for name, expected in manifest['dependencies_sha256'].items():
        assert digest(project / name) == expected, name
    expected_summary = json.loads((packet / 'summary.json').read_text())
    assert expected_summary['coverage']['population_numeric_interval'] is None
    assert expected_summary['independent_experimental_labels'] == 0
    with tempfile.TemporaryDirectory(prefix='metag-terminal-operational-') as temp:
        root = Path(temp) / 'MetaG'
        study = root / 'analysis/modelseed_operational_20261105_01'
        files_verified = 0
        for archive_name, entries in manifest['archive_members'].items():
            with zipfile.ZipFile(packet / archive_name) as archive:
                infos = archive.infolist()
                assert len({e.filename for e in infos}) == len(infos)
                assert {e.filename for e in infos} == set(entries)
                assert sum(e.file_size for e in infos) < 100_000_000
                for info in infos:
                    name = Path(info.filename)
                    assert not name.is_absolute() and '..' not in name.parts
                    assert (info.external_attr >> 16) & 0o170000 != 0o120000
                    if archive_name == 'source_and_inputs.zip':
                        if name.parts[0] == 'analysis':
                            target = study.joinpath(*name.parts[1:])
                        elif name.parts[0] == 'benchmark':
                            target = root / name
                        else:
                            target = study / name
                    elif archive_name == 'raw_records_and_outcomes.zip':
                        target = study / name
                    elif name.parts[0].startswith('assembly_final_'):
                        target = study / 'resource_recovery' / name
                    elif name.parts[0] == 'context':
                        target = study.joinpath(*name.parts[1:])
                    else:
                        target = study / name
                    target.parent.mkdir(parents=True, exist_ok=True)
                    target.write_bytes(archive.read(info))
                    row = entries[info.filename]
                    assert target.stat().st_size == row['bytes'] and digest(target) == row['sha256']
                    files_verified += 1

        sys.dont_write_bytecode = True
        sys.path[:] = [str(root)] + [p for p in sys.path if p and not Path(p).resolve().is_relative_to(project)]
        for name in list(sys.modules):
            if name.split('.')[0] in {'benchmark', 'metag'}:
                del sys.modules[name]
        sys.meta_path.insert(0, NoInference())
        active, denied = [True], []

        def guard(event, args):
            if not active[0]:
                return
            if event.startswith('socket.'):
                denied.append(event)
                raise RuntimeError('Network forbidden in retained-record replay')
            if event == 'open' and isinstance(args[0], (str, bytes)):
                path = Path(os.fsdecode(args[0])).resolve()
                if path.is_relative_to(project):
                    denied.append(str(path))
                    raise RuntimeError('Original project fallback forbidden: ' + str(path))

        sys.addaudithook(guard)
        contexts = [('original', False), ('secondary_water', True)]
        try:
            for label, secondary in contexts:
                output = study / ('replay_' + label)
                script = study / 'analyze_results_v3.py'
                sys.argv = [str(script), '--output', str(output)]
                if secondary:
                    sys.argv += ['--secondary-water']
                with contextlib.redirect_stdout(io.StringIO()):
                    runpy.run_path(str(script), run_name='__main__')
                expected = study / 'resource_recovery' / ('assembly_final_' + label)
                assert {p.name for p in expected.iterdir()} == {p.name for p in output.iterdir()}
                for file in expected.iterdir():
                    wanted = json.loads(file.read_text())
                    actual = json.loads((output / file.name).read_text())
                    if file.name == 'report.json':
                        wanted.pop('utc')
                        actual.pop('utc')
                    assert actual == wanted, (label, file.name)
            assert not denied
        finally:
            active[0] = False
        print(json.dumps(dict(status='passed', contexts_reconstructed=2,
                              archived_files_verified=files_verified,
                              original_project_reads_forbidden=True, network_forbidden=True,
                              inference_imports_forbidden=True, denied_attempts=denied,
                              original_estimates_each_policy=5,
                              secondary_water_estimates_each_policy=9,
                              population_inference_withheld=True,
                              independent_experimental_labels=0), indent=2))


if __name__ == '__main__':
    main()
