"""Replay source triage from its archived inputs; never infer or fetch references."""
from pathlib import Path
import argparse
import contextlib
import gzip
import hashlib
import importlib.abc
import io
import json
import runpy
import sys
import tempfile
import zipfile


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class NoInference(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split('.')[0] in {'torch', 'fairchem', 'cupy', 'pyscf'}:
            raise RuntimeError('Inference import forbidden in structural replay: ' + fullname)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--packet', type=Path,
                        default=Path(__file__).resolve().with_name('source_universe_review'))
    args = parser.parse_args()
    packet = args.packet.resolve()
    manifest = json.loads((packet / 'manifest.json').read_text())
    assert {p.name for p in packet.iterdir() if p.is_file()} == set(manifest['files']) | {'manifest.json'}
    for name, record in manifest['files'].items():
        path = packet / name
        assert path.stat().st_size == record['bytes'] and digest(path) == record['sha256'], name
    expected = {name: json.loads(gzip.decompress((packet / (name + '.json.gz')).read_bytes()))
                for name in ['source_metadata_screen', 'structural_screen']}
    review = json.loads((packet / 'source_review.json').read_text())
    assert review['independent_experimental_labels_admitted'] == 0
    assert not review['model_predictions_performed'] and not review['primary_design_changed']
    assert all(not s['primary_article_body_owned'] for s in review['sources'])
    assert review['source_metadata_counts'] == expected['source_metadata_screen']['counts']

    with tempfile.TemporaryDirectory(prefix='metag-source-universe-') as temp:
        root = Path(temp)
        with zipfile.ZipFile(packet / 'source_replay.zip') as archive:
            members = archive.infolist()
            assert len({e.filename for e in members}) == len(members)
            assert sum(e.file_size for e in members) < 30_000_000
            assert {e.filename for e in members} == set(manifest['replay_members_sha256'])
            for entry in members:
                name = Path(entry.filename)
                assert not name.is_absolute() and '..' not in name.parts
                assert (entry.external_attr >> 16) & 0o170000 != 0o120000
                path = root / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(archive.read(entry))
                assert digest(path) == manifest['replay_members_sha256'][entry.filename]

        project = Path(__file__).resolve().parents[1]
        sys.path[:] = [p for p in sys.path if p and not Path(p).resolve().is_relative_to(project)]
        sys.dont_write_bytecode = True
        sys.meta_path.insert(0, NoInference())
        for name in list(sys.modules):
            if name.split('.')[0] in {'benchmark', 'metag'}:
                del sys.modules[name]

        replay_active = [True]

        def guard(event, args):
            if not replay_active[0]:
                return
            if event.startswith('socket.'):
                raise RuntimeError('Network access forbidden in structural replay')
            if event == 'open' and isinstance(args[0], (str, bytes)):
                path = Path(args[0].decode() if isinstance(args[0], bytes) else args[0]).resolve()
                if path.is_relative_to(project):
                    raise RuntimeError('Original project access forbidden during archived replay: ' + str(path))

        sys.addaudithook(guard)
        study = root / 'MetaG/analysis/reference_universe_scout_20261105_01'
        try:
            for source, result in [('scout.py', 'source_metadata_screen'),
                                   ('structural_screen.py', 'structural_screen')]:
                with contextlib.redirect_stdout(io.StringIO()):
                    runpy.run_path(str(study / source), run_name='__main__')
                actual = json.loads((study / (result + '.json')).read_text())
                # Molecular canonicalization is version-sensitive: an exact replay
                # requires the archived RDKit version, rather than dropping fields.
                assert actual == expected[result], result
        finally:
            # shutil's dir_fd-relative cleanup paths are not reads by the replay.
            replay_active[0] = False
        print(json.dumps(dict(status='passed', reports_reconstructed=2,
                              original_project_reads_allowed=False, network_allowed=False,
                              inference_imports_allowed=False,
                              metadata_counts=actual.get('source_metadata_counts', review['source_metadata_counts']),
                              independently_admitted_references=0), indent=2))


if __name__ == '__main__':
    main()
