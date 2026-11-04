"""Check that the manuscript's referenced historical evidence is present and intact.

Run with a staged manuscript directory whose parent contains the referenced
analysis and artifacts directories::

    python benchmark/check_manuscript_source_manifest.py /path/to/manuscript

This checks the eight named evidence manifests' original source files. It does not
check every numerical claim, prove scientific correctness, or verify current
unpublished calculations. The manuscript's own check_evidence.py is still
needed for macro and outcome checks.
"""
import argparse
import hashlib
import json
from pathlib import Path

MANIFESTS = ('revision_20261008.json', 'dev20_20261022.json')
PHYSICAL_MANIFEST = 'physical_components_20261011.json'
MAPPED_MANIFESTS = (PHYSICAL_MANIFEST, 'completed_domain_20261011.json',
                    'ester_solvent_checks_20261011.json', 'ensemble_refinement_20261011.json')
ROUTE_MANIFEST = 'route_policy_dev13.json'
GEOMETRY_MANIFEST = 'geometry_integrity.json'


def check_sources(manuscript):
    manuscript = Path(manuscript).resolve()
    root = manuscript.parent
    expected = {}
    manifests = [(name, json.loads((manuscript/'evidence'/name).read_text())['sources'])
                 for name in MANIFESTS]
    for name in MAPPED_MANIFESTS:
        files = json.loads((manuscript/'evidence'/name).read_text())['files']
        manifests.append((name, list(files.values())))
    route = json.loads((manuscript/'evidence'/ROUTE_MANIFEST).read_text())
    manifests.append((ROUTE_MANIFEST, [{'path':key,'sha256':digest} for key,digest in route['files'].items()]))
    geometry = json.loads((manuscript/'evidence'/GEOMETRY_MANIFEST).read_text())
    manifests.append((GEOMETRY_MANIFEST, geometry['sources']))
    for section in ('component_replay','exchange_replay'):
        values=geometry['records'][section]['sources_sha256']
        manifests.append((GEOMETRY_MANIFEST, [{'path':key,'sha256':digest} for key,digest in values.items()]))
    for name, rows in manifests:
        for row in rows:
            rel = Path(row['path'])
            # The historical physical-control manifest records an old absolute
            # project root. Resolve only its known project-relative suffix.
            if name in (*MAPPED_MANIFESTS, ROUTE_MANIFEST, GEOMETRY_MANIFEST) and 'MetaG_JCIM' in rel.parts:
                rel = Path(*rel.parts[rel.parts.index('MetaG_JCIM') + 1:])
            if rel.is_absolute() or '..' in rel.parts or not rel.parts:
                raise ValueError(f'Unsafe evidence path: {rel}')
            path = (root/rel).resolve()
            if not path.is_relative_to(root):
                raise ValueError(f'Evidence path escapes release root: {rel}')
            digest = row['sha256']
            if len(digest) != 64 or any(c not in '0123456789abcdef' for c in digest):
                raise ValueError(f'Invalid SHA-256 for {rel}')
            if str(rel) in expected and expected[str(rel)] != digest:
                raise ValueError(f'Conflicting SHA-256 values for {rel}')
            expected[str(rel)] = digest
    missing, changed = [], []
    for rel, digest in sorted(expected.items()):
        path = root/rel
        if not path.is_file():
            missing.append(rel)
        elif hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            changed.append(rel)
    return {'scope':'Original sources listed in eight historical manuscript evidence manifests only.',
            'expected':len(expected), 'verified':len(expected)-len(missing)-len(changed),
            'missing':missing, 'changed':changed}


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('manuscript', type=Path, help='manuscript directory in the staged study archive')
    args=parser.parse_args(argv)
    result=check_sources(args.manuscript)
    print(json.dumps(result, indent=2))
    return 1 if result['missing'] or result['changed'] else 0


if __name__=='__main__':
    raise SystemExit(main())
