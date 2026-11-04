"""Audit candidate chemical identities against supplied training tables, without inference.

This is a source-curation aid, not experimental admission or evidence of prediction
accuracy. Wildcard terms are retained separately; an abstract term must never cause
the other, concrete molecules in its training row to disappear from the audit.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from fractions import Fraction
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
import re

from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.MolStandardize import rdMolStandardize

FP = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048, includeChirality=True)
TERM = re.compile(r'^(\d+(?:\.\d+)?)\s+(.+)$')
UNCHARGER = rdMolStandardize.Uncharger()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def equation_parts(text):
    sides = text.split(' = ')
    if len(sides) > 2 or any(not side.strip() for side in sides):
        raise ValueError('Invalid training equation')
    result = []
    for side, sign in zip(sides, [-1, 1] if len(sides) == 2 else [1]):
        for value in side.split(' + '):
            match = TERM.fullmatch(value.strip())
            coefficient, smiles = (Fraction(match[1]), match[2]) if match else (Fraction(1), value.strip())
            if coefficient <= 0 or not smiles:
                raise ValueError('Empty or nonpositive term')
            result.append((sign*coefficient, smiles))
    return result


@lru_cache(maxsize=None)
def molecular_identity(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError('Unparseable SMILES')
    if any(atom.GetAtomicNum() == 0 for atom in mol.GetAtoms()):
        raise ValueError('Abstract wildcard structure')
    # This operational family key removes ordinary acid/base charge differences,
    # retains covalent connectivity/stereochemistry, and does not infer tautomers,
    # magnesium complexes, hydration pools or chemically unspecified stereocentres.
    mol = Chem.RemoveHs(UNCHARGER.uncharge(mol))
    return Chem.MolToSmiles(mol, isomericSmiles=True), mol


def signature(parts):
    counts = defaultdict(Fraction)
    for coefficient, smiles in parts:
        key, _ = molecular_identity(smiles)
        if key not in ('O', '[H+]', '[H]'):
            counts[key] += Fraction(coefficient)
    values = sorted((k, c) for k, c in counts.items() if c)
    if not values:
        return ()
    denominator = math.lcm(*(v.denominator for _, v in values))
    integers = [(k, int(v*denominator)) for k, v in values]
    divisor = math.gcd(*(abs(v) for _, v in integers))
    primitive = tuple((k, v//divisor) for k, v in integers)
    return min(primitive, tuple((k, -v) for k, v in primitive))


def training_domain(path):
    with Path(path).open(newline='') as stream:
        rows = list(csv.DictReader(stream))
    molecules, reactions, unresolved, templates = {}, set(), [], {}
    for index, row in enumerate(rows):
        try:
            parts = equation_parts(row['reaction'])
        except (KeyError, ValueError) as error:
            unresolved.append(dict(row=index, reason=str(error), description=row.get('description')))
            continue
        complete = True
        for _, smiles in parts:
            try:
                key, mol = molecular_identity(smiles)
                if key not in ('O', '[H+]', '[H]'):
                    entry = molecules.setdefault(key, dict(molecule=mol, rows=[], fingerprint=FP.GetFingerprint(mol)))
                    entry['rows'].append(index)
            except ValueError as error:
                complete = False
                unresolved.append(dict(row=index, smiles=smiles, reason=str(error), description=row.get('description')))
                raw = Chem.MolFromSmiles(smiles)
                if raw is not None and any(a.GetAtomicNum() == 0 for a in raw.GetAtoms()):
                    parameters = Chem.AdjustQueryParameters()
                    parameters.makeDummiesQueries = True
                    parameters.adjustDegree = False
                    query = Chem.AdjustQueryProperties(Chem.RemoveHs(raw), parameters)
                    templates[smiles] = dict(query=query, row=index)
        if complete and ' = ' in row['reaction']:
            reactions.add(signature(parts))
    if not molecules:
        raise ValueError('No concrete training molecules')
    return dict(path=str(path), sha256=sha(path), n_rows=len(rows), molecules=molecules,
                reactions=reactions, unresolved_terms=unresolved, wildcard_templates=templates)


def balance(parts, n_proton=0):
    atoms = defaultdict(Fraction)
    charge = Fraction(str(n_proton))
    for coefficient, smiles in parts:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None or any(a.GetAtomicNum() == 0 for a in mol.GetAtoms()):
            raise ValueError('Candidate must have concrete structures')
        coefficient = Fraction(str(coefficient))
        for atom in Chem.AddHs(mol).GetAtoms():
            atoms[atom.GetSymbol()] += coefficient
        charge += coefficient*Chem.GetFormalCharge(mol)
    atoms['H'] += Fraction(str(n_proton))
    if any(atoms.values()) or charge:
        raise ValueError('Candidate equation is not element/charge balanced')


def audit_candidate(candidate, domain):
    parts = [(Fraction(str(r['coefficient'])), r['smiles']) for r in candidate['species']]
    balance(parts, candidate.get('n_Hplus', 0))
    molecules = domain['molecules']
    keys, fingerprints = list(molecules), [v['fingerprint'] for v in molecules.values()]
    rows = []
    for row in candidate['species']:
        key, mol = molecular_identity(row['smiles'])
        if key in ('O', '[H+]', '[H]'):
            continue
        similarities = DataStructs.BulkTanimotoSimilarity(FP.GetFingerprint(mol), fingerprints)
        best = max(range(len(keys)), key=lambda i: similarities[i])
        plausible = [dict(smiles=s, row=t['row']) for s, t in domain['wildcard_templates'].items()
                     if mol.HasSubstructMatch(t['query'], useChirality=True)]
        rows.append(dict(name=row['name'], supplied_smiles=row['smiles'], family_key=key,
            concrete_family_present=key in molecules,
            matching_training_rows=molecules[key]['rows'] if key in molecules else [],
            nearest_concrete_family=keys[best], nearest_chiral_morgan_similarity=similarities[best],
            compatible_abstract_training_scaffolds=plausible))
    unknown_rows = sorted({r['row'] for r in domain['unresolved_terms']})
    absent = [r for r in rows if not r['concrete_family_present']]
    distant = [r for r in absent if r['nearest_chiral_morgan_similarity'] < .5]
    return dict(candidate=candidate['id'], chemically_balanced=True,
        concrete_training_equation_present=signature(parts) in domain['reactions'], molecules=rows,
        n_absent_concrete_families=len(absent), n_distant_absent_families=len(distant),
        meets_structural_screen_on_concrete_subset=bool(distant),
        training_rows_with_unresolved_terms=unknown_rows,
        primary_benchmark_novelty_certified=False,
        qualification='Concrete-subset chemistry screen only. Abstract labels, source ancestry, state-pool identity, actual scored structures and exact comparator/checkpoint ownership require separate adjudication.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('candidates', type=Path)
    parser.add_argument('--training', action='append', required=True, metavar='NAME=CSV')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    candidates = json.loads(args.candidates.read_text())['candidates']
    if len({r['id'] for r in candidates}) != len(candidates):
        raise ValueError('Duplicate candidate IDs')
    report = dict(schema='metag-independent-candidate-structure-audit-v1',
        source_sha256=sha(__file__), candidates_sha256=sha(args.candidates), rdkit_version=rdBase.rdkitVersion,
        fingerprint=dict(kind='Morgan', radius=2, bits=2048, includeChirality=True, threshold=.5),
        standardization='RDKit Uncharger; RemoveHs; canonical isomeric SMILES. No tautomer or solution-pool inference.',
        inference_performed=False, experimental_references_admitted=0, training={})
    for value in args.training:
        name, filename = value.split('=', 1)
        if not name or name in report['training']:
            raise ValueError('Distinct, nonempty training labels required')
        domain = training_domain(filename)
        report['training'][name] = dict(path=filename, sha256=domain['sha256'], n_rows=domain['n_rows'],
            n_concrete_families=len(domain['molecules']), n_concrete_equations=len(domain['reactions']),
            unresolved_terms=domain['unresolved_terms'], results=[audit_candidate(c, domain) for c in candidates])
    with args.output.open('x') as stream:
        json.dump(report, stream, indent=2, allow_nan=False); stream.write('\n')
    print(json.dumps({name: dict(n_rows=data['n_rows'], n_concrete_families=data['n_concrete_families'],
        n_unresolved_terms=len(data['unresolved_terms']), candidates=[dict(id=r['candidate'],
            absent=r['n_absent_concrete_families'], distant=r['n_distant_absent_families'])
            for r in data['results']]) for name, data in report['training'].items()}, indent=2))


if __name__ == '__main__':
    main()
