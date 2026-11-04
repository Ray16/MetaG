"""Covalent source-identity checks; no solution-state or benchmark admission."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
from rdkit import Chem, rdBase
from rdkit.Chem import rdMolDescriptors

canon = lambda m: Chem.MolToSmiles(m, isomericSmiles=True)
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()


def sugar_roles(molecule):
    rings = [set(r) for r in molecule.GetRingInfo().AtomRings() if len(r) == 5 and
             sorted(molecule.GetAtomWithIdx(i).GetAtomicNum() for i in r) == [6, 6, 6, 6, 8]]
    if len(rings) != 1:
        raise ValueError('One pentofuranose ring required')
    ring = rings[0]
    oxygen = next(i for i in ring if molecule.GetAtomWithIdx(i).GetAtomicNum() == 8)
    first = [i for i in ring if oxygen in [a.GetIdx() for a in molecule.GetAtomWithIdx(i).GetNeighbors()]
             and any(a.GetAtomicNum() == 7 and a.GetIdx() not in ring
                     for a in molecule.GetAtomWithIdx(i).GetNeighbors())]
    if len(first) != 1:
        raise ValueError('One N-glycosidic anomeric carbon required')
    first = first[0]
    second = next(a.GetIdx() for a in molecule.GetAtomWithIdx(first).GetNeighbors()
                  if a.GetIdx() in ring and a.GetIdx() != oxygen)
    third = next(a.GetIdx() for a in molecule.GetAtomWithIdx(second).GetNeighbors()
                 if a.GetIdx() in ring and a.GetIdx() != first)
    fourth = next(a.GetIdx() for a in molecule.GetAtomWithIdx(third).GetNeighbors()
                  if a.GetIdx() in ring and a.GetIdx() != second)
    return ring, (first, second, third, fourth)


def check_sugar(molecule, deoxy):
    if Chem.GetFormalCharge(molecule) != 0 or len(Chem.GetMolFrags(molecule)) != 1:
        raise ValueError('Source figure uses a single neutral reference drawing')
    if any(a.GetIsotope() for a in molecule.GetAtoms()):
        raise ValueError('Isotope-specific source not declared')
    ring, roles = sugar_roles(molecule)
    labels = dict(Chem.FindMolChiralCenters(molecule, includeUnassigned=True))
    if '?' in labels.values():
        raise ValueError('Source sugar stereochemistry is unspecified')
    actual = [labels.get(i, 'achiral') for i in roles]
    expected = ['R', 'achiral' if deoxy else 'R', 'S', 'R']
    if actual != expected:
        raise ValueError('Source betaD sugar configuration mismatch')
    hydroxyls = [a.GetIdx() for a in molecule.GetAtomWithIdx(roles[1]).GetNeighbors()
                 if a.GetIdx() not in ring and a.GetAtomicNum() == 8 and a.GetDegree() == 1]
    if len(hydroxyls) != (0 if deoxy else 1):
        raise ValueError('Source2prime hydroxyl identity mismatch')
    return dict(roles=dict(zip(['C1prime', 'C2prime', 'C3prime', 'C4prime'], roles)), CIP=actual)


def remove_2prime_OH(parent):
    check_sugar(parent, False)
    ring, roles = sugar_roles(parent)
    hydroxyl = [a.GetIdx() for a in parent.GetAtomWithIdx(roles[1]).GetNeighbors()
                if a.GetIdx() not in ring and a.GetAtomicNum() == 8 and a.GetDegree() == 1]
    molecule = Chem.RWMol(parent)
    carbon = molecule.GetAtomWithIdx(roles[1])
    carbon.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    carbon.SetNumExplicitHs(2)
    carbon.SetNoImplicit(True)
    molecule.RemoveAtom(hydroxyl[0])
    result = molecule.GetMol()
    Chem.SanitizeMol(result)
    check_sugar(result, True)
    return result


def substitute(parent, site, group):
    pattern = '[n]1[cH][cH][c](=O)[nH][c]1=O' if site == 'pyrimidine5' else '[cH;r6]([n;r6])[n;r6]'
    if site not in ('pyrimidine5', 'purine2'):
        raise ValueError('Undeclared source substitution site')
    matches = parent.GetSubstructMatches(Chem.MolFromSmarts(pattern))
    positions = {hit[2] if site == 'pyrimidine5' else hit[0] for hit in matches}
    if len(positions) != 1:
        raise ValueError('Source substitution site is not unique')
    position = positions.pop()
    molecule = Chem.RWMol(parent)
    molecule.GetAtomWithIdx(position).SetNumExplicitHs(0)
    if group in ('F', 'Cl', 'Br', 'I', 'N', 'C'):
        attached = molecule.AddAtom(Chem.Atom(group))
        molecule.AddBond(position, attached, Chem.BondType.SINGLE)
    elif group == 'C#C':
        first = molecule.AddAtom(Chem.Atom('C'))
        second = molecule.AddAtom(Chem.Atom('C'))
        molecule.AddBond(position, first, Chem.BondType.SINGLE)
        molecule.AddBond(first, second, Chem.BondType.TRIPLE)
    else:
        raise ValueError('Undeclared source substituent')
    result = molecule.GetMol()
    Chem.SanitizeMol(result)
    return result


def audit(directory):
    directory = Path(directory)
    source = json.loads((directory/'source_compounds.json').read_text())
    if len(source['compounds']) != 24 or len({r['code'] for r in source['compounds']}) != 24:
        raise ValueError('Complete24-compound source panel required')
    molecules, records = {}, []
    for row in source['compounds']:
        path = directory/row['selected_property_file']
        if sha(path) != row['selected_property_sha256']:
            raise ValueError('Cached property identity changed')
        candidates = json.loads(path.read_text())['PropertyTable']['Properties']
        if len(candidates) != 1:
            raise ValueError('Unresolved multiple database identities')
        value = candidates[0]
        molecule = Chem.MolFromSmiles(value.get('SMILES', value.get('IsomericSMILES')))
        if molecule is None:
            raise ValueError('Invalid property structure')
        formula = rdMolDescriptors.CalcMolFormula(molecule)
        if formula != value['MolecularFormula'] or Chem.MolToInchiKey(molecule) != value['InChIKey']:
            raise ValueError('Property formula or InChIKey mismatch')
        sugar = check_sugar(molecule, row['deoxy'])
        molecules[row['code']] = molecule
        records.append(dict(code=row['code'], SI_compound_number=row['SI_compound_number'],
            PubChem_CID=value['CID'], InChIKey=value['InChIKey'], formula=formula,
            canonical_isomeric_smiles=canon(molecule), sugar=sugar,
            reference_drawing_charge=0, solution_state_pool_admitted=False))
    checks = []
    for parent, child in [('Urd','dUrd'), ('Ado','dAdo'), ('Guo','dGuo'), ('Ino','dIno')]:
        if canon(remove_2prime_OH(molecules[parent])) != canon(molecules[child]):
            raise ValueError('Source2prime-deoxy pair mismatch: '+child)
        checks.append(dict(kind='remove_only2primeOH', parent=parent, compound=child, matched=True))
    for row in source['compounds']:
        if 'substitution' not in row:
            continue
        change = row['substitution']
        expected = substitute(molecules[change['parent']], change['site'], change['group'])
        if canon(expected) != canon(molecules[row['code']]):
            raise ValueError('Source substituent or sugar identity mismatch: '+row['code'])
        checks.append(dict(kind='substitution_with_fixed_sugar', compound=row['code'], **change, matched=True))
    return dict(schema='metag-nucleoside-source-identity-v1',
        source_sha256=sha(Path(__file__)), input_sha256=sha(directory/'source_compounds.json'),
        RDKit_version=rdBase.rdkitVersion, n_compounds=24, n_absolute_sugar_checks=24,
        n_relational_checks=len(checks), records=records, relational_checks=checks,
        source_specific_resolution=source['source_specific_resolution'],
        source_numeric_label_issue=source['source_numeric_label_issue'],
        scope='Neutral covalent reference structures matched to sourceFigureS1 substitutions and betaD sugar configurations. Tautomer/protonation/Mg pools, source-condition normalization, development exposure and comparator/checkpoint overlap require separate adjudication.',
        experimental_references_admitted=0, predictions_performed=False, primary_benchmark_launch_allowed=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    report = audit(args.directory)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps({k:report[k] for k in ['n_compounds', 'n_absolute_sugar_checks', 'n_relational_checks', 'experimental_references_admitted']}, indent=2))


if __name__ == '__main__':
    main()
