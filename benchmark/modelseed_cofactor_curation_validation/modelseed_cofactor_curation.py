"""Explicit research curation of four source-owned NAD(P)(H) annotations.

This is not a default scoring route, tautomer canonicalizer or thermodynamic
equivalence rule. A separate annotation profile must own exact source/target
records. Formula/proton accounting preserves a reaction's original defects.
"""
from collections import Counter, defaultdict
from fractions import Fraction
import copy

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import networkx as nx


def canonical(smiles):
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None or not molecule.GetNumAtoms():
        raise ValueError('Invalid molecular identity')
    return Chem.MolToSmiles(molecule, isomericSmiles=True)


def composition(smiles):
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError('Invalid molecular identity')
    atoms = Counter((a.GetAtomicNum(), a.GetIsotope()) for a in Chem.AddHs(molecule).GetAtoms())
    return dict(atoms), Chem.GetFormalCharge(molecule)


def _heavy_graph(molecule):
    Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
    graph = nx.Graph()
    for atom in molecule.GetAtoms():
        if atom.GetAtomicNum() == 1:
            raise ValueError('Explicit/isotopic hydrogens require a separate curation policy')
        graph.add_node(atom.GetIdx(), element=atom.GetAtomicNum(), isotope=atom.GetIsotope(),
                       carbon_stereo=atom.GetProp('_CIPCode') if atom.GetAtomicNum() == 6
                       and atom.HasProp('_CIPCode') else None)
    graph.add_edges_from((b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in molecule.GetBonds())
    return graph


def review_conversion(source_charge, source_smiles, target_charge, target_smiles):
    """Verify heavy connectivity/carbon stereo and exact ordinary-proton delta.

    Bond orders, formal charges and hydrogens are explicitly compared, rather
    than asserted equivalent. No energy or measured equilibrium enters review.
    """
    source, target = Chem.MolFromSmiles(source_smiles), Chem.MolFromSmiles(target_smiles)
    if source is None or target is None:
        raise ValueError('Invalid molecular identity')
    if Chem.GetFormalCharge(source) != source_charge or Chem.GetFormalCharge(target) != target_charge:
        raise ValueError('Declared charge does not match molecular structure')
    sf, sq = composition(source_smiles)
    tf, tq = composition(target_smiles)
    elements = set(sf) | set(tf)
    delta = {key: tf.get(key, 0)-sf.get(key, 0) for key in elements if tf.get(key, 0) != sf.get(key, 0)}
    if any(key != (1, 0) for key in delta):
        raise ValueError('Curation changes nonproton elemental/isotope composition')
    dh = delta.get((1, 0), 0)
    if tq-sq != dh:
        raise ValueError('Curation is not balanced by ordinary protons')
    matcher = nx.algorithms.isomorphism.GraphMatcher(_heavy_graph(source), _heavy_graph(target),
        node_match=nx.algorithms.isomorphism.categorical_node_match(
            ['element', 'isotope', 'carbon_stereo'], [None, None, None]))
    candidates = []
    for index, mapping in enumerate(matcher.isomorphisms_iter()):
        if index >= 4096:
            raise ValueError('Ambiguous graph mapping exceeds declared review limit')
        stereo_changed = False
        for bond in source.GetBonds():
            other = target.GetBondBetweenAtoms(mapping[bond.GetBeginAtomIdx()], mapping[bond.GetEndAtomIdx()])
            declared = {Chem.BondStereo.STEREONONE, Chem.BondStereo.STEREOANY}
            if (bond.GetStereo() not in declared or other.GetStereo() not in declared) and bond.GetStereo() != other.GetStereo():
                stereo_changed = True
                break
        if stereo_changed:
            continue
        atoms = []
        for a in source.GetAtoms():
            b = target.GetAtomWithIdx(mapping[a.GetIdx()])
            if (a.GetFormalCharge(), a.GetTotalNumHs()) != (b.GetFormalCharge(), b.GetTotalNumHs()):
                atoms.append(dict(source_index=a.GetIdx(), target_index=b.GetIdx(), element=a.GetSymbol(),
                    source_charge=a.GetFormalCharge(), target_charge=b.GetFormalCharge(),
                    source_hydrogens=a.GetTotalNumHs(), target_hydrogens=b.GetTotalNumHs()))
        bonds = []
        for a in source.GetBonds():
            b = target.GetBondBetweenAtoms(mapping[a.GetBeginAtomIdx()], mapping[a.GetEndAtomIdx()])
            if a.GetBondType() != b.GetBondType():
                bonds.append(dict(source_atoms=[a.GetBeginAtomIdx(), a.GetEndAtomIdx()],
                    target_atoms=[b.GetBeginAtomIdx(), b.GetEndAtomIdx()],
                    source_order=str(a.GetBondType()), target_order=str(b.GetBondType())))
        candidates.append(((len(atoms)+len(bonds), tuple(mapping[i] for i in range(source.GetNumAtoms()))),
                           dict(atom_mapping={str(i): mapping[i] for i in sorted(mapping)},
                                changed_atoms=atoms, changed_bonds=bonds)))
    if not candidates:
        raise ValueError('Heavy connectivity/isotope/declared carbon stereochemistry mismatch')
    evidence = min(candidates, key=lambda item: item[0])[1]
    carbon_centers = sum(a.GetAtomicNum() == 6 and a.HasProp('_CIPCode') for a in source.GetAtoms())
    return dict(evidence, source_canonical=canonical(source_smiles), target_canonical=canonical(target_smiles),
        source_formula=rdMolDescriptors.CalcMolFormula(source), target_formula=rdMolDescriptors.CalcMolFormula(target),
        source_charge=sq, target_charge=tq, delta_H=dh, delta_charge=tq-sq,
        declared_carbon_stereocenters_preserved=carbon_centers,
        energy_equivalence_established=False,
        scope='Explicit source-identity curation; heavy connectivity, carbon stereo and proton accounting only.')


def validate_registry(registry, compounds):
    """Fail if a pinned source row has changed; never select by a common name."""
    checked = {}
    for cid, entry in registry.items():
        row = compounds[cid]
        expected = entry['source']
        for field in ['smiles', 'charge', 'formula', 'inchikey']:
            if row[field] != expected[field]:
                raise ValueError('Pinned source annotation changed: ' + cid + ':' + field)
        if entry['chebi_id'].removeprefix('CHEBI:') not in row['aliases'].split('ChEBI: ')[-1].split('|')[0].split('; '):
            raise ValueError('Source ChEBI cross-reference does not authorize curation identity')
        target = entry['target']
        review = review_conversion(int(expected['charge']), expected['smiles'], target['charge'], target['smiles'])
        if review['delta_H'] != entry['delta_H']:
            raise ValueError('Registry proton accounting changed')
        checked[cid] = copy.deepcopy(entry)
    return checked


def rebase_equation(terms, registry):
    """Change only explicit local H+ coefficients required by owned conversions.

    Positive coefficients denote products. If target=S+delta_H*H+, the
    replacement column must subtract sum(nu*delta_H) from the H+ coefficient.
    This is a basis conversion; an originally unbalanced equation stays so.
    """
    output = {key: Fraction(value) for key, value in terms.items() if value}
    changes = defaultdict(Fraction)
    for (cid, compartment), coefficient in terms.items():
        if cid in registry:
            delta = registry[cid]['delta_H']
            if isinstance(delta, bool) or not isinstance(delta, int):
                raise ValueError('Registry delta_H must be an exact integer')
            changes[compartment] -= Fraction(coefficient) * delta
    for compartment, correction in changes.items():
        key = ('cpd00067', compartment)
        output[key] = output.get(key, Fraction(0)) + correction
    output = {key: value for key, value in output.items() if value}
    return output, {compartment: str(value) for compartment, value in sorted(changes.items()) if value}


def column_balance(terms, compounds):
    """Exact literal residual; never drop negative counts or infer reservoirs."""
    elements = defaultdict(Fraction)
    charge = Fraction(0)
    for (cid, _), coefficient in terms.items():
        if cid == 'cpd00067':
            formula, q = {(1, 0): 1}, 1
        else:
            row = compounds[cid]
            formula, q = composition(row['smiles'])
            if q != int(row['charge']):
                raise ValueError('Source charge mismatch: ' + cid)
        for key, count in formula.items():
            elements[key] += Fraction(coefficient) * count
        charge += Fraction(coefficient) * q
    return {str(k): str(v) for k, v in sorted(elements.items()) if v}, str(charge)
