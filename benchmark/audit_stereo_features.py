"""Replay the frozen stereo-feature audit without inference or experimental enrollment.

The default replay uses RDKit and independently reconstructs the archived arrays.
--actual-featurizer additionally executes the pinned Torch/PyG feature function on
CPU. Neither mode imports the comparator API or loads its models.
"""
from __future__ import annotations

import argparse
import ast
from collections import Counter, defaultdict
import gzip
import hashlib
import importlib.util
from itertools import combinations
import json
import math
from pathlib import Path
import sys

from rdkit import Chem, rdBase
from rdkit.Chem.rdchem import BondType, HybridizationType


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError('Invalid frozen source SMILES')
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.AddHs(mol)


def stereo(mol):
    mirror, flat = Chem.Mol(mol), Chem.Mol(mol)
    for atom in mirror.GetAtoms():
        atom.InvertChirality()
    Chem.AssignStereochemistry(mirror, cleanIt=True, force=True)
    Chem.RemoveStereochemistry(flat)
    return dict(exact=Chem.MolToSmiles(mol, isomericSmiles=True),
                mirror=Chem.MolToSmiles(mirror, isomericSmiles=True),
                flat=Chem.MolToSmiles(flat, isomericSmiles=True),
                unspecified=[dict(type=str(s.type), center=int(s.centeredOn),
                                  specified=str(s.specified))
                             for s in Chem.FindPotentialStereo(mol)
                             if str(s.specified) != 'Specified'],
                isotopes=any(a.GetIsotope() for a in mol.GetAtoms()))


def one_hot(size, index):
    result = [0.] * size
    result[int(index)] = 1.
    return result


def arrays(mol):
    """Independent reconstruction of the pinned default feature definition.

    Negative charge indices intentionally retain the original Python indexing
    convention. This audit does not repair or change the comparator features.
    """
    atoms = [one_hot(119, a.GetAtomicNum())
             + one_hot(len(HybridizationType.values), a.GetHybridization())
             + one_hot(2, a.GetIsAromatic())
             + one_hot(9, a.GetFormalCharge()) for a in mol.GetAtoms()]
    edges = []
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        feature = one_hot(len(BondType.values) + 1, bond.GetBondType())
        edges.extend([(a, b, feature), (b, a, feature)])
    edges.sort(key=lambda edge: edge[:2])
    return dict(x=atoms, edge_index=[[edge[i] for edge in edges] for i in (0, 1)],
                edge_attr=[edge[2] for edge in edges])


def array_digest(data):
    return hashlib.sha256(json.dumps(data, sort_keys=True, allow_nan=False).encode()).hexdigest()


def actual_function(packet):
    sys.dont_write_bytecode = True
    import torch
    import torch_geometric
    assert not torch.cuda.is_initialized()
    torch.set_num_threads(1)
    torch.set_num_interop_threads(1)
    spec = importlib.util.spec_from_file_location(
        'frozen_comparator_features', packet / 'owned_source/dGbyG/utils/NNFunc.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    def evaluate(mol, **kwargs):
        data = module.mol_to_graph_data(mol, **kwargs)
        assert all(t.device.type == 'cpu' for t in (data.x, data.edge_index, data.edge_attr))
        order = sorted(range(data.edge_index.shape[1]),
                       key=lambda i: tuple(data.edge_index[:, i].tolist()))
        assert not torch.cuda.is_initialized()
        return dict(x=data.x.tolist(), edge_index=data.edge_index[:, order].tolist(),
                    edge_attr=data.edge_attr[order].tolist())
    return evaluate, dict(torch=torch.__version__, PyG=torch_geometric.__version__,
                          CUDA_initialized=False)


def audit(packet, *, actual_featurizer=False):
    packet = Path(packet)
    manifest = read(packet / 'manifest.json')
    for name, identity in manifest['files'].items():
        assert not Path(name).is_absolute() and '..' not in Path(name).parts
        path = packet / name
        assert not path.is_symlink() and path.stat().st_size == identity['bytes'], name
        assert sha(path) == identity['sha256'], name
    protocol, report = read(packet / 'protocol.json'), read(packet / 'report.json')
    assert sha(packet / 'protocol.json') == report['protocol_sha256']
    assert sha(packet / 'source_run.py') == report['source_sha256']
    for name, digest in protocol['inputs'].items():
        assert sha(packet / 'inputs' / Path(name).name) == digest, name
    for name, digest in protocol['sources'].items():
        assert sha(packet / 'owned_source' / name) == digest, name
    parent_protocol = read(packet / 'inputs/protocol.json')
    for name, digest in protocol['sources'].items():
        assert parent_protocol['source_assets'][name]['sha256'] == digest
    # The exact API in the native experiment calls the default feature function.
    tree = ast.parse((packet / 'owned_source/dGbyG/api/utils.py').read_text())
    function = next(n for n in tree.body if isinstance(n, ast.FunctionDef)
                    and n.name == 'predict_standard_dGf_prime')
    calls = [n for n in ast.walk(function) if isinstance(n, ast.Call)
             and isinstance(n.func, ast.Name) and n.func.id == 'mol_to_graph_data']
    assert len(calls) == 1 and len(calls[0].args) == 1 and not calls[0].keywords
    feature_tree = ast.parse((packet / 'owned_source/dGbyG/utils/NNFunc.py').read_text())
    function = next(n for n in feature_tree.body if isinstance(n, ast.FunctionDef)
                    and n.name == 'mol_to_graph_data')
    defaults = [ast.literal_eval(n) for n in function.args.defaults]
    assert defaults == [['atomic number', 'hybridization', 'aromaticity', 'charge'], ['bond type']]
    assert report['default_atom_features'] == defaults[0]
    assert report['default_bond_features'] == defaults[1]

    feature_function, runtime = (actual_function(packet) if actual_featurizer
                                 else (None, dict(mode='independent RDKit reconstruction')))
    stored_arrays = json.loads(gzip.decompress((packet / 'feature_arrays.json.gz').read_bytes()))
    seen_arrays = set()

    def check_arrays(mol):
        result = arrays(mol)
        if feature_function is not None:
            assert feature_function(mol) == result
        return result

    def pair(left_smiles, right_smiles, archived):
        left, right = molecule(left_smiles), molecule(right_smiles)
        a, b = stereo(left), stereo(right)
        assert a == archived['left'] and b == archived['right']
        assert a['flat'] == b['flat'] and not a['isotopes'] and not b['isotopes']
        flat_left, flat_right = Chem.Mol(left), Chem.Mol(right)
        Chem.RemoveStereochemistry(flat_left)
        Chem.RemoveStereochemistry(flat_right)
        mapping = flat_right.GetSubstructMatch(flat_left, useChirality=False)
        assert len(mapping) == left.GetNumAtoms() == right.GetNumAtoms()
        assert list(mapping) == archived['right_to_left_atom_order']
        right = Chem.RenumberAtoms(right, list(mapping))
        left_data, right_data = check_arrays(left), check_arrays(right)
        assert left_data == right_data
        digest = array_digest(left_data)
        assert digest == archived['graph_feature_sha256']
        assert stored_arrays[digest] == left_data
        seen_arrays.add(digest)
        category = ('same_canonical_state' if a['exact'] == b['exact'] else
                    'unresolved_stereo' if a['unspecified'] or b['unspecified'] else
                    'enantiomer_pair' if a['mirror'] == b['exact'] else
                    'fully_specified_non_enantiomer_pair')
        assert archived['category'] == category and archived['feature_arrays_identical']
        return category

    # Counterfactual controls do not use experimental outcome values.
    control = molecule('C[C@H](O)C(=O)O')
    inverted = Chem.Mol(control)
    for atom in inverted.GetAtoms():
        atom.InvertChirality()
    assert stereo(control)['exact'] != stereo(inverted)['exact']
    assert check_arrays(control) == check_arrays(inverted)
    ez, opposite = molecule('Cl/C=C/Cl'), molecule(r'Cl/C=C\Cl')
    assert stereo(ez)['exact'] != stereo(opposite)['exact']
    assert check_arrays(ez) == check_arrays(opposite)
    assert check_arrays(molecule('CCO')) != check_arrays(molecule('CCN'))
    if feature_function is not None:
        optional = defaults[0] + ['chirality']
        assert feature_function(control, atom_features=optional) != feature_function(inverted, atom_features=optional)
    assert all(report['controls'].values())

    formation = read(packet / 'inputs/structural_review.json')['pairs']
    assert len(formation) == len(report['formation_pairs'])
    formation_counts = Counter()
    for source, archived in zip(formation, report['formation_pairs'], strict=True):
        for side in ('left', 'right'):
            assert archived[side + '_name'] == source[side]['name']
            assert archived[side + '_training_row'] == source[side]['row']
        formation_counts[pair(source['left']['smiles'], source['right']['smiles'], archived)] += 1
    assert dict(formation_counts) == report['formation_pair_counts']

    raw = read(packet / 'inputs/raw_outputs.json')
    molecules = raw['molecular_vectors']
    groups = defaultdict(list)
    for key, value in molecules.items():
        groups[stereo(molecule(value['normalized_SMILES']))['flat']].append(key)
    candidates = [(a, b) for group in groups.values() for a, b in combinations(group, 2)
                  if stereo(molecule(molecules[a]['normalized_SMILES']))['exact'] !=
                  stereo(molecule(molecules[b]['normalized_SMILES']))['exact']]
    assert candidates == [(r['left_molecule'], r['right_molecule']) for r in report['native_pairs']]
    native_counts = Counter()
    for (a, b), archived in zip(candidates, report['native_pairs'], strict=True):
        native_counts[pair(molecules[a]['normalized_SMILES'], molecules[b]['normalized_SMILES'], archived)] += 1
        differences = [x - y for x, y in zip(molecules[a]['members_kj'], molecules[b]['members_kj'], strict=True)]
        assert len(differences) == 100 and differences == archived['member_difference_kj']
        assert all(value == 0 for value in differences)
        assert archived['native_mean_difference_kj'] == archived['native_population_SD_kj'] == 0
    assert dict(native_counts) == report['native_pair_counts']
    assert seen_arrays == set(stored_arrays)

    selected = []
    for row in raw['raw_reactions']:
        exact, flat = defaultdict(float), defaultdict(float)
        for term in row['source_terms']:
            info = stereo(molecule(molecules[term['molecule']]['normalized_SMILES']))
            exact[info['exact']] += term['coefficient']
            flat[info['flat']] += term['coefficient']
        if not any(flat.values()) and any(exact.values()):
            selected.append(row)
    assert [r['reaction'] for r in selected] == [r['reaction'] for r in report['graph_cancelled_native_reactions']]
    for source, archived in zip(selected, report['graph_cancelled_native_reactions'], strict=True):
        for field in ('reaction', 'author675_row', 'source_terms'):
            assert source[field] == archived[field]
        vector = [math.fsum(t['coefficient'] * molecules[t['molecule']]['members_kj'][k]
                            for t in source['source_terms']) for k in range(100)]
        assert vector == source['members_native_kj'] == archived['native_neural_members_kj']
        assert all(value == 0 for value in vector)
        assert archived['neural_mean_kj'] == archived['neural_population_SD_kj'] == 0
        assert archived['external_condition_bridge_kj'] == source['condition_bridge_kj']

    design, joined = read(packet / 'retrospective_join_design.json'), read(packet / 'retrospective_cases.json')
    assert sha(packet / 'inputs/native_inputs.json') == design['source_sha256']
    assert design['selected_reactions'] == [r['reaction'] for r in selected]
    assert [r['reaction'] for r in joined['rows']] == design['selected_reactions']
    native_rows = {r['reaction']: r for r in raw['rows']}
    for row in joined['rows']:
        source = native_rows[row['reaction']]
        for key, source_key in [('reference', 'reference'), ('historical_metag', 'metag'), ('comparator', 'comparator')]:
            assert row[key] == source[source_key]
        assert row['native_ensemble_population_SD'] == source['spread'] == 0
        for key in ('historical_metag', 'comparator', 'zero'):
            assert row['absolute_errors'][key] == abs((0 if key == 'zero' else row[key]) - row['reference'])
    maes = {key: math.fsum(r['absolute_errors'][key] for r in joined['rows']) / len(joined['rows'])
            for key in ('historical_metag', 'comparator', 'zero')}
    for key, value in maes.items():
        assert abs(value - joined['mean_absolute_error'][key]) < 1e-12
    assert joined['independent_labels_admitted'] == report['independent_primary_labels_admitted'] == 0
    assert joined['source_or_condition_limits_unchanged']
    assert not report['training'] and report['weights_loaded'] == report['neural_forward_calls'] == 0
    assert report['network_forbidden'] and not report['network_attempts']
    assert report['atom_feature_dimension'] == len(next(iter(stored_arrays.values()))['x'][0]) == 139
    return dict(status='passed', runtime=dict(rdkit=rdBase.rdkitVersion, **runtime),
                authenticated_files=len(manifest['files']),
                formation_pair_counts=dict(formation_counts), native_pair_counts=dict(native_counts),
                graph_cancelled_reactions=[r['reaction'] for r in selected],
                independent_arrays_replayed=len(seen_arrays), actual_featurizer_executed=actual_featurizer,
                retrospective_mean_absolute_error=maes, independent_labels_admitted=0,
                independent_accuracy_advantage_demonstrated=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--packet', type=Path,
                        default=Path(__file__).with_name('stereo_feature_validation'))
    parser.add_argument('--actual-featurizer', action='store_true')
    args = parser.parse_args()
    def no_network(event, _args):
        if event in {'socket.connect', 'socket.getaddrinfo'}:
            raise RuntimeError('Network is forbidden in the feature audit')
    sys.addaudithook(no_network)
    print(json.dumps(audit(args.packet, actual_featurizer=args.actual_featurizer), indent=2))


if __name__ == '__main__':
    main()
