"""Audit dGbyG graph stereochemistry on simple one-to-one reaction pairs.

Run in the dGbyG environment; --infer executes its published API on CPU. This
examines the supplied molecular states, not an equilibrating solution pool.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model-root', type=Path, required=True)
    parser.add_argument('--reactions', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--infer', action='store_true', help='also run the official published ensemble on CPU')
    args = parser.parse_args(argv)
    sys.path.insert(0, str(args.model_root.resolve()))
    import networkx as nx
    import torch
    from rdkit import Chem, rdBase
    from rdkit.Chem import rdChemReactions  # Populate the attribute used in dGbyG API annotations.
    from dGbyG.Chemistry import _Compound
    from dGbyG.utils.NNFunc import mol_to_graph_data
    torch.set_num_threads(1)
    if args.infer:
        # The upstream API otherwise selects cuda:0 merely on import.
        torch.cuda.is_available = lambda: False
        from dGbyG.api.utils import predict_standard_dGf_prime

    def graph(data):
        g = nx.Graph()
        for i, features in enumerate(data.x.tolist()):
            g.add_node(i, features=tuple(features))
        for (i, j), features in zip(data.edge_index.T.tolist(), data.edge_attr.tolist()):
            g.add_edge(i, j, features=tuple(features))
        return g

    reactions = json.loads(args.reactions.read_text())
    records = {}
    for rid, rx in reactions.items():
        parts = list(rx['species'].values())
        if len(parts) != 2 or sorted(float(p[0]) for p in parts) != [-1, 1]:
            raise ValueError(f'{rid}: requires one reactant and one product')
        molecules = [_Compound(Chem.MolFromSmiles(p[2])).mol for p in parts]
        if any(m is None for m in molecules):
            raise ValueError(f'{rid}: invalid molecule')
        node_match = nx.algorithms.isomorphism.categorical_node_match('features', None)
        edge_match = nx.algorithms.isomorphism.categorical_edge_match('features', None)
        equal = nx.is_isomorphic(graph(mol_to_graph_data(molecules[0])), graph(mol_to_graph_data(molecules[1])),
                                 node_match=node_match, edge_match=edge_match)
        rec = {'same_labeled_graph': equal,
               'normalized_isomeric_smiles': [Chem.MolToSmiles(Chem.RemoveHs(m), isomericSmiles=True) for m in molecules]}
        if args.infer:
            a, b = [predict_standard_dGf_prime(m).reshape(-1) for m in molecules]
            delta = float(parts[0][0]) * a + float(parts[1][0]) * b
            rec.update(ensemble_size=len(delta), mean_kJ_mol=float(delta.mean()),
                       sd_kJ_mol=float(delta.std()), max_abs_member_kJ_mol=float(abs(delta).max()))
        records[rid] = rec
    source = args.model_root / 'dGbyG'
    files = [args.reactions, source / 'utils/NNFunc.py', source / 'utils/ChemFunc.py',
             source / 'Chemistry/Compound.py', source / 'network/GNNetwork.py']
    if args.infer:
        files += [source / 'api/utils.py', *sorted((source / 'network/model_params/MPNN_model_params').glob('*.pt'))]
    output = {'scope': __doc__, 'rdkit_version': rdBase.rdkitVersion, 'torch_version': torch.__version__,
              'input_sha256': {str(p): sha(p) for p in files}, 'results': records}
    args.output.write_text(json.dumps(output, indent=2, allow_nan=False) + '\n')
    print(json.dumps({'results': records}, indent=2))


if __name__ == '__main__':
    main()
