from functools import partial, reduce

import rdkit
from rdkit import Chem
from rdkit.Chem.rdchem import Atom, Bond, HybridizationType, ChiralType, BondType, BondStereo

import torch
from torch_geometric.data import Data
from torch_geometric.utils import add_self_loops


# Atom's features
Num_atomic_number = 119 # including the extra mask tokens
Num_atom_hybridization = len(HybridizationType.values)
Num_atom_aromaticity = 2 # Atom's aromaticity (not aromatic or aromactic)
Num_atom_chirality = len(ChiralType.values)
Num_atom_charge = 9
Num_atom_degree = 5

# Bond's features
Num_bond_type = len(BondType.values) # #including aromatic and self-loop edge(value = 22)
# probably necessary
Num_bond_atom_i = Num_atomic_number # soure atom
Num_bond_atom_j = Num_atomic_number # taget atom

atom_funs = {'atomic number':(Num_atomic_number, Atom.GetAtomicNum),
             'hybridization': (Num_atom_hybridization, Atom.GetHybridization),
             'aromaticity': (Num_atom_aromaticity, Atom.GetIsAromatic),
             'charge': (Num_atom_charge, Atom.GetFormalCharge),
             'chirality': (Num_atom_chirality, Atom.GetChiralTag),
             'degree': (Num_atom_degree, Atom.GetDegree),
             }

bond_funs = {'bond type':(Num_bond_type+1, lambda x:x.GetBondType()), # +1 for self loop
             'begin atom num':(Num_atomic_number, lambda x:x.GetBeginAtom().GetAtomicNum()),
             'end atom num':(Num_atomic_number, lambda x:x.GetEndAtom().GetAtomicNum()),
             }


def one_hot(num, idx):
    vector = [0]*num
    vector[idx] = 1
    return vector


def mol_to_graph_data(mol:rdkit.Chem.rdchem.Mol, 
                      atom_features=['atomic number', 'hybridization', 'aromaticity', 'charge'], 
                      bond_features=['bond type']) -> Data:
    # return (x, bond_index, bond_attr)
    # atoms features: including such below(ordered).
    atom_featurizers = [atom_funs[x] for x in atom_features]
    atoms_features = []

    for atom in mol.GetAtoms():
        feature = reduce(lambda x,y:x+y, [one_hot(num, fun(atom).real) for num, fun in atom_featurizers])
        atoms_features.append(feature)
    atoms_features = torch.tensor(atoms_features, dtype=torch.float32)

    
    # bonds index
    bonds_index = torch.empty(size=(2,0), dtype=torch.int64)
    
    # bonds attributes: including such below(ordered)
    # bond_type + bond_begin_atom_features + bond_end_atom_features
    bond_featurizers = [bond_funs[x] for x in bond_features]
    bond_dim = sum([num for num, _ in bond_featurizers])
    bonds_attrs = torch.empty(size=(0, bond_dim), dtype=torch.float32)
        
    for bond in mol.GetBonds():
        begin, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()

        bonds_index = torch.cat((bonds_index, torch.tensor([begin,end]).unsqueeze(1)), dim=1)
        attr = reduce(lambda x,y:x+y, [one_hot(num, fun(bond).real) for num, fun in bond_featurizers])
        bonds_attrs = torch.cat((bonds_attrs, torch.tensor(attr).unsqueeze(0)), dim=0)

        bonds_index = torch.cat((bonds_index, torch.tensor([end,begin]).unsqueeze(1)), dim=1)
        attr = reduce(lambda x,y:x+y, [one_hot(num, fun(bond).real) for num, fun in bond_featurizers])
        bonds_attrs = torch.cat((bonds_attrs, torch.tensor(attr).unsqueeze(0)), dim=0)

    #bonds_index, _ = add_self_loops(bonds_index, num_nodes=mol.GetNumAtoms())
    
    return Data(x=atoms_features, edge_index=bonds_index, edge_attr=bonds_attrs)



def old_mol_to_graph_data(mol:rdkit.Chem.rdchem.Mol) -> Data:
    # return (x, bond_index, bond_attr)
    # atoms features: including such below(ordered).
    atomic_number, hybridization, aromaticity, chirality, charge, degree = [], [], [], [], [], []
    features = [atomic_number, hybridization, aromaticity, chirality, charge, degree]
    atom_features = torch.empty(size=(0,len(features)), dtype=torch.int64)

    for atom in mol.GetAtoms():
        atomic_number = atom.GetAtomicNum()
        hybridization = atom.GetHybridization().real
        aromaticity = atom.GetIsAromatic().real
        chirality = atom.GetChiralTag().real
        charge = atom.GetFormalCharge() + 4
        degree = atom.GetDegree()

        features = torch.tensor([atomic_number, hybridization, aromaticity, chirality, charge, degree])
        atom_features = torch.cat((atom_features, features.unsqueeze(0)), dim=0)

    
    # bonds index
    bond_index = torch.empty(size=(2,0), dtype=torch.int64)
    
    # bonds attributes: including such below(ordered)
    # bond_type + bond_begin_atom_features + bond_end_atom_features
    bond_attr = torch.empty(size=(0, 2*len(features)+1), dtype=torch.int64)
        
    for bond in mol.GetBonds():
        start, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        bond_index = torch.cat((bond_index, torch.tensor([start,end]).view(2,1)), dim=1)
        bond_index = torch.cat((bond_index, torch.tensor([end,start]).view(2,1)), dim=1)

        bond_type = torch.tensor(bond.GetBondType().real, dtype=torch.long).view(1)
        x_i = atom_features[start,:]
        x_j = atom_features[end,:]
        bond_i_j = torch.cat((bond_type, x_i, x_j), dim=0)
        bond_j_i = torch.cat((bond_type, x_j, x_i), dim=0)

        bond_attr = torch.cat((bond_attr, bond_i_j.view(1, 2*len(features)+1)), dim=0)
        bond_attr = torch.cat((bond_attr, bond_j_i.view(1, 2*len(features)+1)), dim=0)
    
    return Data(x=atom_features, edge_index=bond_index, edge_attr=bond_attr)


