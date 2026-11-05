# write by Fan Wenchao
import time
from typing import Any, Dict
import rdkit
import numpy as np
import pandas as pd
from functools import lru_cache

import torch
import torch.nn as nn
from torch_geometric.loader import DataLoader

from dGbyG.config import inference_model_state_dict_dir, train_data_path
from dGbyG.utils.constants import *
from dGbyG.utils.ChemFunc import *
from dGbyG.utils.NNFunc import mol_to_graph_data
from dGbyG.network.GNNetwork import MP_network
from dGbyG.network.data import TrainDataset


TrainingData_df = pd.read_csv(train_data_path)
equation = TrainingData_df.loc[:, 'reaction']
standard_dG_prime = TrainingData_df.loc[:, 'standard_dg_prime']
TrainSet = TrainDataset(equations=equation, dGs=standard_dG_prime)

atom = torch.zeros(size=(139,), dtype=torch.float64)
bond = torch.zeros(size=(23,), dtype=torch.float64)
for data in TrainSet:
    atom += data.x.sum(dim=0)
    bond += data.edge_attr.sum(dim=0)
untrained_atom_attr = (atom == 0)
untrained_bond_attr = (bond == 0)


# GNN models construction
class networks(nn.Module):
    def __init__(self, dir) -> None:
        super().__init__()
        self.nets = nn.ModuleList([])
        for file in os.listdir(dir):
            path = os.path.join(dir, file)
            net = MP_network(atom_dim=TrainSet[0].x.size(1), bond_dim=TrainSet[0].edge_attr.size(1), emb_dim=300, num_layer=2)
            net.load_state_dict(torch.load(path, map_location=torch.device('cpu')))
            self.nets.append(net)
        self.num = len(self.nets)
    
    def forward(self, data, mode='molecule mode'):
        if mode == 'molecule mode':
            outputs = torch.zeros(size=(self.num, 1)) # shape=[number of net, 1]
        elif mode == 'atom mode':
            outputs = torch.zeros(size=(self.num, data.x.shape[0], 1)) # shape=[number of net, atom number, 1]

        for i, net in enumerate(self.nets):
            outputs[i] = net(data, mode) # net.shape = [1,1] or [atom number, 1]

        return outputs#.squeeze()

# read the model parameters
network = networks(inference_model_state_dict_dir)
network.eval()
device = 'cuda:0' if torch.cuda.is_available() else 'cpu'
network.to(device)



def predict_standard_dGf_prime(mol:rdkit.Chem.rdchem.Mol, mode='molecule mode') -> np.ndarray:
    # this function return predictions of a set of models with different baseline, 
    # so the results' mean and std cannot represent the real value of the mol
    data = mol_to_graph_data(mol).to(device)
    with torch.no_grad():
        standard_dGf_prime = network(data, mode=mode).cpu().numpy()#.reshape(-1)

    return standard_dGf_prime # shape=[number of nets, 1]
    


def predict_standard_dGr_prime(rxn_dict:Dict[rdkit.Chem.rdchem.Mol, int | float]) -> Tuple[float]:
    # 
    dGf_mat = []
    for mol, coeff in rxn_dict.items():
        dGr = predict_standard_dGf_prime(mol, mode='molecule mode') # shape=[number of nets, 1]
        dGf_mat.append(dGr * coeff)
    dGf_mat = np.asarray(dGf_mat) # shape=[number of mol, number of nets, 1]
    dGf_mat = dGf_mat.sum(axis=0).squeeze() # shape=[number of nets]

    return dGf_mat.mean(axis=0), dGf_mat.std(axis=0)



def predict_standard_dGr_prime_from_S(S:np.ndarray, mols:List) -> Tuple:
    # S.shape=[number of mols, number of reactions]
    mols_dGf = np.array([predict_standard_dGf_prime(mol) for mol in mols]) # shape=[number of mols, number of nets, 1]
    mols_dGf = mols_dGf.transpose(1,0,2) # shape=[number of nets, number of mols, 1]

    dGr = np.matmul(S.T, mols_dGf).squeeze() # shape=[number of nets, number of reactions]

    return dGr.mean(axis=0), dGr.std(axis=0)


def predict_standard_dGr_prime_from_rxn(rxn:rdkit.Chem.rdChemReactions.ChemicalReaction, radius:int=2):
    # 
    rxn.Initialize() # 
    if not rxn.IsInitialized():
        return np.nan

    reacting_atoms = reacting_map_num_of_rxn(rxn)
    if not reacting_atoms:
        return 0
    
    coeffs = [-1]*rxn.GetNumReactantTemplates() + [1] * rxn.GetNumProductTemplates()
    molecules = list(rxn.GetReactants()) + list(rxn.GetProducts())
    
    molecules_energy = []
    for mol, coeff in zip(molecules, coeffs): # the map numbers of changed atoms
        try:
            mol.UpdatePropertyCache()
        except:
            return np.nan
        mol = Chem.AddHs(mol)
        energy_changed_atom = map_num_to_idx_with_radius(mol, map_num=reacting_atoms, radius=radius)
        if not energy_changed_atom:
            changed_atoms_energy = 0
        else:
            atoms_energy = predict_standard_dGf_prime(mol, mode='atom mode') # shape=[number of nets, atom number, 1]
            changed_atoms_energy = atoms_energy[:,np.array(energy_changed_atom),:].sum(axis=1) # shape=[number of nets, 1]
        molecules_energy.append(changed_atoms_energy * coeff)
    print(molecules_energy)
    molecules_energy = np.asarray(molecules_energy) # shape=[number of molecules, number of nets, 1]
    molecules_energy = molecules_energy.sum(axis=0).squeeze() # shape=[number of nets]

    return molecules_energy.mean(), molecules_energy.std()



'''def predict_standard_dGf_prime(mol:rdkit.Chem.rdchem.Mol, mode='molecule mode') -> np.ndarray:
    graph_data = mol_to_graph_data(mol)
    for data in DataLoader([graph_data], batch_size=1):
        data.to(device)

    with torch.no_grad():
        standard_dGf_prime = network(data, out=mode).cpu().numpy().reshape(-1)

    atom_attr_distribution = graph_data.x.sum(dim=0)
    bond_attr_distribution = graph_data.edge_attr.sum(dim=0)
    if (atom_attr_distribution[untrained_atom_attr]).any() or (bond_attr_distribution[untrained_bond_attr]).any():
        return standard_dGf_prime + np.nan

    return standard_dGf_prime'''
