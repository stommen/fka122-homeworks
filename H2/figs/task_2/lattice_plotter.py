import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

sub_A = np.loadtxt('data/task_2/lattice/sub_A_10.csv', delimiter=',')
sub_B = np.loadtxt('data/task_2/lattice/sub_B_10.csv', delimiter=',')
atom_positions = np.concatenate((sub_A, sub_B), axis=0)

neighbors = np.loadtxt('data/task_2/lattice/neighbors.csv', delimiter=',', dtype=int)
neighbors_A = neighbors[:1000, :]
neighbors_B = neighbors[1000:, :]

def plot_lattice_config(sub_A, sub_B, neighbors, idx=None, its_eq=0, its=0, T=0):
    fig, ax = plt.subplots(subplot_kw={'projection': '3d'})
    colors = np.array(['silver', 'chocolate'])
    if its_eq:
        atoms = np.loadtxt(f'data/task_2/lattice/atoms_{its_eq}_{its}_{T}.csv', delimiter=',', dtype=int)
        color = colors[atoms]
        ax.scatter(atom_positions[:, 0], atom_positions[:, 1], atom_positions[:, 2], c=color, marker='o', alpha=0.7, s=10)
        cu_A = np.sum(atoms[:1000])
        zn_A = 1000 - cu_A
        cu_B = np.sum(atoms[1000:])
        zn_B = 1000 - cu_B
        df = pd.DataFrame({'Cu_A': [cu_A], 'Zn_A': [zn_A], 'Cu_B': [cu_B], 'Zn_B': [zn_B], 'T': [T]})
        print(df)
    else:
        ax.scatter(sub_A[:, 0], sub_A[:, 1], sub_A[:, 2], c='chocolate', marker='o', alpha=0.7, s=10)
        ax.scatter(sub_B[:, 0], sub_B[:, 1], sub_B[:, 2], c='silver', marker='o', alpha=0.7, s=10)

    if idx is not None:
        atom_A = sub_A[idx]
        atom_neighbors_A = neighbors[idx, :]
        ax.scatter(atom_A[0], atom_A[1], atom_A[2], c='chocolate', marker='o')
        for neighbor_idx in atom_neighbors_A:
            neighbor = atom_positions[neighbor_idx, :]
            ax.scatter(neighbor[0], neighbor[1], neighbor[2], c='black', marker='o', alpha=1)

    plt.show()

plot_lattice_config(sub_A, sub_B, neighbors, its_eq=200000, its=1000000, T=400)
