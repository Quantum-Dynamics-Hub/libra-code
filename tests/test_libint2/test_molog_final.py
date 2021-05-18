import numpy as np
from liblibra_core import *
import CP2K_methods
from libra_py import data_conv


nprocs = 10
coord = CP2K_methods.extract_coordinates('cdse13-1.xyz',1)
basis_set_files_path = ['/home/97425008/cp2k-v7/cp2k/data/BASIS_MOLOPT']
unique_atoms = ['Cd','Se']
basis_set_names = ['DZVP-MOLOPT-SR-GTH','DZVP-MOLOPT-SR-GTH']
data = CP2K_methods.find_basis_set(basis_set_files_path,unique_atoms,basis_set_names)
shell_1 = CP2K_methods.make_shell(coord, data, True)
# print(nbasis(shell_1))
AO = compute_overlaps(shell_1,shell_1,nprocs)
AO = data_conv.MATRIX2nparray(AO)
ener, eig = CP2K_methods.read_molog_file('CdSe13-CdSe13-1_0.MOLog')
# l_values
l_vals = CP2K_methods.molog_lvals('CdSe13-CdSe13-1_0.MOLog')
# new indicies based on l_vals
new_indices = CP2K_methods.resort_molog_eigenvectors(l_vals)
eigenvectors = []
for i in range(len(eig[0])):
    # the new and sorted eigenvector
    eigenvector = eig[0][i]
    eigenvector = eigenvector[new_indices]
    # append it to the eigenvectors list
    eigenvectors.append(eigenvector)
eigs = np.array(eigenvectors)
S = np.linalg.multi_dot([eigs, AO ,eigs.T])
print(np.diag(S))

