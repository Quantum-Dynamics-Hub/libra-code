import numpy as np
import time
from liblibra_core import *
from libra_py import CP2K_methods
from libra_py import data_conv

# number of processors
nprocs = 4
# setting up a timer
t1 = time.time()
# extract the coordinates from the trajectory file with the time step
time_step = 1
coord = CP2K_methods.extract_coordinates('cdse13-1.xyz',time_step)
# the full path to different BASIS set files 
basis_set_files_path = ['/home/97425008/cp2k-v7/cp2k/data/BASIS_MOLOPT']
# the unique atoms present in the trajectory
unique_atoms = ['Cd','Se']
# the unique atoms respective basis set names
basis_set_names = ['DZVP-MOLOPT-SR-GTH','DZVP-MOLOPT-SR-GTH']
# the basis set data for the  unique atoms
data = CP2K_methods.find_basis_set(basis_set_files_path,unique_atoms,basis_set_names)
# create the libint2 shell based on the coordinates and the basis set data
# in spherical coordinates (the spherical flag is set to True)
shell_1 = CP2K_methods.make_shell(coord, data, True)
# You can print the number of basis in the shell to make sure
# it is the same as in the MOLog files
print('The number of atomic orbital basis set is:', nbasis(shell_1))
# compute the AO overlap matrix
AO = compute_overlaps(shell_1,shell_1,nprocs)
# turn it into a numpy array for easier dot product multiplication
AO = data_conv.MATRIX2nparray(AO)
# extract the eigenvectors and their energies
ener, eig = CP2K_methods.read_molog_file('CdSe13-CdSe13-1_0.MOLog')
# all the angular momentum values in the MOLog file 
l_vals = CP2K_methods.molog_lvals('CdSe13-CdSe13-1_0.MOLog')
# new indicies based on l_vals
new_indices = CP2K_methods.resort_molog_eigenvectors(l_vals)
# reindexing the eigenvectors
eigenvectors = []
for i in range(len(eig[0])):
    # the new and sorted eigenvector
    eigenvector = eig[0][i]
    eigenvector = eigenvector[new_indices]
    # append it to the eigenvectors list
    eigenvectors.append(eigenvector)
# make it a numoy array
eigs = np.array(eigenvectors)
# compute the MO overlap
S = np.linalg.multi_dot([eigs, AO ,eigs.T])
# print the diagonal elements of the MO overlap matrix
# to make sure that you get 1 on the diagonal
print(np.diag(S))
print('Elapsed time:',time.time()-t1)
