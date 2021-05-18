import numpy as np
import time
from liblibra_core import *
from libra_py import CP2K_methods
from libra_py import molden_methods
from libra_py import data_conv

# number of processors
nprocs = 4
# set up the timer
t1 = time.time()
# creating the shell for the molden file
shells_1, l_vals = molden_methods.molden_file_to_libint_shell('CdSe13-CdSe13-1_0.molden',True)
# extracting the eigenvectors and energies from molden file
eig1, ener1 = molden_methods.eigenvectors_molden('CdSe13-CdSe13-1_0.molden',nbasis(shells_1),l_vals)
# compute the AO overlap matrix
AO = compute_overlaps(shells_1,shells_1,nprocs)
# turn it into a numpy array
AO = data_conv.MATRIX2nparray(AO)
# new indices for the the MOLog eigenvectors 
new_indices = CP2K_methods.resort_molog_eigenvectors(l_vals)
# making all the reindexed eigenvectors
eigenvectors = []
for i in range(len(eig1)):
    # the new and sorted eigenvector
    eigenvector = eig1[i]
    eigenvector = eigenvector[new_indices]
    # append it to the eigenvectors list
    eigenvectors.append(eigenvector)
# make it a numpy array to be able to work with    
eigs = np.array(eigenvectors)
# compute the MO overlap
S = np.linalg.multi_dot([eigs, AO ,eigs.T])
# print out the diagonal element of the MO overlap matrix
# to make sure you get 1 on the diagonal
print(np.diag(S))
print('Elapsed time:',time.time()-t1)
