import numpy as np
from liblibra_core import *
import CP2K_methods
import molden_methods
from libra_py import data_conv

nprocs = 10
# t1 = time.time()
shells_1, l_vals = molden_methods.molden_file_to_libint_shell('CdSe13-CdSe13-1_0.molden',True)
eig1, ener1 = molden_methods.eigenvectors_molden('CdSe13-CdSe13-1_0.molden',nbasis(shells_1),l_vals)
AO = compute_overlaps(shells_1,shells_1,nprocs)
AO = data_conv.MATRIX2nparray(AO)
new_indices = CP2K_methods.resort_molog_eigenvectors(l_vals)
eigenvectors = []
for i in range(len(eig1)):
    # the new and sorted eigenvector
    eigenvector = eig1[i]
    eigenvector = eigenvector[new_indices]
    # append it to the eigenvectors list
    eigenvectors.append(eigenvector)

    
eigs = np.array(eigenvectors)
S = np.linalg.multi_dot([eigs, AO ,eigs.T])
print(np.diag(S))

