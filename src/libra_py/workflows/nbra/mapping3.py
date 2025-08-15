# *********************************************************************************
# * Copyright (C) 2025 Alexey V. Akimov
# *
# * This file is distributed under the terms of the GNU General Public License
# * as published by the Free Software Foundation, either version 3 of
# * the License, or (at your option) any later version.
# * See the file LICENSE in the root directory of this distribution
# * or <http://www.gnu.org/licenses/>.
# *
# *********************************************************************************/
"""

.. module:: mapping3
   :platform: Unix, Windows
   :synopsis: This module implements the methods to map single-electron properties (e.g. KS basis)
       to many-electron ones (e.g. Slater Determinat basis and SAC basis)

.. moduleauthor:: Alexey V. Akimov

"""

import os, sys, math
import numpy as np
from .mapping import sd2indx 

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform == "cygwin":
    from cyglibra_core import *
elif sys.platform == "linux" or sys.platform == "linux2":
    from liblibra_core import *


def ovlp_arb(_SD1, _SD2, S, active_space=None, verbose=False):
    """Compute the overlap of two generic SDs: <SD1|SD2>
    This functions translates the explicit use of two Python `for` loops
    into a more efficient way using numpy
    Args:

        _SD1 ( lists of ints ): first SD, such that:
            SeeAlso: ```inp``` in the ```sd2indx(inp)``` function

        _SD2 ( lists of ints ): second SD, such that:
            SeeAlso: ```inp``` in the ```sd2indx(inp)``` function

        S ( numpy array(N,N) ): is the matrix in the space of 1-el orbital

        active_space ( list of ints ): indices of the orbitals (starting from 1) to 
            include into consideration. If None - all the orbitals will be used [ default: None ]

        verbose ( Bool ): whether to print some extra info [ default: False] 

    Returns:
        complex: the overlap of the two determinants <SD1|SD2>

    """

    nbasis = S.shape[0] #num_of_rows

    SD1, SD2 = [], []

    if active_space == None:
        SD1, SD2 = _SD1, _SD2
    else:
        for a in _SD1:
            if abs(a) in active_space:
                SD1.append(a)
        for a in _SD2:
            if abs(a) in active_space:
                SD2.append(a)

    res = 0.0 + 0j
    if len(SD1) != len(SD2):
        print("Trying to compute an overlap of the SDs with different number of electrons")
        print("Exiting...")
        sys.exit(0)

        
    # Apply the determinant formula
    det_size = len(SD1)
    s = np.zeros( (det_size, det_size), dtype=np.float64); 
    """ Commented on 7/22/2025
    # ============================== The numpy version of the above double-for-loops
    # Find the tensor product of the SD1 and SD2
    SD_tensor_product = np.tensordot(SD1, SD2, axes=0)
    # Next, find the sign of these elementwise multiplications
    SD_tensor_product_sign = np.sign(SD_tensor_product)
    # Now, find where we have alpha-beta indices so that the 
    negative_sign_indices = np.where(SD_tensor_product_sign < 0)
    s[negative_sign_indices] = 0
    """
    # Let's build the matrix related to sd1 and sd2 from the KS orbitals
    # For this, we reuire to turn each element into matrix indices
    sd1 = sd2indx(SD1, nbasis, False, 0) # adopted from the `mapping` module
    sd2 = sd2indx(SD2, nbasis, False, 0)
    # What about beta indices?! We should add `nbasis/2` to them. This is added on 7/22/2025
    # With this simple approach, there is no need for tensor product or the `if else` clause brought
    # previously since the elements corresponding to the alpha-beta orbitals are already
    # zero in the two-spinor format of the KS matrices :)
    beta_indices = np.where(np.array(SD1) < 0)
    sd1 = np.array(sd1)
    sd1[beta_indices] += int(nbasis/2)

    beta_indices = np.where(np.array(SD2) < 0)
    sd2 = np.array(sd2)
    sd2[beta_indices] += int(nbasis/2)


    # Making the numpy array
    s = S[sd1,:][:,sd2]
    
    if verbose==True:
        print(s)
    res = np.linalg.det(s) #* phase

    return res



def ovlp_arb_2(_SD1, _SD2, S, active_space=None, verbose=False):
    """Compute the overlap of two generic SDs: <SD1|SD2>
    The difference of this function with the above one is that it uses explicitly
    two Pythonic `for` loops which are not optimized for very large sets

    Args:

        _SD1 ( lists of ints ): first SD, such that:
            SeeAlso: ```inp``` in the ```sd2indx(inp)``` function

        _SD2 ( lists of ints ): second SD, such that:
            SeeAlso: ```inp``` in the ```sd2indx(inp)``` function

        S ( numpy array(N,N) ): is the matrix in the space of 1-el orbital

        active_space ( list of ints ): indices of the orbitals (starting from 1) to 
            include into consideration. If None - all the orbitals will be used [ default: None ]

        verbose ( Bool ): whether to print some extra info [ default: False] 

    Returns:
        complex: the overlap of the two determinants <SD1|SD2>

    """

    nbasis = S.shape[0] #num_of_rows

    SD1, SD2 = [], []

    if active_space == None:
        SD1, SD2 = _SD1, _SD2
    else:
        for a in _SD1:
            if abs(a) in active_space:
                SD1.append(a)
        for a in _SD2:
            if abs(a) in active_space:
                SD2.append(a)

    res = 0.0 + 0j
    if len(SD1) != len(SD2):
        print("Trying to compute an overlap of the SDs with different number of electrons")
        print("Exiting...")
        sys.exit(0)

        
    # Apply the determinant formula
    det_size = len(SD1)
    s = np.zeros( (det_size, det_size), dtype=np.float64); 
    # ============================== The explicit version with for loops
    # Note: It seems that this explicit vesion doesn't consider the unrestricted spin calculations
    # Forming the overlap of the SDs
    for indx_I, I in enumerate(SD1):
        for indx_J, J in enumerate(SD2):
            # The overlap is non-zero only if the orbitals
            # are occupied with the same-spin electrons.
            if (I * J) > 0:
                i = abs(I) - 1
                j = abs(J) - 1
                s[indx_I, indx_J] = S[i, j] # .real #S.get(i, j).real
            else:
                s[indx_I, indx_J] = 0.0

    if verbose==True:
        print(s)
    res = np.linalg.det(s) #* phase

    return res




def ovlp_mat_arb(SD1, SD2, S, active_space=None):
    """Compute a matrix of overlaps in the SD basis

    Args:
        SD1 ( list of lists of N ints ): a list of N SD determinants, such that:
            SD1[iSD] is a list of integers defining which orbitals are
            occupied in SD with index ```iSD``` and how
            SeeAlso: ```inp``` in the ```sd2indx(inp)``` function

        SD2 ( list of lists of M ints ): a list of M SD determinants, such that:
            SD2[iSD] is a list of integers defining which orbitals are
            occupied in SD with index ```iSD``` and how
            SeeAlso: ```inp``` in the ```sd2indx(inp)``` function

        S ( CMATRIX(K,K) ): is the matrix in the full space of 1-el orbitals. Note - the mapped indices should not
            be larger than K-1

        active_space ( list of ints ): indices of the orbitals (starting from 1) to
            include into consideration. If None - all the orbitals will be used [ default: None ]

    Returns:
        CMATRIX(N, M): overlap matrix composed of elements <SD1(i)|SD2(j)>

    """

    N, M = len(SD1), len(SD2)
    res = np.zeros( (N, M), dtype=np.float64 ) #CMATRIX(N, M)

    for n in range(0, N):
        for m in range(0, M):
            res[n, m] = ovlp_arb(SD1[n], SD2[m], S, active_space, 0) 
            #res.set(n, m, val)
            #res[n, m] = val

    return res

