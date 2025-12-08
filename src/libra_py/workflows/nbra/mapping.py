# *********************************************************************************
# * Copyright (C) 2024-2025 Mohammad Shakiba, Alexey V. Akimov
# * Copyright (C) 2018-2024 Brendan A. Smith, Wei Li, Alexey V. Akimov
# *
# * This file is distributed under the terms of the GNU General Public License
# * as published by the Free Software Foundation, either version 3 of
# * the License, or (at your option) any later version.
# * See the file LICENSE in the root directory of this distribution
# * or <http://www.gnu.org/licenses/>.
# *
# *********************************************************************************/
"""

.. module:: mapping
   :platform: Unix, Windows
   :synopsis: This module implements the methods to map single-electron properties (e.g. KS basis)
       to many-electron ones (e.g. Slater Determinat basis)

.. moduleauthor:: Brendan Smith, Wei Li, and Alexey V. Akimov

"""

import os
import sys
import math
import numpy as np
# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform == "cygwin":
    from cyglibra_core import *
elif sys.platform == "linux" or sys.platform == "linux2":
    from liblibra_core import *

import libra_py.data_conv as data_conv
import libra_py.citools.slatdet as slatdet


def sd2indx(inp, nbasis=0, do_sort=False, user_notation=0):
    """
    This function maps a list of integers defining the occupied spin-orbitals in a SD onto 
    the global indices of the orbitals that are being occupied.

    Args:
        inp ( list of ints ): indices of the occupied spin-orbitals.
            Indexing starts with 1, not 0!

            Positive number ```n``` correspond to an alpha electron occupying the
            orbital (whether it is an alpha-spatial or beta-spatial component) ```n```

            Negative number ```-n``` correspond to a beta electron occupying the
            orbital (whether it is an alpha-spatial or beta-spatial component) ```n```

        nbasis ( int ): the total number of orbitals orbitals (both alpha- and beta-
            spatial components ) in the selected active space

        do_sort ( Boolean ): the flag to tell whether the indices should be sorted in
            a particular order (according to canonic ordering defined in the `citools.slatdet`):
                - True - for new scheme (needed for the SAC) [default]
                - use with False for Pyxaid mapping!

        user_notation (int):
            - 0 : short-hand - in this case, the mapping goes into [0, N) range [default - because it has been used for a while]
            - 1 : extended - in this case, the mapping goes into [0, 2*N) range


    Returns:
        - list of ints: the indices of the orbitals occupied in this SD (different convention)
        - integer: parity - reordering of the returned SD compared to the original ordering of spin-orbitals in the input
        - list of ints: the SD indices in the internal notation

    Example:

        Assume out active space consists of 4 orbitals:  |1^a>, |1^b>, |2^a>, |2^b>.
        Here, the superscript indicates the spin-component of the given orbital.

        In the step2 they are grouped as |1^a>, |2^a>, |1^b>, |2^b>

        And indexed (overall) as:           1     2      3      4  (as used in the input)
                                            0     1      2      3  (as used in the output)
        Note: this indexing convention assumes starting from 1, not from 0!

        So the ```nbasis``` should be = 4

        Input SD = [1, -3] means that the orbital |1^a> is occupied by an alpha electron |a>
        and the orbital |1^b> is occupied by a beta electron |b>, so the SD is:

        SD = 1/sqrt(2) * | |1^a>|a>   |1^b>|b> |

        This function will return the global indices of the occupied orbitals in the
        set of active orbitals. That is the conversion will be:
        [1, -3] -> [0, 2]
        Note that the output indices start from 0

    """

    inp = np.asarray(inp)

    spat = np.where(
        inp > 0,
        inp - 1,
        np.abs(inp) - 1 + (int(nbasis/2) if user_notation == 1 else 0)
    )

    # Rearrange according to a pre-defined canonical ordering, compute
    # the parity of the correponding permutation (compared to the input)
    out = list(spat)
    parity = 1
    if do_sort:
        out = sorted(spat, key=slatdet.canonical_sort_key)
        parity = slatdet.permutation_parity(spat, out)

    return np.array(out), parity



def ovlp_arb(_SD1, _SD2, S, active_space=None, do_sort=False, user_notation=0, verbose=False):
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

        do_sort ( Boolean ): the flag to tell whether the indices should be sorted in
            a particular order (according to canonic ordering defined in the `citools.slatdet`):
                - True - for new scheme (needed for the SAC) [default]
                - use with False for Pyxaid mapping!

        user_notation (int):
            - 0 : short-hand - in this case, the mapping goes into [0, N) range [default - because it has been used for a while]
            - 1 : extended - in this case, the mapping goes into [0, 2*N) range

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
    sd1,_,_ = sd2indx(SD1, nbasis, do_sort, user_notation)
    sd2,_,_ = sd2indx(SD2, nbasis, do_sort, user_notation)


    """ 
    ALEXEY: Instead of this trick, use user_notation = 1

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

    """

    # Making the numpy array
    s = S[sd1,:][:,sd2]
    
    if verbose==True:
        print(s)
    res = np.linalg.det(s)

    return res



def ovlp_mat_arb(SD1, SD2, _S, active_space=None, do_sort=False, user_notation=0, verbose=False):
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

        _S ( numpy array(N,N) ): is the matrix in the space of 1-el orbital. 
            Note - the mapped indices should not be larger than K-1

        active_space ( list of ints ): indices of the orbitals (starting from 1) to
            include into consideration. If None - all the orbitals will be used [ default: None ]

        do_sort ( Boolean ): the flag to tell whether the indices should be sorted in
            a particular order (according to canonic ordering defined in the `citools.slatdet`):
                - True - for new scheme (needed for the SAC) [default]
                - use with False for Pyxaid mapping!

        user_notation (int):
            - 0 : short-hand - in this case, the mapping goes into [0, N) range [default - because it has been used for a while]
            - 1 : extended - in this case, the mapping goes into [0, 2*N) range

        verbose ( Bool ): whether to print some extra info [ default: False]

    Returns:
        numpy.array(N, M): overlap matrix composed of elements <SD1(i)|SD2(j)>

    """

    # For transition compatibility
    S = None
    if type(_S)==CMATRIX or type(_S)==MATRIX:
        S = data_conv.MATRIX2nparray(_S, _dtype=np.complex128)
    else:
        S = _S

    N, M = len(SD1), len(SD2)
    res = np.zeros( (N, M), dtype=np.float64 )

    for n in range(0, N):
        for m in range(0, M):
            res[n, m] = ovlp_arb(SD1[n], SD2[m], S, active_space, do_sort, user_notation, verbose)

    return res



def map_gen_matrix(SD1, SD2, X):
    """
    Args:
        SD1 ( lists of ints ): first SD, such that:
            SeeAlso: ```inp``` in the ```sd2indx(inp,nbasis)``` function

        SD2 ( lists of ints ): second SD, such that:
            SeeAlso: ```inp``` in the ```sd2indx(inp,nbasis)``` function

        X ( CMATRIX(nmo, nmo)) : the amatrix in a single-particple picture

    Returns:
        complex: mapped SD-basis property for the two determinants SD1 and SD2

    """

    nbasis = X.num_of_rows
    # Converting the SDs provided by the user into the internal format to be read by Libra
    sd1 = sd2indx(SD1, nbasis, False, 0)
    sd2 = sd2indx(SD2, nbasis, False, 0)

    res = delta(Py2Cpp_int(sd1), Py2Cpp_int(sd2))

    x = 0.0 + 0.0j
    # The SDs are coupled
    if res[0] == 1:
        x = X.get(res[1], res[2])

    sd1 = sorted(SD1)
    sd2 = sorted(SD2)
    res = delta(Py2Cpp_int(sd1), Py2Cpp_int(sd2))

    spin_ovlp = 1.0
    # Found 1 different entry
    if res[0] == 0:
        spin_ovlp = 0.0
    elif res[0] == 1:
        if res[1] * res[2] > 0.0:
            spin_ovlp = 1.0
        else:
            spin_ovlp = 0.0

    return x * spin_ovlp


def energy_arb(SD, e):
    """Computes the energy of a Slater Determinat (SD) state

    Args:
        SD ( list of ints ): indices of the orbitals occupied in this SD
            SeeAlso: ```inp``` in the ```sd2indx(inp,nbasis)``` function

        e ( MATRIX(nbasis, nbasis) ): energies of all 1-electron levels (alpha- and beta- components)

    Returns:
        double: the sum of the orbital energies for the occupied spin-orbitals -
            This is the approximation of the SD energy

    """
    if isinstance(e, np.ndarray)==True:
        nbasis = e.shape[0]
        sd = sd2indx(SD, nbasis)
        res = 0.0
        for i in sd:
            res = res + e[i, i]


    elif isinstance(e, CMATRIX)==True or isinstance(e, MATRIX)==True:
        nbasis = e.num_of_rows
        sd = sd2indx(SD, nbasis)
        res = 0.0+0.0j  
        for i in sd:
            res = res + e.get(i,i)

    return res


def energy_mat_arb(SD, e, dE):
    """Computes a matrix of the SD energies

    Args:
        SD (list of lists of ints ): a list of SD determinants, such that:
            SD[iSD] is a list of integers defining which orbitals are
            occupied in SD with index ```iSD``` and how
            SeeAlso: ```inp``` in the ```sd2indx(inp,nbasis)``` function

        e ( MATRIX(nbasis, nbasis) ): energies of all 1-electron levels (alpha- and beta- components)
        dE ( list of doubles ): energy corrections added to each SD

    Returns:
        CMATRIX(N,N): the matrix of energies in the SD basis. Here, N = len(SD) - the number of SDs.

    """
    E = None

    if isinstance(e, np.ndarray)==True:
        n = len(SD)
        E = np.zeros((n, n))
        E0 = energy_arb(SD[0], e) + dE[0]
        for i in range(0, n):
            E[i, i] = energy_arb(SD[i], e) + dE[i] - E0

    elif isinstance(e, CMATRIX)==True or isinstance(e, MATRIX)==True:
        n = len(SD)
        E = CMATRIX(n,n)
        E0 = energy_arb(SD[0], e) + dE[0]*(1.0+0.0j)
        for i in range(0, n):
            E.set(i,i, energy_arb(SD[i], e) + dE[i]*(1.0+0.0j) - E0 )

    return E


def orbs2spinorbs(s):
    """
    This function converts the matrices in the orbital basis (e.g. old PYXAID style)
    to the spin-orbital basis.
    Essentially, it makes a block matrix of a double dimension:
           ( s  0 )
    s -->  ( 0  s )

    This is meant to be used for backward compatibility with PYXIAD-generated data

    Args:
        s ( CMATRIX(N,N) ): whatever matrix in the orbital basis, N - is the number of
            doubly-occupied orbitals

    Returns:
        CMATRIX(2N, 2N): whatever matrix in the spin-orbital basis. It has a doubled dimensionality

    """

    sz = s.num_of_cols
    zero = CMATRIX(sz, sz)
    act_sp1 = list(range(0, sz))
    act_sp2 = list(range(sz, 2 * sz))

    S = CMATRIX(2 * sz, 2 * sz)

    push_submatrix(S, s, act_sp1, act_sp1)
    push_submatrix(S, zero, act_sp1, act_sp2)
    push_submatrix(S, zero, act_sp2, act_sp1)
    push_submatrix(S, s, act_sp2, act_sp2)

    return S


