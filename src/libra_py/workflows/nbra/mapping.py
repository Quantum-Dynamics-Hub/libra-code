#*********************************************************************************
#* Copyright (C) 2018-2019 Brendan A. Smith, Wei Li, Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
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

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

#from libra_py import *


def sd2indx(inp,nbasis=0, do_sort=False, user_notation=0): 
    """

    This function maps a list of integers defining the 
    occupied spin-orbitals in a SD onto the global indices of the orbitals
    that are being occupied.     

    Args:
        inp ( list of ints ): indices of the occupied spin-orbitals. 
            Indexing starts with 1, not 0!

            Positive number ```n``` correspond to an alpha electron occupying the 
            orbital (whether it is an alpha-spatial or beta-spatial component) ```n```

            Negative number ```-n``` correspond to a beta electron occupying the 
            orbital (whether it is an alpha-spatial or beta-spatial component) ```n```

        nbasis ( int ): the total number of orbitals orbitals (both alpha- and beta-
            spatial components ) in the selected active space [currently not really used!] 

        do_sort ( Boolean ): the flag to tell whether the indices should be sorted in 
            a particular order:
                - True - for new scheme (needed for the SAC) [default]
                - use with False for Pyxaid mapping!

        user_notation (int): 
            - 0 : short-hand - in this case, the mapping goes into [0, N) range [default - because it has been used for a while]
            - 1 : extended - in this case, the mapping goes into [0, 2*N) range


    Returns: 
        list of ints: the indices of the orbitals occupied in this SD (different convention)

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

    sz = len(inp)

    spat = []
    for i in range(0,sz):
        # alpha
        if inp[i] > 0: 
            spat.append( inp[i] - 1 )

        # beta
        else:
            res = (abs(inp[i])) - 1
            #spat.append ( (abs(inp[i])) - 1  )  #+ nbasis/2 - 1
            if user_notation==1:
                res = res + int(nbasis/2)
            spat.append(res)

    # Rearrange in ascending order: this is needed for 
    # a consistency of the final results among different orbitals
    # Warning! But the reordering messes up the mapping procedure, so it 
    # is better not to have it. Note sure what effect it might have on spin-adaptation 
    # though
    out = list(spat)    
    if do_sort==True:
        out = sorted(spat)   

    return out


def reduce_determinants(_sd1, _sd2):
    """
    This function removes common parts of the two determinants
    
    Args:
        sd1, sd2 (lists of ints): determinants in the user representation
        
    Examples:
        >>  reduce_determinants( [1, -1], [2, -2])
        >>  ([1, -1], [2, -2])

        >>  reduce_determinants( [3,-3, 1, -1], [3, -3, 2, -2])
        >>  ([1, -1], [2, -2])

        >>  reduce_determinants( [1, -1, 2, -2, 3, -3], [1, -1, 2, -2, 3, -4])
        >>  ([-3], [-4])


    Returns:
        tuple of 2 lists: the reduced determinants
    
    """
    
    sd1, sd2 = [], []
    sz1 = len(_sd1)
    sz2 = len(_sd2)
    
    for i in range( sz1 ):
        if _sd1[i] not in _sd2:
            sd1.append( _sd1[i] )
                
    for i in range( sz2 ):
        if _sd2[i] not in _sd1:
            sd2.append( _sd2[i] )
            
    return sd1, sd2




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
    sd1 = sd2indx(SD1,nbasis, False, 0)
    sd2 = sd2indx(SD2,nbasis, False, 0)

    res = delta(Py2Cpp_int(sd1), Py2Cpp_int(sd2))

    x = 0.0+0.0j
    # The SDs are coupled
    if res[0]==1:
        x = X.get( res[1], res[2])


    sd1 = sorted(SD1)
    sd2 = sorted(SD2)
    res = delta(Py2Cpp_int(sd1), Py2Cpp_int(sd2))

    spin_ovlp = 1.0
    # Found 1 different entry
    if res[0]==0:
        spin_ovlp = 0.0
    elif res[0]==1:
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

    nbasis = e.num_of_rows

    sd = sd2indx(SD,nbasis)
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

    n = len(SD)
    E = CMATRIX(n,n)
  
    E0 = 0.0 #energy_arb(SD[0], e) + dE[0]*(1.0+0.0j)
    for i in range(0,n):
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
    act_sp2 = list(range(sz, 2*sz))
    
    S = CMATRIX(2*sz, 2*sz)

    push_submatrix(S, s, act_sp1, act_sp1)
    push_submatrix(S, zero, act_sp1, act_sp2)
    push_submatrix(S, zero, act_sp2, act_sp1)
    push_submatrix(S, s, act_sp2, act_sp2)

    return S


def num_of_perms(x):
    n = len(x)-1 # number of adjacent pairs of numbers
    cnt = 0
    for i in range(n):
        if x[i]>x[i+1]:
            cnt += 1
    return cnt



def ovlp_arb(SD1, SD2, S, use_minimal=False, user_notation=0):
    """Compute the overlap of two generic SDs: <SD1|SD2>

    Args:

        SD1 ( lists of ints ): first SD, such that:
            SeeAlso: ```inp``` in the ```sd2indx(inp,nbasis)``` function

        SD2 ( lists of ints ): second SD, such that:
            SeeAlso: ```inp``` in the ```sd2indx(inp,nbasis)``` function

        S ( CMATRIX(N,N) ): is the matrix in the space of 1-el spin-orbitals (either spin-diabatic or 
            spin-adiabatic or both) , N - is the number of 1-el orbitals.

        use_minimal ( Boolean ): If True, use the minimal subset of Kohn-Sham orbitals needed to describe 
            the overlap of the SDs      

    Returns:
        complex: the overlap of the two determinants <SD1|SD2>

    """

    nbasis = S.num_of_rows

    if use_minimal == True:
        SD1, SD2 = reduce_determinants(SD1, SD2)

    # Converting the SDs provided by the user into the internal format to be read by Libra
    sd1 = sd2indx(SD1,nbasis, False, user_notation)
    sd2 = sd2indx(SD2,nbasis, False, user_notation)

    # Compute the phase using the original determinants in the internal notation
    phase = (-1)**(num_of_perms(sd1) + num_of_perms(sd2))

    # Now reduce the determinants for faster calculations
    if use_minimal == True:
        SD1, SD2 = reduce_determinants(SD1, SD2)  

        # Convert the SDs to the internal notation again, but this time we'd be using the reduced ones
        sd1 = sd2indx(SD1,nbasis, False, user_notation)
        sd2 = sd2indx(SD2,nbasis, False, user_notation)


    res = 0.0+0j
    if len(sd1)>0 and len(sd2)>0:
      if len(sd1)==len(sd2):
    
          # Now apply the determinant formula
          s = CMATRIX(len(sd1),len(sd2))
          # Forming the overlap of the SDs
          for i in range(0,len(sd1)):
              for j in range(0,len(sd2)):
                  # The overlap is non-zero only if the orbitals
                  # are occupied with the same-spin electrons. 
                  if (SD1[i] * SD2[j]) > 0:          
                      s.set(i,j,S.get(sd1[i],sd2[j]))
                  else:
                      s.set(i,j,0.0,0.0)

          res = det(s) * phase
      else:
          # Checking if the matrix is square
          print("\nWARNING: the matrix of Kohn-Sham orbitial overlaps is not a square matrix") 
          print("\nExiting now ..")
          print("sd1 = ", sd1)
          print("sd2 = ", sd2)
          print("len(sd1) = ", len(sd1))
          print("len(sd2) = ", len(sd2))
          sys.exit(0)

    return res



def ovlp_arb_mo(SD1, SD2, S, user_notation=0):
    """Compute the overlap of two generic SDs: <SD1|SD2>

    This function relies on the methodology described in the Pyxaid paper 
    and later in the one by Izmaylov which showed that the NACs in the SD basis 
    can be directly mapped to the NACs in the orbital basis

    Args:

        SD1 ( lists of ints ): first SD, such that:
            SeeAlso: ```inp``` in the ```sd2indx(inp,nbasis)``` function

        SD2 ( lists of ints ): second SD, such that:
            SeeAlso: ```inp``` in the ```sd2indx(inp,nbasis)``` function

        S ( CMATRIX(N,N) ): is the matrix in the space of 1-el spin-orbitals (either spin-diabatic or 
            spin-adiabatic or both) , N - is the number of 1-el orbitals.

    Returns:
        complex: the overlap of the two determinants <SD1|SD2>

    """

    nbasis = S.num_of_rows
    # Converting the SDs provided by the user into the internal format to be read by Libra
    sd1 = sd2indx(SD1,nbasis, False, user_notation)
    sd2 = sd2indx(SD2,nbasis, False, user_notation)

    res = delta(Py2Cpp_int(sd1), Py2Cpp_int(sd2))

    phase = (-1)**(num_of_perms(sd1) + num_of_perms(sd2))

    s = 0.0
    # The SDs are coupled
    if res[0]==1:
        s = S.get( res[1], res[2])  


    # For similar SDs
    if sd1==sd2:
        #         TEST on June 16, 2022 - AVA
        # Forming the overlap of the SDs
        x = CMATRIX(len(sd1),len(sd2))
        for i in range(0,len(sd1)):
            for j in range(0,len(sd2)):
                if (SD1[i] * SD2[j]) > 0:
                    x.set(i,j,S.get(sd1[i],sd2[j]))
        
        s = det(x)
        
        #s = 1.0+0.0j 
    return s * phase



def ovlp_mat_arb(SD1, SD2, S, use_minimal=True, use_mo_approach=False, user_notation=0):
    """Compute a matrix of overlaps in the SD basis 

    Args:
        SD1 ( list of lists of N ints ): a list of N SD determinants, such that:
            SD1[iSD] is a list of integers defining which orbitals are 
            occupied in SD with index ```iSD``` and how 
            SeeAlso: ```inp``` in the ```sd2indx(inp,nbasis)``` function

        SD2 ( list of lists of M ints ): a list of M SD determinants, such that:
            SD2[iSD] is a list of integers defining which orbitals are 
            occupied in SD with index ```iSD``` and how 
            SeeAlso: ```inp``` in the ```sd2indx(inp,nbasis)``` function

        S ( CMATRIX(K,K) ): is the matrix in the space of 1-el spin-orbitals (either spin-diabatic or 
            spin-adiabatic or both) , K - is the number of 1-el orbitals.

        use_minimal ( Boolean ): If True, use the minimal subset of Kohn-Sham orbitals needed to describe 
            the overlap of the SDs. This is passed to the funciton that computes the overlaps [ default: True ]
 
        use_mo_approach ( Boolean ): If True, the computation of the time-overlaps of SDs will be based on the 
            computing the determinant of the matrix of the corresponding MO overlaps. This is a general
            approach, so by default, this value is set to False. In fact, it doesn't make sense for anything 
            but computing the NACs (eventually). However, as the size of the SDs increases (in terms of how many
            MOs are included), the complexity of this function increases cubically, so it may be very expensive.
            For computing the NACs via the HST formula, it is advisable (as shown in the Pyxaid paper and later
            by Izmaylov) to compute NACs between SDs in terms of the NACs of the MOs - this is WAAAY faster approach.
            But again, this is not a valid method if one needs just the overlaps of the SDs. [ default: False ]
            If `use_mo_approach==True`, the parameter `use_minimal` is not used (see the pun)
            


    Returns:
        CMATRIX(N, M): overlap matrix composed of elements <SD1(i)|SD2(j)>

    """

    N, M = len(SD1), len(SD2)
    res = CMATRIX(N,M)

    for n in range(0,N):
        for m in range(0,M):
            val = 0.0
            if use_mo_approach==True:
                val = ovlp_arb_mo(SD1[n], SD2[m], S, user_notation)
            else: 
                val = ovlp_arb(SD1[n], SD2[m], S, use_minimal, user_notation)

            res.set(n,m, val)

    return res



