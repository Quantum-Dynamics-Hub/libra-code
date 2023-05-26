#*********************************************************************************
#* Copyright (C) 2018-2019 Brendan A. Smith, Wei Li, Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 3 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
"""

.. module:: mapping2
   :platform: Unix, Windows
   :synopsis: This module implements the methods to map single-electron properties (e.g. KS basis)
       to many-electron ones (e.g. Slater Determinat basis). Revised version to simplify the user's experience

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


def sd2indx(sd): 
    """

    This function maps integers from the user that appear in the SD definition to internal indices

    Args:
        sd ( list of ints ): ids of the occupied orbitals in a given SD configuration. Indexing starts with 1, not 0!

            Positive numbers correspond to an alpha electron occupying the orbital

            Negative numbers correspond to a beta electron occupying the orbital

    Returns: 
        3 list of ints: the indices of the orbitals occupied in this SD in internal notation, the indexing starts from 0
             The first list - is for all orbitals, second - for alpha, and third - for beta

    Example: 

        Assume our active space consists of 2 orbitals: |1>,|2>

        One can then define 2-electron configurations like:  | 1a, 1b | or |1a, 2b|, etc.

        Here a and b refer to spin of electron occupying these orbitals

        On the user side, these configurations would be encoded as [1, -1] and [1, -2], respectively

        This function will map these numbers into [1, -1] -> [0, 0] and [1, -2] -> [0, 1]

    """
    sz = len(sd)

    res, resa, resb = [], [], []
    for i in range(0,sz):
        res.append( abs( sd[i]) - 1 )

        if sd[i]>0:
            resa.append(  abs( sd[i]) - 1 )
        if sd[i]<0:
            resb.append(  abs( sd[i]) - 1 )

    return res, resa, resb


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


def num_of_perms(x):
    """
    This function computes the number of incorrectly ordered pairs in a sequence of numbers `x`

    Args:
        x ( list ): a sequence of numeric numbers

    Returns:
        int: the number of incorrectly ordered adjacent pairs of numbers

    """
    n = len(x)-1 # number of adjacent pairs of numbers
    cnt = 0
    for i in range(n):
        if x[i]>x[i+1]:
            cnt += 1
    return cnt



def ovlp_arb(SD1, SD2, S, reduce_det=False):
    """Compute the overlap of two generic SDs: <SD1|SD2>

    Args:

        SD1 ( lists of ints ): first SD, such that:
            SeeAlso: ```inp``` in the ```sd2indx(inp)``` function

        SD2 ( lists of ints ): second SD, such that:
            SeeAlso: ```inp``` in the ```sd2indx(inp)``` function

        S ( CMATRIX(N,N) ): is the matrix in the space of 1-el orbital

        reduce_det ( Boolean ): If True, use the minimal distinct subset of orbitals needed to describe 
            the overlap of the SDs

    Returns:
        complex: the overlap of the two determinants <SD1|SD2>

    """

    nbasis = S.num_of_rows

    # Converting the SDs provided by the user into the internal format to be read by Libra
    sd1, sd1_a, sd1_b = sd2indx(SD1)
    sd2, sd2_a, sd2_b = sd2indx(SD2)

    # Compute the phase using the original determinants in the internal notation
    phase = (-1)**(num_of_perms(sd1) + num_of_perms(sd2))

    # Now reduce the determinants for faster calculations
    if reduce_det == True:
        SD1, SD2 = reduce_determinants(SD1, SD2)  

        # Convert the SDs to the internal notation again, but this time we'd be using the reduced ones
        sd1, _, _ = sd2indx(SD1)
        sd2, _, _ = sd2indx(SD2)

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


def ovlp_mat_arb(SD1, SD2, S, reduce_det=False):
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

        reduce_det ( Boolean ): If True, use the minimal subset of orbitals needed to describe 
            the overlap of the SDs. This is passed to the funciton that computes the overlaps [ default: True ]

    Returns:
        CMATRIX(N, M): overlap matrix composed of elements <SD1(i)|SD2(j)>

    """

    N, M = len(SD1), len(SD2)
    res = CMATRIX(N,M)
    N_orb_max = S.num_of_cols

    for i in range(N):
        if max(sd2indx(SD1[i])[0] ) > N_orb_max or min( sd2indx(SD1[i])[0] )<0:
            print(F"Error: Configurations {SD1[i]} is not consistent with the 1-electron matrix of dimensions {N_orb_max} x {N_orbs_max}\n")
            print("Exiting...\n")
            sys.exit(0)

    for i in range(M):
        if max( sd2indx(SD2[i])[0] ) > N_orb_max or min(sd2indx(SD2[i])[0] ) <0:
            print(F"Error: Configurations {SD2[i]} is not consistent with the 1-electron matrix of dimensions {N_orb_max} x {N_orbs_max}\n")
            print("Exiting...\n")
            sys.exit(0)



    for n in range(0,N):
        for m in range(0,M):
            val = ovlp_arb(SD1[n], SD2[m], S, reduce_det)
            res.set(n,m, val)

    return res



