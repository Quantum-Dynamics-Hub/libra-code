#*********************************************************************************
#* Copyright (C) 2018-2019 Brendan A. Smith, Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/

import os
import sys
import math

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *

def energy_arb(SD, e):
    # Computes the energy of the SD
    # SD = [int_list] the list of integers represents the occupied orbitals
    # e - is a diagonal matrix with 1-electron orbital energies

    nbasis = e.num_of_rows

    sd = sd2indx(SD,nbasis)
    res = 0.0+0.0j

    for i in sd:
        res = res + e.get(i,i) 
   
    return res 


def energy_mat_arb(SD, e, dE):
    # Computes a matrix of the SD energies 
    # SD - [ [list_of_ints], [list_of_ints], ... ]    # a list of SDs 
    # e - is a diagonal matrix with 1-electron orbital energies
    # dE - a list of energy corrections added to each SD

    n = len(SD)
    E = CMATRIX(n,n)

    E0 = energy_arb(SD[0], e) + dE[0]*(1.0+0.0j)
    for i in xrange(n):
        E.set(i,i, energy_arb(SD[i], e) + dE[i]*(1.0+0.0j) - E0 )

    return E

def sd2indx(inp,nbasis):
    """
    This function maps a list of integers defining a given spin
    configuration on the list of indices of the corresponding 
    spin-orbitals
    """

    sz = len(inp)
    out = [0] * sz

    for i in xrange(sz):
        # alpha
        if inp[i] > 0: 
            out[i] = inp[i] - 1
        # beta
        else:
            out[i] = (abs(inp[i])) + nbasis/2 - 1

    # Rearrange in ascending order
    out = sorted(out)

    return out


def ovlp_arb(SD1, SD2, S):
    """
    Overlap of two generic SDs: <SD1|SD2>
    SD1, SD2 - are the lists of indices of the spin-orbitals (in the corresponding
    active spaces) that are included in the two configurations
    S - is the matrix in the spaces of 1-el spin-orbitals (either spin-diabatic or 
    spin-adiabatic or both) 
    """

    nbasis = S.num_of_rows

    sd1 = sd2indx(SD1,nbasis)
    sd2 = sd2indx(SD2,nbasis)
    s = CMATRIX(len(sd1),len(sd2))

    #print "\nNOW WE ARE COMPUTING THE OVERLAP"
    #print SD1
    #print SD2
    #print sd1
    #print sd2

    for i in xrange(len(sd1)):
        for j in xrange(len(sd2)):
             s.set(i,j,S.get(sd1[i],sd2[j]))

    #s.show_matrix()
    #print "\n"
    #print det(s)

    return det(s) # Eq. 16


def ovlp_mat_arb(SD1, SD2, S):
    """
    A matrix composed of the SD overlaps
    SD1, SD2 - lists of lists, representing two bases of Slater Determinant functions
    e.g. SD1 = [[1,2], [1,3], [2,3]]
    S - is the elementary overlap of the 1-electron functions
    """

    N, M = len(SD1), len(SD2)
    res = CMATRIX(N,M)

    #print S.show_matrix()
    #print SD1
    #print SD2
    #print N
    #print M

    for n in xrange(N):
        for m in xrange(M):
            res.set(n,m, ovlp_arb(SD1[n], SD2[m], S))

    return res



