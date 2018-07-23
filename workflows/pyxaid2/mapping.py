#*********************************************************************************
#* Copyright (C) 2017 Alexey V. Akimov
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


def elementary_overlap(psi1, psi2):
    """
    This function creates a matrix of elementary overlaps - the overlaps
    of the 1-electron 2-component spin-orbitals. The spin-orbitals can have an
    arbitrary meaning: either spin-adiabatic (then both alpha and beta spatial
    parts are non-zero) or spin-diabatic (then one of the components is zero)
     
    psi1, psi2 - are the lists of 2 elements, each element is CMATRIX(npw, ndim) with the
    plane-wave coefficients for the spatial part for all spin-orbitals
    ndim may be different for psi1 and psi2

    psi1[0], psi2[0] - alpha
    psi1[1], psi2[1] - beta

    """
    n, m = psi1[0].num_of_cols, psi2[0].num_of_cols

    ovlp = CMATRIX(n,m)   # <psi1 | psi2 >

    ovlp = psi1[0].H() * psi2[0] +  psi1[1].H() * psi2[1]   # Eq. 17

    return ovlp    
    
    
def dia2indx(inp):
# This function maps a list of integers defining a diabatic spin-orbital configuration
# on the list of indices of the corresponding spin-orbitals
# Example: diabatic spin-orbital
# [1, -1, 2, -2, 3, -3, etc.] -> [0, 0, 1, 1, 2, 2, etc.] -> [ [0,1,2] , [0,1,2] ]
#                                                                alp       bet
# Change becuase the index of electron cannot be greater than the number of bands! 
# Ex) if 4 bands total: let's say user defines [1,-1,2,-3] = S1 
# dia2indx_old would give [0,1,2,5] - we cannot have the index 5 if only 4 bands!, due to how 
# SD overlaps are calculated herein. 

    sz = len(inp)

    alp = []
    bet = []
    for i in inp:
        if i > 0:
            alp.append(abs(i)-1)
        else:
            bet.append(abs(i)-1)       

    alp = sorted(alp)
    bet = sorted(bet)             

    return [alp,bet]


def ovlp_dd(SD1, SD2, S):
# Overlap of two diabatic SDs: <SD1|SD2>
# SD1, SD2 - are the lists of indices of the spin-orbitals (in the corresponding
# active spaces) that are included in the two configurations
# S - is the matrix in the spaces of 1-el spin-orbitals (either spin-diabatic or 
# spin-adiabatic or both) 
# See Eq. 16 in the write up

    # new formulation - sd1 and sd2 now a list of 2 lists
    sd1 = dia2indx(SD1)
    sd2 = dia2indx(SD2)

    print "sd1 = ", sd1
    print "sd2 = ", sd2

    # now, need to make submatrices
    s_alp = CMATRIX(len(sd1[0]),len(sd2[0]))
    pop_submatrix(S, s_alp, sd1[0], sd2[0])
    s_bet = CMATRIX(len(sd1[1]),len(sd2[1]))
    pop_submatrix(S, s_bet, sd1[1], sd2[1])

    print "\nJust computed (popped) s_alp and s_bet seperately in ovlp_dd\n"
    print "s_alp = ", s_alp.show_matrix()
    print "s_bet = ", s_bet.show_matrix()

    # Now, we are going to combine the two matrices via a direct sum

    s = CMATRIX(len(sd1[0])+len(sd1[1]),len(sd2[0])+len(sd2[1]))

    push_submatrix(s, s_alp, range(0,len(sd1[0])), range(0,len(sd2[0])))
    push_submatrix(s, s_bet, range(len(sd1[0]),len(sd1[0])+len(sd1[1])), range((len(sd2[0])),len(sd2[0])+len(sd2[1])))

    print "\nShowing s matrix\n"
    s.show_matrix()
    print "det(s) = ", det(s)
    return det(s)

def ovlp_mat_dd(SD1, SD2, S):
    #  A matrix composed of the SD overlaps
    #  SD1, SD2 - lists of lists, representing two bases of Slater Determinant functions
    #  e.g. SD1 = [[1,2], [1,3], [2,3]]
    #  S - is the elementary overlap of the 1-electron functions

    N, M = len(SD1), len(SD2)
    res = CMATRIX(N,M)

    for n in xrange(N):
        for m in xrange(M):
            res.set(n,m, ovlp_dd(SD1[n], SD2[m], S))

    print "\nJust finished calculing ovlp_mat_dd - printing res"
    res.show_matrix()
    print "\n"

    return res



def energy_dd(sd, e):
    # Computes the energy of an arbitrary SD
    # SD = [int_list] the list of integers represents the occupied orbitals
    # e - is a diagonal matrix with 1-electron orbital energies
   
    n = len(sd)
    res = 0.0+0.0j

    for i in xrange(n):
        for j in xrange(len(sd[i])):

            print "i=", i, "j=", j, "energy=", e.get(sd[i][j],sd[i][j])
            res = res + e.get(sd[i][j],sd[i][j])
    return res


def energy_mat_dd(SD, e, dE):
    # Computes a matrix of the SD energies 
    # SD - [ [list_of_ints], [list_of_ints], ... ]    # a list of SDs 
    # e - is a diagonal matrix with SD-electron orbital energies
    # dE - a list of energy corrections added to each SD

    n = len(SD)
    E = CMATRIX(n,n)

    E0 = 0.0 #energy_dd(SD[0], e) + dE[0]*(1.0+0.0j)
    for i in xrange(n):

        E.set(i,i, energy_dd(dia2indx(SD[i]), e) + dE[i]*(1.0+0.0j) - E0)

    return E

#========== Diabatic Above ==========#




#========== Adiabatic Below =========#

def adi2indx(inp,nks):
# This mapping is almost 1-to-1, just change the indexing convention
# inp - start indexing from 1
# out - start indexing from 0

    sz = len(inp)
    out = [0] * sz

    for i in xrange(sz):
        out[i] = inp[i] - 1
    # Add beta index basd on nks
    for i in xrange(sz):
        out.append( -out[i] - nks )    

    # Now seperate into alp and bet sublists
    alp = []
    bet = []
    for i in out:
        if i >= 0:
            alp.append(i)
        else:
            bet.append(abs(i))

    alp = sorted(alp)
    bet = sorted(bet)

    return [alp,bet]

    return out




def ovlp_aa(SD1, SD2, S):
    # Overlap of two spin-adiabatic SDs
    # SD1, SD2 - list of ints, representing SDs, e.g. [1,2]
    # S - is the elementary overlap of the 1-electron functions

    nks = S.num_of_cols / 2

    sd1 = adi2indx(SD1,nks)
    sd2 = adi2indx(SD2,nks)

    print "sd1", sd1
    print "sd2", sd2

    # now, need to make submatrices
    s_alp = CMATRIX(len(sd1[0]),len(sd2[0]))
    pop_submatrix(S, s_alp, sd1[0], sd2[0])
    s_bet = CMATRIX(len(sd1[1]),len(sd2[1]))
    pop_submatrix(S, s_bet, sd1[1], sd2[1])

    print "\nJust computed (popped) s_alp and s_bet seperately in ovlp_aa\n"
    print "s_alp = ", s_alp.show_matrix()
    print "s_bet = ", s_bet.show_matrix()

    s = CMATRIX(len(SD1),len(SD2))
    alp_row = s_alp.num_of_rows; alp_col = s_alp.num_of_cols
    bet_row = s_bet.num_of_rows; bet_col = s_bet.num_of_cols


    s = s_alp + s_bet

    print "\nShowing s matrix\n"
    s.show_matrix()
    print "det(s) = ", det(s)

    return det(s)


def ovlp_mat_aa(SD1, SD2, S):
    #  A matrix composed of the SD overlaps 
    #  SD1, SD2 - lists of lists, representing two bases of Slater Determinant functions
    #  e.g. SD1 = [[1,2], [1,3], [2,3]]
    #  S - is the elementary overlap of the 1-electron functions

    N, M = len(SD1), len(SD2)
    res = CMATRIX(N,M)

    for n in xrange(N):
        for m in xrange(M):
            res.set(n,m, ovlp_aa(SD1[n], SD2[m], S))

    print "\nJust finished calculing ovlp_mat_aa - printing res"
    res.show_matrix()
    print "\n"

    return res


      
def energy_aa(SD, e, nks):
    # Computes the energy of the SD
    # SD = [int_list] the list of integers represents the occupied orbitals
    # e - is a diagonal matrix with 1-electron orbital energies

    
    sd = adi2indx(SD, nks)
    res = 0.0+0.0j

    print "SD=", SD
    print "sd=", sd

    for i in sd[0]:

        print "sd[0]=", sd[0], "i=", i, "energy=", e.get(i,i)

        res = res + e.get(i,i) 
   
    return res 
   


def energy_mat_aa(SD, e, dE):
    # Computes a matrix of the SD energies 
    # SD - [ [list_of_ints], [list_of_ints], ... ]    # a list of SDs 
    # e - is a diagonal matrix with 1-electron orbital energies
    # dE - a list of energy corrections added to each SD

    nks = e.num_of_cols / 2
    n = len(SD)
    
    print "\nComputing H_el\n"
    E = CMATRIX(n,n)
    E0 = 0.0 #energy_aa(SD[0], e, nks) + dE[0]*(1.0+0.0j)
    for i in xrange(n):

        print "n=", i
        print "SD=", SD
        print "SD["+str(i)+"]=", SD[i]

        E.set(i,i, energy_aa(SD[i], e, nks) + dE[i]*(1.0+0.0j) - E0)  

    return E
              
