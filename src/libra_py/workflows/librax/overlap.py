#*********************************************************************************
#* Copyright (C) 2016 Kosuke Sato, Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/

## \file overlap.py
# This module implements the functions that calculates
# the overlap matrixes of atomic and molecular orbitals with different time steps.
# This returns the overlap matrix of molecular orbitals like  <MO(t)|MO(t+dt)>.


import os
import sys
import math
import numpy


if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

def AO_overlap(ao_i, ao_j):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] ao_i, ao_j : atomic orbital basis at different time step.
    # This function returns overlap matrix of atomic orbitals with different time step
    # like <AO(t)|AO(t+dt)>.
    #
    # Used in: overlap.py/overlap

    Norb = len(ao_i)

    S = MATRIX(Norb,Norb)

    # overlap matrix of S
    for i in range(0,Norb): # all orbitals
        for j in range(0,Norb):
            S.set(i,j,gaussian_overlap(ao_i[i],ao_j[j]))

    return S

def MO_overlap(S,ao_i, ao_j, Ci, Cj, basis_option):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] S : overlap matrix of atomic orbitals
    # \param[in] ao_i, ao_j : atomic orbital basis at different time step.
    # \param[in] Ci, Cj : molecular coefficients at different time step.
    # \param[in] basis_option : "= 1" means ab initio and "= 2" semi empirical .
    # This function returns overlap matrix of molecular orbitals with different time step
    # like <MO(t)|MO(t+dt)>.
    #
    # Used in: overlap.py/overlap

    Norb = len(ao_i)
    P = MATRIX(Norb, Norb)

    if basis_option == 1: # ab initio calculation
#        S = AO_overlap(ao_i, ao_j)
        P = Ci.T() * S * Cj
    elif basis_option == 2: # semi empirical calculation
        P = Ci.T() * Cj
    else:
        print "basis_sets has an illegal value, so stop"
        sys.exit()
    
    return P


def overlap(ao1,ao2,C1,C2,basis_sets):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] ao1, ao2 : atomic orbital basis at different time step.
    # \param[in] C1, C2 : molecular coefficients at different time step.
    # \param[in] basis_option : "= 1" means ab initio and "= 2" semi empirical .
    # This function returns overlap matrix of atomic orbitals with different time step
    # like <MO(t)|MO(t+dt)>.
    #
    # Used in: gamess_to_libra.py/gamess_to_libra
    # this is mostly a test function

    N = len(ao1)
    S11 = MATRIX(N,N); S12 = MATRIX(N,N); S21 = MATRIX(N,N); S22 = MATRIX(N,N);

    S11 = AO_overlap(ao1,ao1)
    S22 = AO_overlap(ao2,ao2)
    S12 = AO_overlap(ao1,ao2)
    S21 = AO_overlap(ao2,ao1)

    P11 = MO_overlap(S11,ao1,ao1,C1,C1,basis_sets)
    P22 = MO_overlap(S22,ao2,ao2,C2,C2,basis_sets)
    P12 = MO_overlap(S12,ao1,ao2,C1,C2,basis_sets)
    P21 = MO_overlap(S21,ao2,ao1,C2,C1,basis_sets)
    
    return P11, P22, P12, P21


def overlap_sd(MO1, MO2):
##
# \brief Compute the overlap of two determinants SD1 and SD2
# \param[in] MO1 A first set of orbitals (determinant) - the CMATRIX object of dimensions: Npw x Norb
# \param[in] MO2 A first set of orbitals (determinant) - the CMATRIX object of dimensions: Npw x Norb
# Here, Norb - the number of MOs in the set, should be the same for each object
# Npw - the number of basis functions (plane waves) used to represent the MO
#

    if MO1.num_of_rows != MO2.num_of_rows:
        print "The vertical dimensions of the two matrices do not match"
        sys.exit(0)

    if MO1.num_of_cols != MO2.num_of_cols:
        print "The horizontal dimensions of the two matrices do not match"
        sys.exit(0)

    Norb = MO1.num_of_cols  # the number of 1-el orbitals in each determinant

    ovlp = CMATRIX(Norb, Norb)
    ovlp = MO1.H()*MO2     # the overlap of the 1-el MOs

    det = []
    for i in xrange(Norb):
        row_i = []
        for j in xrange(Norb):
            row_i.append(ovlp.get(i,j))
        det.append(row_i)

    OVLP = numpy.linalg.det(det)/FACTORIAL(Norb)

    return OVLP # this is a complex number, in general


def overlap_sd_basis(sd_basis1, sd_basis2):
##
# This function constructs an overlap matrix of different Slater determinants
# from two sets: sd_basis1 and sd_basis2, which can be the same, but may be 
# different (e.g. at different time steps)
#
# \param[in] sd_basis1  A list of length n, containing different Slater determinants
#                       (each is represented by a CMATRIX object)
# \param[in] sd_basis2  -- // -- similar to the above
#

    n = len(sd_basis1)    
    m = len(sd_basis2)

    ovlp = CMATRIX(n,m)
 
    for i in xrange(n):
        for j in xrange(m):
            ovlp.set( i,j, overlap_sd(sd_basis1[i], sd_basis2[j]) )

    return ovlp  # this is a complex-valued (in general) matrix




