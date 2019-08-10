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
## \file hamiltonian_el.py
# This module defines the functions that return time-averaged energy and
# the Non-Adiabatic couplings (NACs).

import os
import sys
import math

from overlap import *

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *


def compute_nac_sd(MO_old, MO_cur, dt):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] MO_old A list of MO sets, each defining SD for a given electronic state: at time t-dt
    # \param[in] MO_cur A list of MO sets, each defining SD for a given electronic state: at time t
    # \param[in] dt nuclear time step
    # NAC - non-adiabatic coupling matrix
    #
    # Used in: x_to_libra.py/gamess(g09)_to_libra

    # Although the sets MO1 and MO2 are not mutually-orthogonal, so there would be a dS/dt contribution,
    # we compute only the Hermitian part, since the non-Hermitian will cancel out in the solving TD-SE


    nstates = len(MO_cur)  # the number of electronic states
    NAC = CMATRIX(nstates,nstates)

    for i in xrange(nstates):
        for j in xrange(nstates):
            s_01 = overlap_sd(MO_old[i],MO_cur[j])   # <SD_i(t-dt)|SD_j(t)>
            s_10 = overlap_sd(MO_cur[i],MO_old[j])   # <SD_i(t)|SD_j(t-dt)>
            NAC.set(i,j,(s_01 - s_10)/(2.0*dt))

    return NAC


def NAC(P12,P21,dt_nucl):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] P12, P21 : overlap matrix of molecular orbitals at different time step.
    # \param[in] dt_nucl  : time step width of nuclear motion
    # This function returns Non-Adiabatic Couplings(NACs)
    #
    ##### Used in: x_to_libra_gms.py/gamess_to_libra

    Norb = P12.num_of_rows
    D = MATRIX(Norb,Norb)

    D = 0.50/dt_nucl * ( P12 - P21 )

    return D

def average_E(E1,E2):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] E1, E2 : molecular energies at different time step.
    # This function returns the time-averaged molecular energies.
    #
    # Used in: x_to_libra.py/gamess_to_libra

    Norb = E1.num_of_rows
    E = MATRIX(Norb,Norb)

    E = 0.50 * (E1 + E2)

    return E

def average_S(S1,S2):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] S1, S2 : Overlap at different time step.
    # This function returns the time-averaged overlap matrix S(t+dt/2).
    #
    # Used in: x_to_libra.py/gamess_to_libra

    Norb = S1.num_of_rows
    S = MATRIX(Norb,Norb)

    S = 0.50 * (S1 + S2)

    return S

