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

## \file states.py
# This module defines the function which creates a list of ground and excited states.
# It outputs the key parameter named "excitation"

# ************************ CAUTION ***************************************
# This function creates only singlet type excitation states.
# ************************************************************************

import os
import sys
import math

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *


def pyxaid_states(states, min_shift, max_shift):
    ##
    # This function converts our Libra type states into Pyxaid convention
    # E.g.
    # HOMO = 3, min_shift = -1, max_shift = 1
    # [0,1,2,3,4,5,6,7]  <- original set
    #     [0,1,2]        <- reduced set
    # nocc = 2
    # active_space_sz = 3
    # so: gs = [1,-1,2,-2]
    #
    # Used in: hamiltonian_vib.py/update_vibronic_hamiltonian
    
    active_space_sz = max_shift - min_shift + 1
    nstates = len(states)

    nocc = -min_shift + 1 # the number of double occupied orbitals

    gs = [] # ground states
    for i in xrange(nocc):
        indx = i + 1  # +1 so we can distinguish 0 and -0
        gs.append(indx)
        gs.append(-indx)

    # Now lets create all excitations from the reference state
    res = []
    print "All generated configurations (pyxaid indexing)\n"
    for i in xrange(nstates):
        es = list(gs)
        print "Configuration ", i, es

        # Here +1 is needed to be consistent with Pyxiad indexing convention
        a = states[i].from_orbit[0] - min_shift + 1
        a_s = states[i].from_spin[0]  # +1 for alp, -1 for bet
        b = states[i].to_orbit[0] - min_shift + 1
        b_s = states[i].to_spin[0]  # +1 for alp, -1 for bet

        # orbital indices with sign = spin-orbit indices
        A = a * a_s
        B = b * b_s

        # Replace A with B
        indx = es.index(A)
        es[indx] = B

        #print "Replace orbital ",A, " at position ",indx," with orbital ",B

        res.append(es)

    return res


def create_states(Nmin,HOMO,Nmax,spin,flip):
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in]  Nmin  lowest molecular orbital taken for TD-SE calculation
    # \param[in]  HOMO  Highest Occupied Molecular Orbital
    # \param[in]  Nmin  Highest molecular orbital taken for TD-SE calculation
    # \param[in]  spin  spin is considered : option 0 -> no, 1 -> yes
    # \param[in]  flip  spin flip is considered if spin = 1: option 0 -> no, 1 -> yes
    # excitations - returned list of "excitation" objects in active space.
    #
    # Used in:  run.py

    LUMO = HOMO + 1
    
    excitations = []

    if spin == 0:
        sp_st = [1]
    elif spin == 1:
        sp_st = [1,-1]
    else:
        print "Error in create_states: spin must be 0 or 1"
        print "Value given = ", spin
        print "Exiting..."
        sys.exit(0)

    # note "excitation" object.
    # "excitation" has 4 index : (from orbital, from spin, to orbital, to spin)
    # \orbital index  0 -> HOMO, 1 -> LUMO, 2 -> LUMO+1, -1 -> HOMO-1, etc....
    # \spin index 1 -> alpha -1 -> beta
    # e.g. (0,1,0,1) means ground state (no excitation)
    #      (0,1,1,1)       alpha electron in HOMO is excited to LUMO without spin flip
    #      (0,1,1,-1)      alpha electron in HOMO is excited to LUMO with spin flip (in this case, spin-orbital coupling should be included)
    
    print "GS: HOMO-0,alpha -> HOMO-0,alpha"
    excitations.append(excitation(0,1,0,1)) # Add a Ground State

    if flip == 0: # Excited States without spin flip
        for sp in sp_st:
            if sp == 1:
                e_spin = "alpha"
            else:
                e_spin = "beta"
            icount = 0
            for om in range(HOMO,Nmin-1,-1):
                
                for uom in range(LUMO,Nmax+1,1):

                    print "SE%i: HOMO-%i,%s -> LUMO+%i,%s" %(icount,abs(om-HOMO),e_spin,uom-HOMO-1,e_spin)
                    excitations.append(excitation(om-HOMO,sp,uom-HOMO,sp))
                    icount += 1

    elif flip == 1 and spin == 1: # Excited states with spin flip
        for sp1 in sp_st:
            for sp2 in sp_st:
                for om in range(Nmin,HOMO+1):
                    for uom in range(LUMO,Nmax+1):

                        excitations.append(excitation(om-HOMO,sp1,uom-HOMO,sp2))

    else:
        print "Error in create_states: flip must be 0 or 1 and spin 0 or 1"
        print "And,when spin = 0, flip should not be 1"
        print "Value flip and spin given = ", flip, spin
        print "Exiting..."
        sys.exit(0)

    print "excitations length is ",len(excitations)
    #print excitations

    return excitations
