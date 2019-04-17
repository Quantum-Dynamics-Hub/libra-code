#*********************************************************************************
#* Copyright (C) 2016 Ekadashi Pradhan, Alexey V. Akimov
#* 
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/

## \file spin_index.py 
# This module implements the functions that extract spin index of excitation object

import os
import sys
import math

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *


def index_spin(ex,active_space,homo):
##
# This function creates alpha (up-spin) and beta (down spin) index list
# which are later used to create SD (Slater Determinant) data type
# \param[in] params Parameter dictionary, particularly, params["excitations"] and params["nel"]
#                   is used in this function
# \param[in] active_space List of active space orbital number
# \param[out] alp list of alpha spin index
# \param[out] bet list of beta spin index 
#
    #print "in index_spin..."

    # In case of non-spin-polarized calculation
    # Index of HOMO orbital in the active space, e.g., if active space is [4,5,6,7,8] and HOMO
    # orbital is 6, then the h_idx of HOMO in the active space is 2, new index in the active space is
    # [0,1,2,3,4]

    #print "homo = ", homo
    #print "active_space = ", active_space

    h_idx = active_space.index(homo)
    alp,bet = [],[]

        # use excitation object to create proper SD object for different excited state
        # Currently, excitation only from HOMO is included
        # Consider 2 electrons in the active space:
        # Ground state: 2 electrons are sitting on HOMO (HOMO-1 is not included into our active space)
        # Schematics: U - spin-up, D - spin-down, I - inactive electrons
        #
        #       SD0            SD1           SD2
        #
        #  3  ---------     ---------      ---------
        #  2  ---------     --   D --      -- U   --
        #  1  -- U D --     -- U   --      -- U   --
        #  0  -- I I --     -- I I --      -- I I --

    if ex.size>1:
        print "Only single excitations are currently allowed. Exiting...\n"
        sys.exit(0)

    # Ground state "excitation"
    if ex.from_orbit[0] == ex.to_orbit[0]:
        alp.append(ex.to_orbit[0] + h_idx)
        bet = alp

    elif ex.from_orbit[0] != ex.to_orbit[0]:

        # spin flip excitations
        if ex.from_spin[0] != ex.to_spin[0]:
            if ex.to_spin[0] == -1:
                bet.append(ex.from_orbit[0] + h_idx)
                bet.append(ex.to_orbit[0] + h_idx)
            elif ex.to_spin[0] == 1:
                alp.append(ex.from_orbit[0] + h_idx)
                alp.append(ex.to_orbit[0] + h_idx)

        # same spin excitations
        else:
            if ex.to_spin[0] == 1:
                alp.append(ex.to_orbit[0] + h_idx)
                bet.append(ex.from_orbit[0] + h_idx)
            elif ex.to_spin[0] == -1:
                bet.append(ex.to_orbit[0] + h_idx)
                alp.append(ex.from_orbit[0] + h_idx)


    #print "alp = ", alp
    #print "bet = ", bet
    return alp, bet


def truncate_spin(ex,alp,bet,nel_SD,homo_idx):
##
# This function further reduce the mo_pool and selects the number of
# electrons in the SD
#
    #print "begining spin index=",alp,bet,"homo idx=",homo_idx
    if nel_SD % 2 == 0: #Closed Shell System
        #print "closed shell system"
        if nel_SD == 2:
            #print "Number of Electron in SD =",nel_SD
            #print "from orb=",ex.from_orbit[0],"from spin=",ex.from_spin[0],"to orbit=",ex.to_orbit[0],"to spin=",ex.to_spin[0]
            alp,bet = [],[]
            if ex.from_spin[0] != ex.to_spin[0]:   #Triplet state
                if ex.to_spin[0] == 1:
                    #print "alp=",alp,"homo_idx=",homo_idx,ex.from_orbit[0]
                    alp.append(homo_idx + ex.from_orbit[0])
                    alp.append(homo_idx + ex.to_orbit[0])
                    #print "alp loop 1=",alp,bet
                if ex.to_spin[0] == -1:
                    bet.append(homo_idx + ex.from_orbit[0])
                    bet.append(homo_idx + ex.to_orbit[0])
                #print "working"
            else:    # Singlet state
                if ex.to_spin[0] == 1:
                    alp.append(homo_idx + ex.to_orbit[0])
                    bet.append(homo_idx + ex.from_orbit[0])
                if ex.to_spin[0] == -1:
                    alp.append(homo_idx + ex.from_orbit[0])
                    bet.append(homo_idx + ex.to_orbit[0])
    #print alp,bet
    return alp,bet



        


def index_spin_gen(ex,active_space,nel,nel_SD):
##
# This function creates alpha (up-spin) and beta (down spin) index list
# which are later used to create SD (Slater Determinant) data type
# \param[in] ex is params["excitations"] 
# \param[in] is params["nel"], used to find out if its and open shell or closed shell system
#            and also gives HOMO orbital number
# \param[in] active_space List of active space orbital number
# \param[out] alp list of alpha spin index
# \param[out] bet list of beta spin index 
#
    #print "in index_spin..."

    # In case of non-spin-polarized calculation
    # Index of HOMO orbital in the active space, e.g., if active space is [4,5,6,7,8] and HOMO
    # orbital is 6, then the h_idx of HOMO in the active space is 2, new index in the active space is
    # [0,1,2,3,4]
    homo = nel/2 +  nel % 2
    open_shell = nel % 2   # if 0 then closed shell, if 1 then open shell
    #print "homo = ", homo
    #origin_ex = homo + ex.from_orbit[0]
    #print "active_space = ", active_space

    homo_idx = active_space.index(homo)
    h_idx = homo_idx + 1
    #h_idx = active_space.index(origin_ex) + 1 # For python indexing
    alp,bet = [],[]

        # use excitation object to create proper SD object for different excited state
        # Currently, single excitation only from all possible active space orbitals are included
        # Consider 2 electrons in the active space:
        # Ground state: 2 electrons are sitting on HOMO (HOMO-1 is not included into our active space)
        # Schematics: U - spin-up, D - spin-down, I - inactive electrons
        #
        #       SD0            SD1           SD2
        #
        #  3  ---------     ---------      ---------
        #  2  ---------     --   D --      -- U   --
        #  1  -- U D --     -- U   --      -- U   --
        #  0  -- I I --     -- I I --      -- I I --

    if ex.size>1:
        print "Only single excitations are currently allowed. Exiting...\n"
        sys.exit(0)

    if open_shell == 0:
        # Closed shell system
        # Ground state "excitation"
        #___________________________________
        for i in xrange(h_idx):
            alp.append(i)
            bet.append(i)
            #-----------------------------------
        if ex.from_orbit[0] != ex.to_orbit[0]:
            # spin flip excitations
            if ex.from_spin[0] != ex.to_spin[0]:
                if ex.to_spin[0] == -1:
                    alp.remove(homo_idx + ex.from_orbit[0])
                    bet.append(homo_idx + ex.to_orbit[0])
                elif ex.to_spin[0] == 1:
                    bet.remove(homo_idx + ex.from_orbit[0])
                    alp.append(homo_idx + ex.to_orbit[0])
            # same spin excitations
            else:
                if ex.to_spin[0] == 1:
                    alp.remove(homo_idx + ex.from_orbit[0])
                    alp.append(homo_idx + ex.to_orbit[0])
                elif ex.to_spin[0] == -1:
                    bet.remove(homo_idx + ex.from_orbit[0])
                    bet.append(homo_idx + ex.to_orbit[0])
        #print alp,bet
        alp,bet = truncate_spin(ex,alp,bet,nel_SD,homo_idx)

    elif open_shell == 1:
        # Closed shell system
        # Ground state "excitation"
        #___________________________________
        for i in xrange(h_idx-1):
            alp.append(i)
            bet.append(i)
        alp.append(h_idx-1)  # Right now, the unpaired electron is in alpha spin
            #-----------------------------------

        if ex.from_orbit[0] != ex.to_orbit[0]:

            if ex.to_orbit[0] == 0 and ex.to_spin[0] == 1:
                print "Currently, the GS is snglet (spin up), no two electron in same spin in an orbital ... exiting"
                sys.exit(0)

            # spin flip excitations
            elif ex.from_spin[0] != ex.to_spin[0]:
                if ex.to_spin[0] == -1:
                    alp.remove(homo_idx + ex.from_orbit[0])
                    bet.append(homo_idx + ex.to_orbit[0])
                elif ex.to_spin[0] == 1:
                    bet.remove(homo_idx + ex.from_orbit[0])
                    alp.append(homo_idx + ex.to_orbit[0])

            # same spin excitations
            else:
                if ex.to_spin[0] == 1:
                    alp.remove(homo_idx + ex.from_orbit[0])
                    alp.append(homo_idx + ex.to_orbit[0])
                elif ex.to_spin[0] == -1:
                    bet.remove(homo_idx + ex.from_orbit[0])
                    bet.append(homo_idx + ex.to_orbit[0])

    #print "alp = ", alp
    #print "bet = ", bet
    return alp, bet


def index_spin_test1():
    params = {}
    params["excitations"] = [ excitation(0,1,0,1), 
                              excitation(0,1,1,1), 
                              excitation(0,-1,1,-1),
                              excitation(0,1,1,-1),
                              excitation(0,-1,1,1)
                            ]
    params["nel"] = 12
    homo = params["nel"]/2 +  params["nel"] % 2

    nel = 12 # for using index_spin_gen
    active_space = [6,7]  # HOMO, LUMO
    for ex in params["excitations"]:
        index_spin(ex, active_space, homo)


def index_spin_test2():
    params = {}
    params["excitations"] = [ excitation(0,1,0,1),
                              excitation(0,1,1,1),
                              excitation(-1,1,1,1),
                              excitation(-2,1,1,1),
                              excitation(-2,-1,1,1)
                            ]
    nel = 12 # for using index_spin_gen
    nel_SD = 2
    active_space = [4,5,6,7]  # HOMO, LUMO
    for ex in params["excitations"]:
        index_spin_gen(ex, active_space, nel, nel_SD)


# For testing and debug purposes
#index_spin_test1()
#index_spin_test2()


