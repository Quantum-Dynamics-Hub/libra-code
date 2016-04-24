#*********************************************************************************                     
#* Copyright (C) 2016 Kosuke Sato, Alexey V. Akimov                                                   
#*                                                                                                     
#* This file is distributed under the terms of the GNU General Public License                          
#* as published by the Free Software Foundation, either version 2 of                                   
#* the License, or (at your option) any later version.                                                 
#* See the file LICENSE in the root directory of this distribution   
#* or <http://www.gnu.org/licenses/>.          
#***********************************************************************************
## \file tsh.py 
# This module implements the generic function for TSH calculations as well as some
# customized versions of TSH


import os
import sys
import math
import copy

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *


def surface_hopping(mol, el, ham, rnd, params):
    ## This function perform generic surface hopping
    # \param[in,out] mol a list containing Nuclear objects
    # \param[in,out] el  a list containing Electronic objects
    # \param[in,out] ham a list containing Hamiltonian objects
    # \param[in] rnd     a random number generator object. It is important that
    # we use the same (global) object every time this function is called. If we
    # create it here and use - the statistical properties will be very poor, because 
    # every new "random" number will be not far from the common seed value
    # \param[in] params  a dictionary containing control parameters, specifically:
    #
    # - params["tsh_method"] : choose hopping probability calculation scheme
    #    1 - FSSH, 2 - GFSH, 3 - MSSH
    # - params["rep"] : choose the representation for velocity rescaling
    #    0 - diabatic, 1 - adiabatic
    # - params["do_rescaling"] : choose how to account for detailed balance
    #    0 - don't do explicit rescaling, so it is used with Boltzmann factor hopping probability scaling
    #    1 - do excplicit velocity rescaling
    # - params["do_reverse"] : how to handle the nuclear velocities when hop is frustrated
    #    0 - keep as they are, don't change
    #    1 - reverse the velocities
    # - params["dt_nucl"] : nuclear time step in fs
    # - params["Temperature"] : Temperature of the environment, K
    # - params["print_tsh_probabilities"] : Print hopping probabilities matrix
    #    0 - don't print,  1 - print
    # - params["check_tsh_probabilities"] : Run-time sanity test of the dt_nucl, 
    #    0 - don't check,  1 - check
    # - params["use_boltz_factor"] : whether to scale the hopping probabilities by the Boltzmann factor
    #    0 - don't scale,  1 - scale


    # Parameters to internal variables - for convenience
    tsh_method = params["tsh_method"]
    rep = params["rep"]
    do_rescaling = params["do_rescaling"]
    do_reverse = params["do_reverse"]
    dt_nucl = params["dt_nucl"]
    Temperature = params["Temperature"]
    print_prob = params["print_tsh_probabilities"]
    check_prob = params["check_tsh_probabilities"]
    use_boltz_factor = params["use_boltz_factor"]
    

    # Parameters characterizing the system and the ensemble
    if len(el)!=len(mol):
        print "Error in surface_hopping: The size of the ensemble of Electronic objects (",len(el),") should be the same\
              as the length of the ensemble of Nuclear objects (", len(mol),")\n"             
    if len(el)!=len(ham):
        print "Error in surface_hopping: The size of the ensemble of Electronic objects (",len(el),") should be the same\
              as the length of the ensemble of Hamiltonian objects (", len(mol),")\n"             

    ntraj = len(el)          # this is the number of trajectories in an ensemble
    nstates = el[0].nstates  # how many electronic DOF


    #Compute hopping probabilities
    g = MATRIX(nstates,nstates) # initialize a matrix of hopping probability

    if tsh_method == 1: # FSSH
        compute_hopping_probabilities_fssh(mol, el, ham, g, dt_nucl, use_boltz_factor, Temperature)
    elif tsh_method == 2: # GFSH
        compute_hopping_probabilities_gfsh(mol, el, ham, g, dt_nucl, use_boltz_factor, Temperature)
    elif tsh_method == 3: # MSSH
        compute_hopping_probabilities_mssh(mol, el, ham, g, dt_nucl, use_boltz_factor, Temperature)
    else:
        print "Warning in surface_hopping: tsh_method can be 1, 2, or 3. Other values are not defined"


    # output hopping probability
    if print_prob == 1:
        print "hopping probability matrix is:"
        print g.show_matrix()

    # check elements of g matrix are less than 1 or not.
    if check_prob == 1:
        for st in xrange(nstates):
            for st1 in xrange(nstates):
                if g.get(st,st1) > 1:
                    print "g(%d,%d) is %f, larger than 1; better to decrease dt_nucl" %(st,st1,g.get(st,st1))

    # Attempt to hop
    ksi = rnd.uniform(0.0,1.0) # generate random number for every trajectory   

    # Everything else - change of electronic state and velocities rescaling/reversal happens here     
    el.istate = hop(el.istate, mol, ham, ksi, g, do_rescaling, rep, do_reverse)


    # Nothing to return - mol, ham, and el objects are modified accordingly



def surface_hopping_cpa(mol, el, ham, rnd, params):
    ## This function performs surface hopping with Boltzmann factor used to 
    # rescale hopping probabilities in lieu of explicit velocity rescaling

    # Update parameters
    params["do_rescaling"] = 0      # No explicit velocity rescaling

    params["use_boltz_factor"] = 1  # we don't need to use Boltzmann factor, since  
                                    # we are using velocity rescaling in the hopping procedure.
                                    # Although the rescaling doesn't account for the direction, but it
                                    # still accounts for energy partitioning between electronic and
                                    # nuclear DOFs

    # rep and do_reverse are irrelevant


    # Call actual calculations 
    surface_hopping(mol, el, ham, rnd, params)



def surface_hopping_cpa2(mol, el, ham, rnd, params):
    ## This function performs surface hopping with velocity rescaling according to 
    # total energy conservation (but not along the derivative coupling vectors)

    # Update parameters
    params["do_rescaling"] = 1      # Explicit velocity rescaling

    params["use_boltz_factor"] = 0  # we don't need to use Boltzmann factor, since  
                                    # we are using velocity rescaling in the hopping procedure.
                                    # Although the rescaling doesn't account for the direction, but it
                                    # still accounts for energy partitioning between electronic and
                                    # nuclear DOFs

    params["rep"] = 0               # velocity rescaling will be done based on the total energy conservation,
                                    # no derivative couplings will be needed - we don't have them
                                    # !!! This option makes do_reverse not relevant - so
                                    # we can set it to any value
    # do_reverse becomes irrelevant


    # Call actual calculations 
    surface_hopping(mol, el, ham, rnd, params)





