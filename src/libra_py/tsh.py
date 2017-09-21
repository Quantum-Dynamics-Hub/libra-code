#*********************************************************************************                     
#* Copyright (C) 2016-2017 Kosuke Sato, Alexey V. Akimov                                                   
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
#
# The module contain the following functions:
#
#   hop_py(initstate, g, ksi)
#   set_random_state(prob, ksi)
#   compute_sh_statistics(nstates, istate)
#   avarage_populations(el)
#   surface_hopping(mol, el, ham, rnd, params)
#   surface_hopping_cpa(mol, el, ham, rnd, params)
#   surface_hopping_cpa2(mol, el, ham, rnd, params)
#   ida_py(Coeff, old_st, new_st, E_old, E_new, T, ksi, do_collapse)


import os
import sys
import math
import copy

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *


def hop_py(initstate, g, ksi):
    ##
    # This function implements a simple surface hopping procedure
    # \param[in]   initstate [ integer ] The state index before hop  
    # \param[in]           g [ MATRIX ] The surface hopping matrix, the element g(i,j) contains the probability for a i->j transition
    # \param[in]         ksi [ float ] A random number uniformly distributed in the range of (0.0, 1.0) 

    # The function returns:
    # finstate  [ integer ] The index of the final state after hop


    nstates = g.num_of_cols
    finstate = initstate;

    left, right = 0.0, 0.0

    for i in xrange(nstates):
      if i==0:
          left = 0.0
          right = g.get(initstate,i)
      else:
          left = right
          right = right + g.get(initstate,i)
 
      if((left<ksi) and (ksi<=right)):
          finstate = i

    return finstate


def set_random_state(prob, ksi):
    ##
    # This function implements a simple random state selection procedure. Each state is selected with a given probability
    # \param[in]   prob [ list of floats ] The probabilities of all states 
    # \param[in]   ksi [ float ] A random number uniformly distributed in the range of (0.0, 1.0) 

    # The function returns:
    # finstate  [ integer ] The index of the selected state


    nstates = len(prob)
    finstate = 0;

    left, right = 0.0, 0.0

    for i in xrange(nstates):
      if i==0:
          left = 0.0
          right = prob[i]
      else:
          left = right
          right = right + prob[i]
 
      if((left<ksi) and (ksi<=right)):
          finstate = i

    return finstate


def compute_sh_statistics(nstates, istate):
    ##
    # This function computes the SH statistics for an ensemble of trajectories
    # \param[in]   nstates [ integer ] The number of allowed quantum state
    # \param[in]   istate [ list of integers ] The list containing the info about the index of quantum state in which each trajectory is found
    # The length of the list is equal to the number of trajectories. Each element of the list is the state index (integer)

    # The function returns:
    # coeff_sh [ list of integers ] The list containing the average population of each quantum state. The length of the list is equal to the 
    # total number of quantum states considered


    num_sh_traj = len(istate)
    f = 1.0/float(num_sh_traj)

    coeff_sh = MATRIX(nstates, 1)

    for i in xrange(num_sh_traj):
        st = istate[i]
        coeff_sh.set(st, coeff_sh.get(st) + f)
 
    return coeff_sh


def update_sh_pop( states , nstates):
    ##
    # states - is a vector of state index of each trajectory
    #  so len(states) - the number of trajectories 
    #  and states[j] - the state of the trajectory j
    # nstates - the number of the states possible
    #
    # Returns the SH-based population of all states

    pops = [0.0] * nstates
    ntraj = len(states)

    incr = 1.0/float(ntraj)

    for j in xrange(ntraj): # for all trajectories
        pops[ states[j] ] += incr

    return pops



    

def avarage_populations(el):
    ##
    # This function computes the SH statistics for an ensemble of trajectories
    # \param[in]   el [ list of Electronic ] The list containing electronic DOF variables for all trajectories in ensemble.
    # The length of the list determines the number of trajectories in ensemble

    # The function returns:
    # sh_pops [ list of float ] The list containing the average population of each quantum state based on the statistics of the
    # discrete states in which each trajectory resides. 
    # se_pops [ list of float ] The list containing the average population of each quantum state based on the amplitudes of all quantum states
    # as obtained from the TD-SE solution. 
    # rho [CMATRIX] The matrix containing trajectory-averaged SE populations and coherences.

    #The length of the list is equal to the total number of quantum states considered

    ntraj = len(el)        # the total number of trajectories
    nstat = el[0].nstates  # the number of quantum states, assume that all objects in the "el" list 
                           # are similar

    sh_pops = [0.0] * nstat   # average SH populations of all states
    se_pops = [0.0] * nstat   # average SE populations of all states
    rho = CMATRIX(nstat, nstat) # trajectory-averaged density matrix

    f = 1.0/float(ntraj)
    for traj in xrange(ntraj): # for all trajectories
        sh_pops[ el[traj].istate ] += f

        for st1 in xrange(nstat):
            se_pops[ st1 ] += f * el[traj].rho(st1,st1).real

            for st2 in xrange(nstat):
                rho.set(st1, st2, rho.get(st1,st2) + f * el[traj].rho(st1,st2) )

    return sh_pops, se_pops, rho


def ave_pop(denmat_sh, denmat_se):
# 
# \param[in] denmat_sh (list of CMATRIX(nst_in, nst_in)) Vector with the density matrix (diagonal in SH) for each trajectory
# \param[in] denmat_se (list of CMATRIX(nst_in,nst_in)) Vector with the SE density matrix for each trajectory 
#
#  Returns: Ensemble averaged SH and SE density matrices (CMATRIX(nst_out, nst_out) each) 
#  

    ntraj = len(denmat_sh)
    nst_out = denmat_sh[0].num_of_cols

    ave_pop_sh = CMATRIX(nst_out, nst_out)
    ave_pop_se = CMATRIX(nst_out, nst_out)
    den = 1.0/float(ntraj)

    for i in xrange(ntraj):
        ave_pop_se = ave_pop_se + den * denmat_se[i]   # SE
        ave_pop_sh = ave_pop_sh + den * denmat_sh[i]   # SH


    return ave_pop_sh, ave_pop_se




def amplitudes2denmat(coeffs):
# \param[in] coeffs (list of CMATRIX(nstates, 1)) wavefunction amplitudes for all trajectories

    ntraj = len(coeffs)
    denmat = []

    for tr in xrange(ntraj):
        denmat.append( coeffs[tr] * coeffs[tr].H() )

    return denmat


def denmat2prob(P):
# \param[in] P (CMATRIX) Density matrix
#
    nst = P.num_of_cols
    prob = [0.0] * nst

    for i in xrange(nst):
        prob[i] = P.get(i,i).real

    return prob





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




    g = MATRIX(nstates,nstates) # initialize a matrix of hopping probability

    for i in xrange(ntraj):

        #Compute hopping probabilities
        if tsh_method == 1: # FSSH
            compute_hopping_probabilities_fssh(mol[i], el[i], ham[i], g, dt_nucl, use_boltz_factor, Temperature)
        elif tsh_method == 2: # GFSH
            compute_hopping_probabilities_gfsh(mol[i], el[i], ham[i], g, dt_nucl, use_boltz_factor, Temperature)
        elif tsh_method == 3: # MSSH
            compute_hopping_probabilities_mssh(mol[i], el[i], ham[i], g, dt_nucl, use_boltz_factor, Temperature)
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
        el[i].istate = hop(el[i].istate, mol[i], ham[i], ksi, g, do_rescaling, rep, do_reverse)


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



def hopping(Coeff, Hvib, istate, sh_method, do_collapse, ksi, ksi2, dt, T):
    """
    A simplified version for the CPA-like hopping

    Coeff (CMATRIX(nstates, 1) ) object with the amplitudes of all states
    Hvib (CMATRIX(nstates, nstates) )  object containing the vibronic Hamiltonian 
    istate (int) the index of the initial state
    sh_method (int) selector of the TSH method: 0 - MSSH, 1 - FSSH
    do_collapse (int) flag to turn on the decoherence via ID-A: 0 - no decoherence, 1 - decoherence via ID-A
    ksi, ksi2 (float in [0, 1]) random numbers cotrolling the execution of SH
    dt (float) time interval for the surface hopping (in a.u.)
    T (float) temperature in K

    Returns: the index (int) of a new state 

    """
    g = 0.0
    if sh_method==0:
        g = compute_hopping_probabilities_mssh(Coeff)
    elif sh_method==1:
        g = compute_hopping_probabilities_fssh(Coeff, Hvib, dt)

    old_st = istate
    new_st = hop(istate, g, ksi)

    if new_st != old_st:
        E_old = Hvib.get(old_st,old_st).real
        E_new = Hvib.get(new_st,new_st).real

        # ID-A decoherence                
        istate, Coeff1 = tsh.ida_py(Coeff, old_st, new_st, E_old, E_new, T, ksi2, do_collapse)


    return istate #, Coeff1
    






def ida_py(Coeff, old_st, new_st, E_old, E_new, T, ksi, do_collapse):

    ##
    # This function implements the decoherence correction 
    # \param[in]       Coeff [ CMATRIX or Electronic ] An object containig electronic DOFs. 
    # \param[in]      old_st [ integer ] The state index before hop  
    # \param[in]      new_st [ integer ] The state index after hop 
    # \param[in]       E_old [ float ] The energy of the initial state, before hop 
    # \param[in]       E_new [ float ] The energy of the final state, after hop 
    # \param[in]           T [ float ] The Temperature of nuclear DOF  
    # \param[in]         ksi [ float ] A random number uniformly distributed in the range of (0.0, 1.0) 
    # \param[in] do_collapse [ 0 or 1 ] The flag turning the decoherence (at the IDA level on/off). 1 - include decoherence, 0 - do not include decoherence 

    # The function returns:
    # res [ integer ] - index of the final state, after IDA is applied (or not)
    # C [CMATRIX or Electronic] - the updated state of the electronic DOF, in the same data type as the input

    kb = 3.166811429e-6  # Hartree/K
    res = old_st
    dE = (E_new - E_old)

    # Compute the Boltzmann scaling factor, but only if we consider a hop up in energy
    boltz_f = 1.0   
    if dE > 0.0:
        argg = dE/(kb*T)
        if argg > 50.0:
            boltz_f = 0.0
        else:
            boltz_f = math.exp(-argg)


    # In case the electronic DOF are given in the form of CMATRIX
    if type(Coeff).__name__ == "CMATRIX":

        C = CMATRIX(Coeff)
        
        if dE>0.0:
        
            if ksi<boltz_f:
                res = new_st  # accepted hop
                
                # Collapse the wavefunction to the new state 
                if do_collapse:
                    C *= 0.0; C.set(new_st, 1.0+0.0j)
            else:
                # Unsuccessful hop - collapse wfc back to the original state^M
                if do_collapse:
                    C *= 0.0; C.set(old_st, 1.0+0.0j)
        else:
            res = new_st
        
        return res, C


    # In case the electronic DOF are given in the form of Electronic object
    elif type(Coeff).__name__ == "Electronic":
        
        C = Electronic(Coeff)
                
        if dE>0.0:
        
            if ksi<boltz_f:
                res = new_st  # accepted hop
                # Collapse the wavefunction to the new state                                                                                                           
                for st in xrange(C.nstates):
                    C.q[st], C.p[st] = 0.0, 0.0
                C.q[new_st], C.p[new_st] = 1.0, 0.0
            else:
                # Unsuccessful hop - collapse wfc back to the original state^M                                                                                         
                for st in xrange(C.nstates):
                    C.q[st], C.p[st] = 0.0, 0.0
                C.q[old_st], C.p[old_st] = 1.0, 0.0
        else:
            res = new_st
        
        C.istate = res # return res
        
        return res, C


def sdm_py(Coeff, dt, act_st, En, Ekin, C_param = 1.0, eps_param = 0.1):

    ##
    # This function implements the simplified decay of mixing algorithm for decoherence correction
    # Reference: Granucci, G.; Persico, M. J. Chem. Phys. 2007, 126, 134114
    #
    # \param[in]       Coeff [ CMATRIX or Electronic ] An object containig electronic DOFs. 
    # \param[in]          dt [ float ] The integration timestep. Units = a.u. of time
    # \param[in]      act_st [ integer ] The active state index
    # \param[in]       En    [ list of floats ] Energies of the states. Units = Ha
    # \param[in]        Ekin [ float ] The classical kinetic energy of nuclei. Units = Ha
    # \param[in]     C_param [ float ] The method parameter, typically set to 1.0 Ha
    # \param[in]   eps_param [ float ] The method parameter, typically set to 0.1 Ha

    # The function returns:
    # C [CMATRIX or Electronic] - the updated state of the electronic DOF, in the same data type as the input

    kb = 3.166811429e-6  # Hartree/K
   

    # In case the electronic DOF are given in the form of CMATRIX
    if type(Coeff).__name__ == "CMATRIX":

        # The results will be stored here
        C = CMATRIX(Coeff)


        # First - update all the coefficients for the non-active states        
        N = Coeff.num_of_elts 
        for i in xrange(N):
            if i != act_st:    
                itau = ( En[i] - En[act_st] ) / ( C_param + (eps_param/Ekin) )
                sclf = math.exp(-dt*itau)
                C.scale(i, 0, sclf)

        # Population of the active state
        p_aa_old = (C.get(act_st,act_st).conjugate * C.get(act_st,act_st)).real 

        new_norm = (C.H() * C).get(0,0).real - p_aa_old  # total population of all inactive states
                                                         # after rescaling
        p_aa_new = 1.0 - new_norm

        sclf = 1.0
        if p_aa_old > 0.0:
            sclf = math.sqrt( p_aa_new / p_aa_old )  # scaling factor for the active state
        

        # Rescale the active state
        C.scale(act_st, 0, sclf)
        
        return C

    # In case the electronic DOF are given in the form of Electronic object
    elif type(Coeff).__name__ == "Electronic":
        
        C = Electronic(Coeff)

        # First - update all the coefficients for the non-active states        
        N = C.nstates 
        new_norm = 0.0
        for i in xrange(N):
            if i != act_st:    
                itau = ( En[i] - En[act_st] ) / ( C_param + (eps_param/Ekin) )
                sclf = math.exp(-dt*itau)
                C.q[i] = C.q[i] * sclf
                C.p[i] = C.p[i] * sclf

                new_norm += C.rho(i, i).real 

        # new_norm now contains the total population of all inactive states after rescaling
        # How much of population is left for the new active state
        p_aa_new = 1.0 - new_norm

        sclf = 1.0
        if p_aa_old > 0.0:
            sclf = math.sqrt( p_aa_new / p_aa_old )  # scaling factor for the active state

        # Rescale the active state
        C.q[act_st] = C.q[act_st] * sclf
        C.p[act_st] = C.p[act_st] * sclf


        return C

        

