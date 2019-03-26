#*********************************************************************************                     
#* Copyright (C) 2016-2019 Kosuke Sato, Alexey V. Akimov                                                   
#*                                                                                                     
#* This file is distributed under the terms of the GNU General Public License                          
#* as published by the Free Software Foundation, either version 2 of                                   
#* the License, or (at your option) any later version.                                                 
#* See the file LICENSE in the root directory of this distribution   
#* or <http://www.gnu.org/licenses/>.          
#***********************************************************************************
"""
.. module:: tsh
   :platform: Unix, Windows
   :synopsis: This module implements the generic function for TSH calculations as well as some
       customized versions of TSH recipes
.. moduleauthor:: Kosuke Sato, Alexey V. Akimov

"""

__author__ = "Alexey V. Akimov, Kosuke Sato"
__copyright__ = "Copyright 2016-2019 Kosuke Sato, Alexey V. Akimov"
__credits__ = ["Alexey V. Akimov", "Kosuke Sato"]
__license__ = "GNU-3"
__version__ = "1.0"
__maintainer__ = "Alexey V. Akimov"
__email__ = "alexvakimov@gmail.com"
__url__ = "https://quantum-dynamics-hub.github.io/libra/index.html"



import os
import sys
import math
import copy

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

import units
import probabilities


def sample(x, mean_x, sigma_x, rnd):  
    """
    This function generates ntraj ndof-dimensional vectors sampled from a 
    normal distribution with a given mean and variance

    Args: 
        x ( MATRIX(ndof, ntraj) ): Each column of the matrix corresponds to 
            a vector of certain properties (e.g. coordinates, momenta, of all DOFs) for 
            a given trajectory (element of ensemble)
        mean_x ( MATRIX(ndof, 1) ):  The mean of the ndof-dimensional vector (component-wise)
        sigma_x ( MATRIX(ndof, 1) ): The variance width for each component
        rnd ( Random ): The random number generator object

    Returns:
        None: but changes the matrix ```x```

    """
    nr, nc = x.num_of_rows, x.num_of_cols
    for i in range(nr):
        for j in range(nc):    
            x.set(i,j, mean_x.get(i,0) + sigma_x.get(i,0) * rnd.normal() )



def hop_py(initstate, g, ksi):
    """
    The Python implementation of the stochastic hopping procedure

    Args:
        initstate ( int ): The index of the initial state, before hop  
        g ( MATRIX(N, N) ): The surface hopping matrix, the element g(i,j) 
            contains the probability for a i->j transition. Here, N - is the 
            total number of states considered in the transitions
        ksi ( double ): A random number uniformly distributed in the range of (0.0, 1.0) 
            Essentially, it determines the outcome of the procedure

    Returns: 
        int: finstate: The index of the final state after hop

    """

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
    """
    This function implements a simple random state selection procedure. 
    Each state is selected with a given probability

    Args:
        prob ( list of N doubles ): The probabilities of all N states 
        ksi ( double ): A random number uniformly distributed in the range of (0.0, 1.0).
            It determines the outcome of this function.

    Returns:
        integer: finstate: The index of the selected state

    """

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

def surface_hopping_nbra(mol, el, ham, rnd, params):
    ## This function performs surface hopping under neglect of back-reaction approximation.
    # Here, the number of nuclear objects and that of electronic ones are not equal.
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

    ninit = params["nconfig"]
    nstates_init = len(params["excitations_init"])
    nstates = el[0].nstates  # how many electronic DOF
    num_SH_traj = params["num_SH_traj"]

    g = MATRIX(nstates,nstates) # initialize a matrix of hopping probability

    for iconf in xrange(ninit):
        for i_ex in xrange(nstates_init):
            i = iconf*nstates_init + i_ex
            for itraj in xrange(num_SH_traj): # all stochastic SH realizations
                iel = iconf*nstates_init*num_SH_traj + i_ex*num_SH_traj + itraj
                #Compute hopping probabilities
                if tsh_method == 1: # FSSH
                    compute_hopping_probabilities_fssh(mol[i], el[iel], ham[i], g, dt_nucl, use_boltz_factor, Temperature)
                elif tsh_method == 2: # GFSH
                    compute_hopping_probabilities_gfsh(mol[i], el[iel], ham[i], g, dt_nucl, use_boltz_factor, Temperature)
                elif tsh_method == 3: # MSSH
                    compute_hopping_probabilities_mssh(mol[i], el[iel], ham[i], g, dt_nucl, use_boltz_factor, Temperature)
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
                el[iel].istate = hop(el[iel].istate, mol[i], ham[i], ksi, g, do_rescaling, rep, do_reverse)

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
    surface_hopping_nbra(mol, el, ham, rnd, params)



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




def boltz_factor(E_new, E_old, T, boltz_opt):
    # Compute the Boltzmann scaling factor, but only if we consider a hop up in energy

    dE = (E_new - E_old)
    boltz_f = 1.0 

    if boltz_opt==0:
        boltz_f = 1.0

    elif boltz_opt==1:
        if dE > 0.0:
            argg = dE/(units.kB*T)
            if argg > 50.0:
                boltz_f = 0.0
            else:
                boltz_f = math.exp(-argg)

    elif boltz_opt==2:
        if dE > 0.0:
            boltz_f = probabilities.Boltz_cl_prob_up(dE, T)

    elif boltz_opt==3:
        if dE > 0.0:
            boltz_f = probabilities.Boltz_quant_prob([0.0, dE], T)[1]

    return boltz_f



def ida_py(Coeff, old_st, new_st, E_old, E_new, T, ksi, do_collapse, boltz_opt=1):

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
    # \param[in] boltz_opt [0, 1, or 2] How to determine if the hop may be frustrated:
    #              0 - all proposed hops are accepted - no rejection based on energies
    #              1 - proposed hops are accepted with exp(-E/kT) probability - the old (hence the default approach)
    #              2 - proposed hops are accepted with the probability derived from Maxwell-Boltzmann distribution - more rigorous
    #              3 - generalization of "1", but actually it should be changed in case there are many degenerate levels

    # The function returns:
    # res [ integer ] - index of the final state, after IDA is applied (or not)
    # C [CMATRIX or Electronic] - the updated state of the electronic DOF, in the same data type as the input

    res = old_st

    dE = E_new - E_old
    boltz_f = boltz_factor(E_new, E_old, T, boltz_opt)

    # In case the electronic DOF are given in the form of CMATRIX
    if type(Coeff).__name__ == "CMATRIX":

        C = CMATRIX(Coeff)
        
        if dE>0.0:
        
            if ksi<boltz_f:
                res = new_st  # we've got enough kinetic energy - accept the hop
                
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
                if do_collapse:
                    for st in xrange(C.nstates):
                        C.q[st], C.p[st] = 0.0, 0.0
                    C.q[new_st], C.p[new_st] = 1.0, 0.0
            else:
                # Unsuccessful hop - collapse wfc back to the original state^M                                                                                         
                if do_collapse:
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

        



def hopping(Coeff, Hvib, istate, sh_method, do_collapse, ksi, ksi2, dt, T, boltz_opt=1):
    """
    A simplified version for the CPA-like hopping

    Coeff (CMATRIX(nstates, 1) ) object with the amplitudes of all states
    Hvib (CMATRIX(nstates, nstates) )  object containing the vibronic Hamiltonian 
    istate (int) the index of the initial state
    sh_method (int) selector of the TSH method: 0 - MSSH, 1 - FSSH, 2 - GFSH
    do_collapse (int) flag to turn on the decoherence via ID-A: 0 - no decoherence, 1 - decoherence via ID-A
    ksi, ksi2 (float in [0, 1]) random numbers cotrolling the execution of SH
    dt (float) time interval for the surface hopping (in a.u.)
    T (float) temperature in K

    boltz_opt [0, 1, 2, or 3] How to determine if the hop may be frustrated:
               0 - all proposed hops are accepted - no rejection based on energies
               1 - proposed hops are accepted with exp(-E/kT) probability - the old (hence the default approach)
               2 - proposed hops are accepted with the probability derived from Maxwell-Boltzmann distribution - more rigorous
               3 - generalization of "1", but actually it should be changed in case there are many degenerate levels

    Returns: the index (int) of a new state 

    """
    g = 0.0
    if sh_method==0:
        g = compute_hopping_probabilities_mssh(Coeff)
    elif sh_method==1:
        g = compute_hopping_probabilities_fssh(Coeff, Hvib, dt)
    elif sh_method==2:
        g = compute_hopping_probabilities_gfsh(Coeff, Hvib, dt)


    old_st = istate
    new_st = hop(istate, g, ksi)

    Coeff1 = CMATRIX(Coeff)

    if new_st != old_st:
        E_old = Hvib.get(old_st,old_st).real
        E_new = Hvib.get(new_st,new_st).real

        # ID-A decoherence                
        istate, Coeff1 = ida_py(Coeff, old_st, new_st, E_old, E_new, T, ksi2, do_collapse, boltz_opt)


    return istate, Coeff1
    


def project_out(Coeff, i):
    """
      Projects the state i out of a coherent superposition of states
      i   [int]  The index of the state to be projected out

    """
    nstates = Coeff.num_of_rows
    
    ci = Coeff.get(i,0)
    pi = (ci.conjugate() * ci).real
    nrm = 1.0 - pi
    if nrm<0.0:
        nrm = 0.0
    if nrm>0.0:
        nrm = 1.0/math.sqrt(nrm)

    Coeff.scale(-1, 0, nrm)
    Coeff.set(i, 0, 0.0+0.0j)



def collapse(Coeff, i):   
    """ 
    Collapse the wfc but such that to preserve the phase!
    """
    Coeff *= 0.0

    ci = Coeff.get(i,0)
    pi = (ci.conjugate() * ci).real
    if pi>0.0:
        Coeff.set(i, 0, ci/math.sqrt(pi))     
    else:
        Coeff.set(i, 0, 1.0+0.0j)


def dish_py(Coeff, istate, t_m, tau_m, Hvib, boltz_opt, T, ksi1, ksi2):
    """
    Decoherence-induced surface hopping (DISH)
    Reference: Jaeger, H. M.; Fischer, S.; Prezhdo, O. V. Decoherence-Induced Surface Hopping. J. Chem. Phys. 2012, 137, 22A545.

    Coeff  [CMATRIX, N x 1]   Amplitudes of electronic states
    istate [int]          Initial state
    t_m    [MATRIX, N x 1]    Matrix of the times each state resides in a coherence interval since the last decoherence event
    tau_m  [MATRIX, N x 1]    Matrix of the coherence intervals for each electronic state
    Hvib   [CMATRIx, N x N]   Energies 
    boltz_opt [0, 1, 2, or 3] How to determine if the hop may be frustrated:
               0 - all proposed hops are accepted - no rejection based on energies
               1 - proposed hops are accepted with exp(-E/kT) probability - the old (hence the default approach)
               2 - proposed hops are accepted with the probability derived from Maxwell-Boltzmann distribution - more rigorous
               3 - generalization of "1", but actually it should be changed in case there are many degenerate levels
    T      [double]  The temperature of the system [K]
    ksi1   [double]  A random number from a uniform distribution on the [0,1] interval
    ksi2   [double]  Another random number from a uniform distribution on the [0,1] interval

    The function modifies  Coeff and t_m variables
    Returns: 
     - the index of the electronic state after the "hop" (old state index or a new one)

    """

    fstate = istate

    nstates = Coeff.num_of_rows
    dm = Coeff * Coeff.H()     # density matrix

    i = 0
    has_decoherence = False   # set to True if we have encountered a decoherence event

    while i<nstates and not has_decoherence:

        # The state i has evolved coherently for longer than the coherence interval
        # so it has to experience a decoherence event 
        if t_m.get(i) >= tau_m.get(i):
            # There are essentially two outcomes when the decoherence takes place:

            if ksi1 < dm.get(i,i).real:  
                # One: we collapse the wavefunction onto the state i with the probability 
                # given by the population of that state

                # Now, lets determine if the hop is possible based on the energy conservation
                # considerations 
                boltz_f = boltz_factor(Hvib.get(i,i).real, Hvib.get(istate, istate).real, T, boltz_opt)
 
                if ksi2 < boltz_f:
                    #  Now, decide about decoherence         
                    collapse(Coeff, i)      # Here is the actuall collapse
                    fstate = i
                else:
                    project_out(Coeff, i)       # Project out of the state
                    fstate = istate

            # The below section is not present in the DISH implementation in PYXAID
            # but this is not how it should be done according to the paper of Jaeger et al.
            # so this implementation does follow the algorithm outlined in the paper
            else:
                # Second: project the system out of that state
                project_out(Coeff, i)   # Project out of the state
                fstate = istate

            # Reset the time axis for state i (only for this state)
            # other states still reside in a coherent superposition
            t_m.set(i, 0.0)

            # Set the flag that we have attempted a decoherence event
            # so we done with DISH at this point in time
            has_decoherence = True

        i = i + 1

    return fstate;

