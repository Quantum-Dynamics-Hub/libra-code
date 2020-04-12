#*********************************************************************************                     
#* Copyright (C) 2019 Alexey V. Akimov                                                   
#*                                                                                                     
#* This file is distributed under the terms of the GNU General Public License                          
#* as published by the Free Software Foundation, either version 3 of                                   
#* the License, or (at your option) any later version.                                                 
#* See the file LICENSE in the root directory of this distribution   
#* or <http://www.gnu.org/licenses/>.          
#***********************************************************************************
"""
.. module:: dynamics
   :platform: Unix, Windows
   :synopsis: This module implements a wrapper function for doing Ehrenfest/TSH/Verlet/etc. dynamics
       List of functions:  
           * init_nuclear_dyn_var(Q, P, M, params, rnd)
           * init_electronic_dyn_var(params, rnd)
           * init_amplitudes(q, Cdia, Cadi, dyn_params, compute_model, model_params, transform_direction=0)
           * run_dynamics(_q, _p, _iM, _Cdia, _Cadi, _states, _dyn_params, compute_model, _model_params, rnd)
           * generic_recipe(q, p, iM, _dyn_params, compute_model, _model_params, _init_elec, rnd)
           * run_multiple_sets(init_cond, _dyn_params, compute_model, _model_params, _init_nucl, _init_elec, rnd)

.. moduleauthor:: Alexey V. Akimov

"""

__author__ = "Alexey V. Akimov"
__copyright__ = "Copyright 2019 Alexey V. Akimov"
__credits__ = ["Alexey V. Akimov"]
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

import util.libutil as comn
import libra_py.units as units
import libra_py.data_outs as data_outs
import libra_py.data_savers as data_savers
import libra_py.tsh as tsh
import libra_py.tsh_stat as tsh_stat
#import libra_py.dynamics as dynamics_io

from . import save


def init_nuclear_dyn_var(Q, P, M, params, rnd):
    """
    Args:
        Q ( list of doubles ): the mean values of coordinates for all DOFs [ units: a.u.]
        P ( list of doubles ): the mean values of momenta for all DOFs [ units: a.u. ]
        M ( list of doubles ): masses of all nuclear DOFs [ units: a.u. ]

        params ( dictionary ): control parameters
 
            * **params["init_type"]** ( int ): the type of sampling of nuclear DOFs
     
                - 0 : initialize ```ntraj``` identical copies of coordinates and momenta
 
                - 1 : keep all coordinates identical, but sample momenta from the normal 
                    distribution with a given width in each dimension

                - 2 : keep all momenta identical, but sample coordinates from the normal 
                    distribution with a given width in each dimension

                - 3 : sample both coordinates and momenta from the normal 
                    distributions with given widths in each dimension

            * **params["force_constant"]** ( list of double ): force constants involved in the Harmonic
                oscillator model: U = (1/2) * k * x^2, and omega = sqrt( k / m )
                These parameters define the Harmonic oscillator ground state wavefunctions from which
                the coordinates and momenta are sampled. [ units: a.u. = Ha/Bohr^2, default: [0.001] ]
                The length should be consistent with the length of Q, P and M arguments

            * **params["ntraj"]** ( int ): the number of trajectories - the parameter defines the
                number of columns the output matrices will have [ default: 1 ]



        rnd ( Random ): random numbers generator object


    Returns:
        q, p, iM:  where:

            * q ( MATRIX(ndof, ntraj) ) : coordinates for all trajectories
            * p ( MATRIX(ndof, ntraj) ) : momenta for all trajectories
            * iM ( MATRIX(ndof, 1) ) : inverse masses of all DOFs (same across the trajectories)

    """


    critical_params = [ ]
    default_params = { "init_type":0, "force_constant":[0.001], "ntraj":1  }
    comn.check_input(params, default_params, critical_params)
        
    init_type = params["init_type"]  
    force_constant = params["force_constant"]
    ntraj = params["ntraj"]    


    if init_type not in [0, 1, 2, 3]:
        print(F"WARNINIG in init_nuclear_dyn_var: \
              the init_type = {init_type} is not known\
              Allowed values are: [0, 1, 2, 3]" )

    if len(Q)!=len(P):
        print(F"ERROR in init_nuclear_dyn_var: \
              the length of input Q is = {len(Q)}, \
              the length of input P is = {len(P)}, but they should be equal to each other" )
        sys.exit(0)

    if len(Q)!=len(M):
        print(F"ERROR in init_nuclear_dyn_var: \
              the length of input Q is = {len(Q)}, \
              the length of input M is = {len(M)}, but they should be equal to each other" )
        sys.exit(0)

    if len(Q)!=len(force_constant):
        print(F"ERROR in init_nuclear_dyn_var: \
              the length of input Q is = {len(Q)}, \
              the length of input P is = {len(force_constant)}, but they should be equal to each other" )
        sys.exit(0)

    if ntraj < 1:
        print(F"ERROR in init_nuclear_dyn_var: \
              the ntraj is= {ntraj}, should be at least 1" )
        sys.exit(0)



    # At this point, it is safe to define ndof:
    ndof = len(Q)
    q, p, iM = MATRIX(ndof,ntraj), MATRIX(ndof,ntraj), MATRIX(ndof,1)
    for dof in range(ndof):
        iM.set(dof, 0, 1.0/M[dof])


    # Mean values
    mean_q, mean_p = MATRIX(ndof,1), MATRIX(ndof, 1)
    for dof in range(ndof):
        mean_q.set(dof, 0, Q[dof])
        mean_p.set(dof, 0, P[dof])

    # Deviations
    sigma_q, sigma_p = MATRIX(ndof,1), MATRIX(ndof, 1)
    if init_type == 0:
        for dof in range(ndof):        
            sigma_q.set(dof, 0, 0.0)
            sigma_p.set(dof, 0, 0.0)

    elif init_type == 1:
        for dof in range(ndof):        
            sigma_q.set(dof, 0, 0.0)
            sigma_p.set(dof, 0, math.sqrt( 0.5*math.sqrt((force_constant[dof]*M[dof])) ))

    elif init_type == 2:
        for dof in range(ndof):        
            sigma_q.set(dof, 0, math.sqrt( 0.5*math.sqrt(1.0/(force_constant[dof]*M[dof])) ))
            sigma_p.set(dof, 0, 0.0 )

    elif init_type == 3:
        for dof in range(ndof):        
            sigma_q.set(dof, 0, math.sqrt( 0.5*math.sqrt(1.0/(force_constant[dof]*M[dof])) ))
            sigma_p.set(dof, 0, math.sqrt( 0.5*math.sqrt((force_constant[dof]*M[dof])) ))

    # Now sample the values
    tsh.sample(q, mean_q, sigma_q, rnd)
    tsh.sample(p, mean_p, sigma_p, rnd)

    return q, p, iM



def init_electronic_dyn_var(params, rnd):
    """
    Args:

        params ( dictionary ): control parameters
 
            * **params["init_type"]** ( int ): the type of sampling of electronic DOFs
     
                - 0 : initialize all states according to "istate" and sets all
                    amplitudes to 1.0  [ default ]
 
                - 1 : initialize all states according to "istate" but sets 
                    amplitudes to exp(-2*pi*i*rnd), where rnd is a random 
                    number uniformly distributed on the [0, 1] interval

                - 2 : initialize all states according to "istates" - the 
                    integer indices are selected randomly according to populations provided in
                    variable "istates", the amplitudes are set to be sqrt(istates[i]), but
                    their phases are identical (like in the option 0)

                - 3 : initialize all states according to "istates" - the 
                    integer indices are selected randomly according to populations provided in
                    variable "istates", the amplitudes are set to be sqrt(istates[i]), but
                    their phases are set to exp(-2*pi*i*rnd), where rnd is a random 
                    number uniformly distributed on the [0, 1] interval (line in option 1)

            * **params["nstates"]** ( int ): the number of electronic states in the basis
                [ default: 1 ]

            * **params["istate"]** ( int ): the index of the initial electronic state, used 
                only when **params["init_type"]** is 0 or 1, in which case it defines on 
                which state the amplitudes will be initialized [ default: 0]

            * **params["istates"]** ( list of ints ): the list of the populations on all electronic
                states included in the basis, used only when **params["init_type"]** is 2 or 3, 
                in which case it defines on which states the amplitudes will be initialized.
                The length of this list should be consistent with ```nstates``` variable. And the sum
                of all entries should be 1.0 [ default: [1.0], meaning that only the lowest state
                is occupied ]

            * **params["rep"]** ( int ): defines for which repersentation we generate the amplitudes
 
                - 0 : diabatic, the corresponding matrix is initialized, the adiabatic amplitudes are zeroes
                - 1 : adiabatic, the corresponding matrix is initialized, the diabatic amplitudes are zeroes [ default ]

            * **params["ntraj"]** ( int ): the number of trajectories - the parameter defines the
                number of columns the output matrices will have [ default: 1]

        rnd ( Random ): random numbers generator object


    Returns:
        Cdia, Cadi, states:  where:

            * Cdia ( CMATRIX(nstates, ntraj) ) : amplitudes on all diabatic states for all trajectories
            * Cadi ( CMATRIX(nstates, ntraj) ) : amplitudes on all adiabatic states for all trajectories
            * states ( list of ints ) : state indices for each trajectory

    """

    # Read the parameters
    critical_params = [  ]
    default_params = { "init_type":0, "nstates":1, "istate":0, "istates":[1.0], "rep":1,  "ntraj":1  }
    comn.check_input(params, default_params, critical_params)

    init_type = params["init_type"]
    nstates = params["nstates"]  
    istate = params["istate"]          
    istates = params["istates"]
    rep = params["rep"]  
    ntraj = params["ntraj"]

    # Sanity check
    if rep not in [0, 1]:
        print(F"WARNINIG in init_electronic_dyn_var: \
              the rep = {rep} is not known\
              Allowed values are: [0, 1]" )

    if init_type not in [0, 1, 2, 3]:
        print(F"WARNINIG in init_electronic_dyn_var: \
              the init_type = {init_type} is not known\
              Allowed values are: [0, 1, 2, 3]" )

    if ntraj < 1:
        print(F"ERROR in init_electronic_dyn_var: \
              the ntraj is= {ntraj}, should be at least 1" )
        sys.exit(0)
    
    if init_type in [0, 1]:
        if istate >= nstates:
            print(F"ERROR in init_electronic_dyn_var: \
                  the istate is= {istate}, but should be less than {nstates}" )
            sys.exit(0)
             
    if init_type in [2, 3]:
        if len(istates)!= nstates:
            print(F"ERROR in init_electronic_dyn_var: \
                  the istates array is of length {len(istates)}, but should \
                  be of length {nstates}")
            sys.exit(0)
        
        summ = sum(istates)
        if math.fabs( summ - 1.0 ) > 1e-5:
            print(F"ERROR in init_electronic_dyn_var: \
                  the sum of the entries in the istates array is {summ}, but should be 1.0")
            sys.exit(0)

               
                
    # Dynamical variables
    Cdia, Cadi = CMATRIX(nstates, ntraj), CMATRIX(nstates, ntraj)
    states = intList() 


    for traj in range(ntraj):

        if init_type==0:
 
            if rep==0:
                Cdia.set(istate, traj, 1.0+0.0j);  
            elif rep==1:
                Cadi.set(istate, traj, 1.0+0.0j);  
            states.append(istate) 

        elif init_type==1:

            ksi = rnd.uniform(0.0, 1.0)
            ampl = math.cos(2*math.pi*ksi) + 1.0j*math.sin(2.0*math.pi*ksi)

            if rep==0:
                Cdia.set(istate, traj, ampl);  
            elif rep==1:
                Cadi.set(istate, traj, ampl);  
            states.append(istate) 


        elif init_type==2:

            for state, pop in enumerate(istates):
                ampl = math.sqrt( pop )

                if rep==0:
                    Cdia.set(state, traj, ampl);  
                elif rep==1:
                    Cadi.set(state, traj, ampl);  

            ksi = rnd.uniform(0.0, 1.0)
            states.append(tsh.set_random_state(istates, ksi)) 


        elif init_type==3:

            for state, pop in enumerate(istates):

                ksi = rnd.uniform(0.0, 1.0)
                ampl = math.cos(2*math.pi*ksi) + 1.0j*math.sin(2.0*math.pi*ksi)
                ampl = ampl * math.sqrt( pop )

                if rep==0:
                    Cdia.set(state, traj, ampl);  
                elif rep==1:
                    Cadi.set(state, traj, ampl);  

            ksi = rnd.uniform(0.0, 1.0)
            states.append(tsh.set_random_state(istates, ksi)) 

    projections = CMATRIXList()
    for traj in range(ntraj):
        projections.append( CMATRIX(nstates, nstates) )
        projections[traj].identity()


    return Cdia, Cadi, projections, states



def init_amplitudes(q, Cdia, Cadi, dyn_params, compute_model, model_params, transform_direction=0):
    """
    Args:
        q ( MATRIX(ndof, ntraj) ): coordinates 
        Cdia ( MATRIX(ndia, ntraj) ): amplitudes of diabatic states
        Cadi ( MATRIX(nadi, ntraj) ): amplitudes of adiabatic states

        dyn_params ( Python dictionary ): control of the dynamics
        compute_model ( Python function ): the function that does the calculations
        model_params ( Python dictionary ): parameters of the computational model

        transform_direction ( int ): type of transformation

            - 0: diabatic to adiabatic
            - 1: adiabatic to diabatic


    Returns:
        None: but changes Cadi or Cdia, according to the transform_direction

    """

    # Dimensions
    ndia = Cdia.num_of_rows
    nadi = Cadi.num_of_rows
    nnucl= q.num_of_rows
    ntraj= q.num_of_cols


    # Prepare the Hamiltonian's hierarchy
    ham = nHamiltonian(ndia, nadi, nnucl)
    ham.add_new_children(ndia, nadi, nnucl, ntraj)
    ham.init_all(2,1)

    # Compute the Hamiltonian transformation
    update_Hamiltonian_q(dyn_params, q, ham, compute_model, model_params)

    # Do the transformations
    if transform_direction==0:  # dia -> adi
        Cadi = transform_amplitudes(dyn_params, Cdia, ham)
    elif rep_tdse==1:  # adi -> dia
        Cdia = transform_amplitudes(dyn_params, Cadi, ham)



def run_dynamics(_q, _p, _iM, _Cdia, _Cadi, _projectors, _states, _dyn_params, compute_model, _model_params, rnd):
    """
    
    Args: 
        _q ( MATRIX(nnucl, ntraj) ): coordinates of the "classical" particles [units: Bohr]
        _p ( MATRIX(nnucl, ntraj) ): momenta of the "classical" particles [units: a.u. of momenta]
        _iM ( MATRIX(nnucl, 1) ): masses of classical particles [units: a.u.^-1]
        _Cdia ( CMATRIX(ndia, ntraj) ): amplitudes of the diabatic basis states
        _Cadi ( CMATRIX(nadi, ntraj) ): amplitudes of the adiabatic basis states
        _projectors ( list of CMATRIX(nadi, nadi) ): cumulative phase correction and state tracking matrices
        _states ( intList, or list of ntraj ints ): the quantum state of each trajectory

        dyn_params ( dictionary ): parameters controlling the execution of the dynamics
            Can contain:
      
            * **dyn_params["rep_tdse"]** ( int ): selects the representation in which 
                nuclear/electronic (Ehrenfest core) dynamics is executed. The representation 
                used to integrate TD-SE

                - 0: diabatic representation
                - 1: adiabatic representation [ default ]


            * **dyn_params["rep_ham"]** (int): The representation of the Hamiltonian update: 

                - 0: diabatic   [ default ]
                - 1: adiabatic
                              
                This is the representation in which the computed properties are assumed to be
                For instance, we may have set it to 1, to read the adiabatic energies and couplings,
                to bypass the diabatic-to-adiabatic transformation, which may be useful in some atomistic
                calculations, or with the NBRA


            * **dyn_params["rep_sh"]** ( int ): selects the representation which is 
                used to perform surface hopping

                - 0: diabatic representation
                - 1: adiabatic representation [ default ]


            * **dyn_params["rep_lz"]** ( int ): The representation to compute LZ probabilitieis:
 
                - 0: diabatic [ default ] 
                - 1: adiabatic 


            * **dyn_params["tsh_method"]** ( int ): Formula for computing SH probabilities: 
   
                - -1: adiabatic - no hops [ default ]
                - 0: FSSH
                - 1: GFSH 
                - 2: MSSH
                - 3: DISH


            * **dyn_params["force_method"]** ( int ): How to compute forces in the dynamics: 

                - 0: don't compute forces at all - e.g. we do not really need them
                - 1: state-specific  as in the TSH or adiabatic (including adiabatic excited states) [ default ]
                - 2: Ehrenfest


            * **dyn_params["nac_update_method"]** ( int ): How to update NACs and vibronic Hamiltonian before
               electronic TD-SE propagation:

                - 0: don't update them (e.g. for simplest NAC)
                - 1: update according to changed momentum and existing derivative couplings [ default ]


            * **dyn_params["rep_force"]** ( int ): In which representation to compute forces: 

                - 0: diabatic
                - 1: adiabatic [ default ]


            * **dyn_params["use_boltz_factor"]** ( int ): Whether to scale the SH probabilities
                by the Boltzmann factor: 
 
                - 0: do not scale [ default ]
                - 1: scale 


            * **dyn_params["Temperature"]** ( double ): Temperature of the system [ default: 300 K ]


            * **dyn_params["do_reverse"]** ( int ): what to do with velocities on frustrated hops:

                - 0: do not revert momenta at the frustrated hops
                - 1: do revert the momenta [ default ]

            * **dyn_params["vel_rescale_opt"]** ( int ): How to rescale momenta if the hops are successful:

                - -1: do not rescale, as in the NBRA [ default ]
                - 0: rescale in the diabatic basis - don't care about the
                    velocity directions, just a uniform rescaling 
                - 1: rescale along the directions of derivative couplings


            * **dyn_params["do_phase_correction"]** ( int ): the algorithm to correct phases on adiabatic states
 
                - 0: no phase correction
                - 1: according to our phase correction algorithm [ default: 1 ]


            * **dyn_params["tol"]** ( double ): the minimal value of the time-overlap for considering 
                phase corrections - no correction applied if the time-overlap is smaller in magnitude than 
                this parameter [ default: 1e-3 ]
 

            * **dyn_params["state_tracking_algo"]** ( int ): the algorithm to keep track of the states' identities
 
                - 0: no state tracking
                - 1: Sato
                - 2: using the mincost, Munkres-Kuhn [ default: 2 ]


            * **dyn_params["MK_alpha"]** ( double ): Munkres-Kuhn alpha 
                (selects the range of orbitals included in reordering) [default: 0.0]


            * **dyn_params["MK_verbosity"]** ( double ): Munkres-Kuhn verbosity: 

                - 0: no extra output [ default ]
                - 1: print details



            * **dyn_params["entanglement_opt"]** ( int ): A selector of a method to couple the trajectories in this ensemble:

                - 0: no coupling across trajectories [ default ]
                - 1: ETHD
                - 2: ETHD3 (experimental)
                - 22: another flavor of ETHD3 (experimental)

            * **dyn_params["ETHD3_alpha"]** ( double ): Gaussian exponents that dresses up the trajectories in the ETHD3 method
                in the coordinate space, that is   ~exp(-alpha*(R-R0)^2 ) [ default: 0.0 ]


            * **dyn_params["ETHD3_beta"]** ( double ): Gaussian exponents that dresses up the trajectories in the ETHD3 method
                in the coordinate space, that is   ~exp(-alpha*(P-P0)^2 ) [ default: 0.0 ]


            * **dyn_params["decoherence_algo"]** ( int ): selector of the method to incorporate decoherence:

                - -1: no decoherence [ default ]
                - 0: MSDM
                - 1: ID-A


            * **dyn_params["decoh_rates"]** ( MATRIX(ntates, nstates) ): the matrix of dephasing rates [ units: a.u. of time ^-1 ]


            * **dyn_params["dt"]** ( double ): the nuclear and electronic integration
                timestep [ units: a.u. of time, default: 41.0 ]


            * **dyn_params["nsteps"]** ( int ): the number of NA-MD steps to do [ default: 1 ]


            * **dyn_params["prefix"]** ( string ): the name of the folder to be created by this function. 
                All the data files will be created in that folder


            * **dyn_params["hdf5_output_level"]** ( int ): controls what info to return as the result of this function

                - -1: all return values are None [ default ]
                - 0: also, all return values are none 
                - 1: 1-D info - energies, SE and SH populations averaged over ensemble vs. time
                - 2: 2-D info - coordinates, momenta, SE amplitudes, and SH states populations 
                    for each individual trajectory vs. time 
                - 3: 3-D info - St, Hvib_adi, Hvib_dia matrices (energies, couplings, etc.) for each
                    individual trajectory vs. time 

            * **dyn_params["mem_output_level"]** ( int ): controls what info to return as the result of this function

                Same meaning and output as with hdf5_output_level, except all the variables are first stored in memory
                (while the calcs are running) and then they are written into the HDF5 file at the end of the calculations.
                This is a much faster version of hdf5 saver. 


        compute_model ( PyObject ): the pointer to the Python function that performs the Hamiltonian calculations

        _model_params ( dictionary ): contains the selection of a model and the parameters 
            for that model Hamiltonian

        rnd ( Random ): random numbers generator object



    Returns:
        tuple: with the elements of the tuple listed below.

            * obs_T ( list of `nsteps` doubles ): time [units: a.u.]
            * obs_q ( list of `nsteps` MATRIX(nnucl, ntraj) ): coordinates of all trajectories [ units: Bohr ]
            * obs_p ( list of `nsteps` MATRIX(nnucl, ntraj) ): momenta of all trajectories [ units: a.u. ]
            * obs_Ekin ( list of `nsteps` doubles ): average kinetic energy of an ensemble of trajectories [units: a.u.]
            * obs_Epot ( list of `nsteps` doubles ): average potential energy of an ensemble of trajectories [units: a.u.]
            * obs_Etot ( list of `nsteps` doubles ): average total energy of an ensemble of trajectories [units: a.u.]
            * obs_dEkin ( list of `nsteps` doubles ): standard deviation of kinetic energy of an ensemble of trajectories [units: a.u.]
            * obs_dEpot ( list of `nsteps` doubles ): standard deviation of potential energy of an ensemble of trajectories [units: a.u.]
            * obs_dEtot ( list of `nsteps` doubles ): standard deviation of total energy of an ensemble of trajectories [units: a.u.]
            * obs_Cadi ( list of `nsteps` CMATRIX(nadi, ntraj) ): amplitudes of adiabatic electronic states of all trajectories 
            * obs_Cdia ( list of `nsteps` CMATRIX(ndia, ntraj) ): amplitudes of diabatic electronic states of all trajectories 
            * obs_dm_adi ( list of `nsteps` CMATRIX(nadi, nadi) ): ensemble-averaged density matrix in adiabatic basis
            * obs_dm_dia ( list of `nsteps` CMATRIX(ndia, ndia) ): ensemble-averaged density matrix in diabatic basis
            * obs_pop ( list of `nsteps` MATRIX(nadi, 1) ): ensemble-averaged TSH populations of adiabatic states
            * obs_states ( list of `nsteps` of lists of `ntraj` ints):  # indices of the quantum states of each trajectory
            * obs_hvib_adi ( list of `ntraj` lists, each being a list of `nsteps` objects CMATRIX(nadi, nadi) ): trajectory-resolved
                vibronic Hamiltonians for each timestep in the adiabatic representation
            * obs_hvib_dia ( list of `ntraj` lists, each being a list of `nsteps` objects CMATRIX(ndia, ndia) ): trajectory-resolved
                vibronic Hamiltonians for each timestep in the diabatic representation
            * obs_St ( list of `ntraj` lists, each being a list of `nsteps` objects CMATRIX(nadi, nadi) ): trajectory-resolved
                time-overlaps of the adiabatic wavefunctions
            * obs_U ( list of `ntraj` lists, each being a list of `nsteps` objects CMATRIX(ndia, nadi) ): trajectory-resolved
                diabatic-to-adiabatic transformation matrix


        .. note::
          the elements are None, if they are excluded by the input options

    """

        
    # Create copies of the input dynamical variables, so we could run several such
    # functions with the same input variables without worries that they will be altered
    # inside of each other

    model_params = dict(_model_params)
    dyn_params = dict(_dyn_params)

    q = MATRIX(_q)
    p = MATRIX(_p)
    iM = MATRIX(_iM)
    Cdia = CMATRIX(_Cdia)
    Cadi = CMATRIX(_Cadi)
    states = intList()
    projectors = CMATRIXList()
    
    
    for i in range(len(_states)):
        states.append(_states[i])
        projectors.append(CMATRIX(_projectors[i]))

    DR = MATRIX(len(_states), len(_states))
    AG = MATRIX(len(_states), len(_states))


    # Parameters and dimensions
    properties_to_save = [ "timestep", "time", "Ekin_ave", "Epot_ave", "Etot_ave", 
                           "dEkin_ave", "dEpot_ave", "dEtot_ave", "states", "SH_pop", "SH_pop_raw",
                           "D_adi", "D_adi_raw", "D_dia", "D_dia_raw", "q", "p", "Cadi", "Cdia", 
                           "hvib_adi", "hvib_dia", "St", "basis_transform", "projector"
                         ] 

    critical_params = [  ] 
    default_params = { "rep_tdse":1, "rep_ham":0, "rep_sh":1, "rep_lz":0, "tsh_method":-1,
                       "force_method":1, "nac_update_method":1, "rep_force":1,
                       "use_boltz_factor":0, "Temperature":300.0, "do_reverse":1, "vel_rescale_opt":-1,
                       "do_phase_correction":1, "tol":1e-3,
                       "state_tracking_algo":2, "MK_alpha":0.0, "MK_verbosity":0, 
                       "entanglement_opt":0, "ETHD3_alpha":0.0, "ETHD3_beta":0.0, 
                       "decoherence_algo":-1, "decoherence_rates":DR,
                       "decoherence_times_type":0, "decoherence_C_param":1.0, 
                       "decoherence_eps_param":0.1, "dephasing_informed":0,
                       "ave_gaps":AG, "instantaneous_decoherence_variant":1, "collapse_option":0,
                       "ensemble":0, "thermostat_params":{},
                       "dt":1.0*units.fs2au, "nsteps":1, 
                       "hdf5_output_level":-1, "prefix":"out", "properties_to_save":properties_to_save,
                       "mem_output_level":-1,
                       "use_compression":0, "compression_level":[0,0,0], 
                       "txt_output_level":0,
                       "progress_frequency":0.1
                     }

    comn.check_input(dyn_params, default_params, critical_params)
        
    prefix = dyn_params["prefix"]    
    rep_tdse = dyn_params["rep_tdse"]
    nsteps = dyn_params["nsteps"]
    dt = dyn_params["dt"]
    tol = dyn_params["tol"]
    hdf5_output_level = dyn_params["hdf5_output_level"]
    mem_output_level = dyn_params["mem_output_level"]
    do_phase_correction = dyn_params["do_phase_correction"]
    state_tracking_algo = dyn_params["state_tracking_algo"]
    force_method = dyn_params["force_method"]
    properties_to_save = dyn_params["properties_to_save"]
    use_compression = dyn_params["use_compression"]
    compression_level = dyn_params["compression_level"]
    ensemble = dyn_params["ensemble"]
    
    ndia = Cdia.num_of_rows
    nadi = Cadi.num_of_rows
    nnucl= q.num_of_rows
    ntraj= q.num_of_cols


    # Initialize savers
    _savers = save.init_tsh_savers(dyn_params, model_params, nsteps, ntraj, nnucl, nadi, ndia)


    # ======= Hierarchy of Hamiltonians =======
    ham = nHamiltonian(ndia, nadi, nnucl)
    ham.add_new_children(ndia, nadi, nnucl, ntraj)
    ham.init_all(2,1)
    model_params.update({"timestep":0})
    
    update_Hamiltonian_q(dyn_params, q, projectors, ham, compute_model, model_params)
    update_Hamiltonian_p(dyn_params, ham, p, iM)  


    U = []
    for tr in range(ntraj):
        U.append(ham.get_basis_transform(Py2Cpp_int([0, tr]) ))


    therm = ThermostatList();
    if ensemble==1:
        for traj in range(ntraj):
            therm.append( Thermostat( dyn_params["thermostat_params"] ) )
            therm[traj].init_nhc();

                
    # Do the propagation
    for i in range(nsteps):
    
        #============ Compute and output properties ===========        
        # Amplitudes, Density matrix, and Populations
        if rep_tdse==0:
            # Diabatic to raw adiabatic
            ham.ampl_dia2adi(Cdia, Cadi, 0, 1) 
            # Raw adiabatic to dynamically-consistent adiabatic
            Cadi = dynconsyst_to_raw(Cadi, projectors)

        elif rep_tdse==1:
            ham.ampl_adi2dia(Cdia, Cadi, 0, 1)

        dm_dia, dm_adi, dm_dia_raw, dm_adi_raw, = tsh_stat.compute_dm(ham, Cdia, Cadi, projectors, rep_tdse, 1)        
        pops, pops_raw = tsh_stat.compute_sh_statistics(nadi, states, projectors)


        # Energies 
        Ekin, Epot, Etot, dEkin, dEpot, dEtot = 0.0, 0.0, 0.0,  0.0, 0.0, 0.0
        Etherm, E_NHC = 0.0, 0.0
        if force_method in [0, 1]:
            Ekin, Epot, Etot, dEkin, dEpot, dEtot = tsh_stat.compute_etot_tsh(ham, p, Cdia, Cadi, projectors, states, iM, rep_tdse)
        elif force_method in [2]:
            Ekin, Epot, Etot, dEkin, dEpot, dEtot = tsh_stat.compute_etot(ham, p, Cdia, Cadi, projectors, iM, rep_tdse)

        for bath in therm:
            Etherm += bath.energy()
        E_NHC = Etot + Etherm

        

        save.save_tsh_data_123(_savers, dyn_params, 
                       i, dt, Ekin, Epot, Etot, dEkin, dEpot, dEtot, Etherm, E_NHC, states,
                       pops, pops_raw, dm_adi, dm_adi_raw, dm_dia, dm_dia_raw, q, p, Cadi, Cdia  )

    
        for tr in range(ntraj):
            x = ham.get_basis_transform(Py2Cpp_int([0, tr]) )            
            St = U[tr] * x.H()             
            U[tr] = CMATRIX(x)

            if hdf5_output_level>=4: 
                hvib_adi = ham.get_hvib_adi(Py2Cpp_int([0, tr])) 
                hvib_dia = ham.get_hvib_dia(Py2Cpp_int([0, tr])) 
                save.save_hdf5_4D(_savers["hdf5_saver"], i, tr, hvib_adi, hvib_dia, St, U[tr], projectors[tr])

            if mem_output_level>=4: 
                hvib_adi = ham.get_hvib_adi(Py2Cpp_int([0, tr])) 
                hvib_dia = ham.get_hvib_dia(Py2Cpp_int([0, tr])) 
                save.save_hdf5_4D(_savers["mem_saver"], i, tr, hvib_adi, hvib_dia, St, U[tr], projectors[tr])



        #============ Propagate ===========        
        model_params.update({"timestep":i})        
        if rep_tdse==0:
            compute_dynamics(q, p, iM, Cdia, projectors, states, ham, compute_model, model_params, dyn_params, rnd, therm)
        elif rep_tdse==1:
            compute_dynamics(q, p, iM, Cadi, projectors, states, ham, compute_model, model_params, dyn_params, rnd, therm)


    if _savers["mem_saver"]!=None:
        _savers["mem_saver"].save_data( F"{prefix}/mem_data.hdf", properties_to_save, "w")
        return _savers["mem_saver"]




def run_dynamics_old(_q, _p, _iM, _Cdia, _Cadi, _projectors, _states, _dyn_params, compute_model, _model_params, rnd):
    """

    This is a deprecated function - it returns/saves the results of the calculations via either
    text files or memory, both ways are pretty bad. But if for some strange reason you can't use HDF5
    you may want to use this version. So, I keep int here, just in case.
  
    This function is similar to the one above, but it stores the computed properties in files
   
    Args: 
        _q ( MATRIX(nnucl, ntraj) ): coordinates of the "classical" particles [units: Bohr]
        _p ( MATRIX(nnucl, ntraj) ): momenta of the "classical" particles [units: a.u. of momenta]
        _iM ( MATRIX(nnucl, 1) ): masses of classical particles [units: a.u.^-1]
        _Cdia ( CMATRIX(ndia, ntraj) ): amplitudes of the diabatic basis states
        _Cadi ( CMATRIX(nadi, ntraj) ): amplitudes of the adiabatic basis states
        _projectors ( list of CMATRIX(nadi, nadi) ): cumulative phase correction and state tracking matrices
        _states ( intList, or list of ntraj ints ): the quantum state of each trajectory

        dyn_params ( dictionary ): parameters controlling the execution of the dynamics
            Can contain:
      
            * **dyn_params["rep_tdse"]** ( int ): selects the representation in which 
                nuclear/electronic (Ehrenfest core) dynamics is executed. The representation 
                used to integrate TD-SE

                - 0: diabatic representation
                - 1: adiabatic representation [ default ]


            * **dyn_params["rep_ham"]** (int): The representation of the Hamiltonian update: 

                - 0: diabatic   [ default ]
                - 1: adiabatic
                              
                This is the representation in which the computed properties are assumed to be
                For instance, we may have set it to 1, to read the adiabatic energies and couplings,
                to bypass the diabatic-to-adiabatic transformation, which may be useful in some atomistic
                calculations, or with the NBRA


            * **dyn_params["rep_sh"]** ( int ): selects the representation which is 
                used to perform surface hopping

                - 0: diabatic representation
                - 1: adiabatic representation [ default ]


            * **dyn_params["rep_lz"]** ( int ): The representation to compute LZ probabilitieis:
 
                - 0: diabatic [ default ] 
                - 1: adiabatic 


            * **dyn_params["tsh_method"]** ( int ): Formula for computing SH probabilities: 
   
                - -1: adiabatic - no hops [ default ]
                - 0: FSSH
                - 1: GFSH 
                - 2: MSSH
                - 3: DISH


            * **dyn_params["force_method"]** ( int ): How to compute forces in the dynamics: 

                - 0: don't compute forces at all - e.g. we do not really need them
                - 1: state-specific  as in the TSH or adiabatic (including adiabatic excited states) [ default ]
                - 2: Ehrenfest


            * **dyn_params["nac_update_method"]** ( int ): How to update NACs and vibronic Hamiltonian before
               electronic TD-SE propagation:

                - 0: don't update them (e.g. for simplest NAC)
                - 1: update according to changed momentum and existing derivative couplings [ default ]


            * **dyn_params["rep_force"]** ( int ): In which representation to compute forces: 

                - 0: diabatic
                - 1: adiabatic [ default ]


            * **dyn_params["use_boltz_factor"]** ( int ): Whether to scale the SH probabilities
                by the Boltzmann factor: 
 
                - 0: do not scale [ default ]
                - 1: scale 


            * **dyn_params["Temperature"]** ( double ): Temperature of the system [ default: 300 K ]


            * **dyn_params["do_reverse"]** ( int ): what to do with velocities on frustrated hops:

                - 0: do not revert momenta at the frustrated hops
                - 1: do revert the momenta [ default ]

            * **dyn_params["vel_rescale_opt"]** ( int ): How to rescale momenta if the hops are successful:

                - -1: do not rescale, as in the NBRA [ default ]
                - 0: rescale in the diabatic basis - don't care about the
                    velocity directions, just a uniform rescaling 
                - 1: rescale along the directions of derivative couplings


            * **dyn_params["do_phase_correction"]** ( int ): the algorithm to correct phases on adiabatic states
 
                - 0: no phase correction
                - 1: according to our phase correction algorithm [ default: 1 ]


            * **dyn_params["tol"]** ( double ): the minimal value of the time-overlap for considering 
                phase corrections - no correction applied if the time-overlap is smaller in magnitude than 
                this parameter [ default: 1e-3 ]
 

            * **dyn_params["state_tracking_algo"]** ( int ): the algorithm to keep track of the states' identities
 
                - 0: no state tracking
                - 1: Sato
                - 2: using the mincost, Munkres-Kuhn [ default: 2 ]


            * **dyn_params["MK_alpha"]** ( double ): Munkres-Kuhn alpha 
                (selects the range of orbitals included in reordering) [default: 0.0]


            * **dyn_params["MK_verbosity"]** ( double ): Munkres-Kuhn verbosity: 

                - 0: no extra output [ default ]
                - 1: print details



            * **dyn_params["entanglement_opt"]** ( int ): A selector of a method to couple the trajectories in this ensemble:

                - 0: no coupling across trajectories [ default ]
                - 1: ETHD
                - 2: ETHD3 (experimental)
                - 22: another flavor of ETHD3 (experimental)

            * **dyn_params["ETHD3_alpha"]** ( double ): Gaussian exponents that dresses up the trajectories in the ETHD3 method
                in the coordinate space, that is   ~exp(-alpha*(R-R0)^2 ) [ default: 0.0 ]


            * **dyn_params["ETHD3_beta"]** ( double ): Gaussian exponents that dresses up the trajectories in the ETHD3 method
                in the coordinate space, that is   ~exp(-alpha*(P-P0)^2 ) [ default: 0.0 ]


            * **dyn_params["decoherence_algo"]** ( int ): selector of the method to incorporate decoherence:

                - -1: no decoherence [ default ]
                - 0: MSDM
                - 1: ID-A


            * **dyn_params["decoh_rates"]** ( MATRIX(ntates, nstates) ): the matrix of dephasing rates [ units: a.u. of time ^-1 ]


            * **dyn_params["dt"]** ( double ): the nuclear and electronic integration
                timestep [ units: a.u. of time, default: 41.0 ]


            * **dyn_params["nsteps"]** ( int ): the number of NA-MD steps to do [ default: 1 ]


            * **dyn_params["prefix"]** ( string ): the name of the folder to be created by this function. 
                All the data files will be created in that folder


            * **dyn_params["output_level"]** ( int ): controls what info to return as the result of this function

                - -1: all return values are None [ default ]
                - 0: also, all return values are none 
                - 1: 1-D info - energies, SE and SH populations averaged over ensemble vs. time
                - 2: 2-D info - coordinates, momenta, SE amplitudes, and SH states populations 
                    for each individual trajectory vs. time 
                - 3: 3-D info - St, Hvib_adi, Hvib_dia matrices (energies, couplings, etc.) for each
                    individual trajectory vs. time 


            * **dyn_params["file_output_level"]** ( int ): controls what info to print out into files

                - -1: print nothing at all [ default ]
                - 0: 0-D info - just the parameters of the simulation
                - 1: 1-D info - energies, SE and SH populations averaged over ensemble vs. time 
                - 2: 2-D info - coordinates, momenta, SE amplitudes, and SH states populations 
                    for each individual trajectory vs. time 
                - 3: 3-D info - St, Hvib_adi, Hvib_dia matrices (energies, couplings, etc.) for each
                    individual trajectory vs. time 
            

        compute_model ( PyObject ): the pointer to the Python function that performs the Hamiltonian calculations

        _model_params ( dictionary ): contains the selection of a model and the parameters 
            for that model Hamiltonian

        rnd ( Random ): random numbers generator object



    Returns:
        tuple: with the elements of the tuple listed below.

            * obs_T ( list of `nsteps` doubles ): time [units: a.u.]
            * obs_q ( list of `nsteps` MATRIX(nnucl, ntraj) ): coordinates of all trajectories [ units: Bohr ]
            * obs_p ( list of `nsteps` MATRIX(nnucl, ntraj) ): momenta of all trajectories [ units: a.u. ]
            * obs_Ekin ( list of `nsteps` doubles ): average kinetic energy of an ensemble of trajectories [units: a.u.]
            * obs_Epot ( list of `nsteps` doubles ): average potential energy of an ensemble of trajectories [units: a.u.]
            * obs_Etot ( list of `nsteps` doubles ): average total energy of an ensemble of trajectories [units: a.u.]
            * obs_dEkin ( list of `nsteps` doubles ): standard deviation of kinetic energy of an ensemble of trajectories [units: a.u.]
            * obs_dEpot ( list of `nsteps` doubles ): standard deviation of potential energy of an ensemble of trajectories [units: a.u.]
            * obs_dEtot ( list of `nsteps` doubles ): standard deviation of total energy of an ensemble of trajectories [units: a.u.]
            * obs_Cadi ( list of `nsteps` CMATRIX(nadi, ntraj) ): amplitudes of adiabatic electronic states of all trajectories 
            * obs_Cdia ( list of `nsteps` CMATRIX(ndia, ntraj) ): amplitudes of diabatic electronic states of all trajectories 
            * obs_dm_adi ( list of `nsteps` CMATRIX(nadi, nadi) ): ensemble-averaged density matrix in adiabatic basis
            * obs_dm_dia ( list of `nsteps` CMATRIX(ndia, ndia) ): ensemble-averaged density matrix in diabatic basis
            * obs_pop ( list of `nsteps` MATRIX(nadi, 1) ): ensemble-averaged TSH populations of adiabatic states
            * obs_states ( list of `nsteps` of lists of `ntraj` ints):  # indices of the quantum states of each trajectory
            * obs_hvib_adi ( list of `ntraj` lists, each being a list of `nsteps` objects CMATRIX(nadi, nadi) ): trajectory-resolved
                vibronic Hamiltonians for each timestep in the adiabatic representation
            * obs_hvib_dia ( list of `ntraj` lists, each being a list of `nsteps` objects CMATRIX(ndia, ndia) ): trajectory-resolved
                vibronic Hamiltonians for each timestep in the diabatic representation
            * obs_St ( list of `ntraj` lists, each being a list of `nsteps` objects CMATRIX(nadi, nadi) ): trajectory-resolved
                time-overlaps of the adiabatic wavefunctions
            * obs_U ( list of `ntraj` lists, each being a list of `nsteps` objects CMATRIX(ndia, nadi) ): trajectory-resolved
                diabatic-to-adiabatic transformation matrix


        .. note::
          the elements are None, if they are excluded by the input optinos

    """

        
    # Create copies of the input dynamical variables, so we could run several such
    # functions with the same input variables without worries that they will be altered
    # inside of each other

    model_params = dict(_model_params)
    dyn_params = dict(_dyn_params)

    q = MATRIX(_q)
    p = MATRIX(_p)
    iM = MATRIX(_iM)
    Cdia = CMATRIX(_Cdia)
    Cadi = CMATRIX(_Cadi)
    states = intList()
    projectors = CMATRIXList()
    
    
    for i in range(len(_states)):
        states.append(_states[i])
        projectors.append(CMATRIX(_projectors[i]))

    DR = MATRIX(len(_states), len(_states))
    AG = MATRIX(len(_states), len(_states))


    # Parameters and dimensions
    critical_params = [  ] 
    default_params = { "rep_tdse":1, "rep_ham":0, "rep_sh":1, "rep_lz":0, "tsh_method":-1,
                       "force_method":1, "nac_update_method":1, "rep_force":1,
                       "use_boltz_factor":0, "Temperature":300.0, "do_reverse":1, "vel_rescale_opt":-1,
                       "do_phase_correction":1, "tol":1e-3,
                       "state_tracking_algo":2, "MK_alpha":0.0, "MK_verbosity":0, 
                       "entanglement_opt":0, "ETHD3_alpha":0.0, "ETHD3_beta":0.0, 
                       "decoherence_algo":-1, "decoherence_rates":DR,
                       "decoherence_times_type":0, "decoherence_C_param":1.0, 
                       "decoherence_eps_param":0.1, "dephasing_informed":0,
                       "ave_gaps":AG, "instantaneous_decoherence_variant":1, "collapse_option":0,
                       "ensemble":0, "thermostat_params":{},
                       "dt":1.0*units.fs2au, "nsteps":1, 
                       "output_level":-1, "file_output_level":-1, "prefix":"tmp"
                     }

    comn.check_input(dyn_params, default_params, critical_params)
        
    prefix = dyn_params["prefix"]    
    rep_tdse = dyn_params["rep_tdse"]
    nsteps = dyn_params["nsteps"]
    dt = dyn_params["dt"]
    tol = dyn_params["tol"]
    output_level = dyn_params["output_level"]
    file_output_level = dyn_params["file_output_level"]
    do_phase_correction = dyn_params["do_phase_correction"]
    state_tracking_algo = dyn_params["state_tracking_algo"]
    force_method = dyn_params["force_method"]
    
    ndia = Cdia.num_of_rows
    nadi = Cadi.num_of_rows
    nnucl= q.num_of_rows
    ntraj= q.num_of_cols



    obs_T = None
    obs_q, obs_p = None, None
    obs_Ekin, obs_Epot, obs_Etot = None, None, None
    obs_dEkin, obs_dEpot, obs_dEtot = None, None, None
    obs_Cadi, obs_Cdia = None, None 
    obs_dm_adi, obs_dm_dia = None, None
    obs_pop, obs_states = None, None
    obs_hvib_adi, obs_hvib_dia, obs_St, obs_U = None, None, None, None
         
    if output_level>=1:
        obs_T = [] # time
        obs_Ekin = []  # average kinetic energy 
        obs_Epot = []  # average potential energy 
        obs_Etot = []  # average total energy 
        obs_dEkin = []  # kinetic energy fluctuation
        obs_dEpot = []  # potential energy fluctuation
        obs_dEtot = []  # total energy fluctuation
        obs_dm_adi = []  # average SE-based density matrix in adiabatic basis
        obs_dm_dia = []  # average SE-based density matrix in diabatic basis
        obs_pop = []  # average SH-based populations adiabatic basis

    if output_level>=2:
        obs_q = [] # coordinates of all trajectories
        obs_p = [] # momenta of all trajectories
        obs_Cadi = []  # average TD-SE amplitudes in the adiabatic basis
        obs_Cdia = []  # average TD-SE amplitudes in the diabatic basis
        obs_states = []  # indices of the quantum states of each trajectory

    if output_level>=3:
        obs_hvib_adi = []  # vibronic Hamiltonians for each trajectory in adiabatic rep.
        obs_hvib_dia = []  # vibronic Hamiltonians for each trajectory in diabatic rep.
        obs_St = [] # time-overlaps for each trajectory in the adiabatic rep.
        obs_U = []  # dia-to-adia transformation matrices for each trajectory

        for tr in range(ntraj):
            obs_hvib_adi.append([])
            obs_hvib_dia.append([])
            obs_St.append([])
            obs_U.append([])


    # Create an output directory, if not present    
    if file_output_level>=0:
        if not os.path.isdir(prefix):
            os.mkdir(prefix)
            
        # Simulation parameters                    
        f = open("%s/_dyn_params.txt" % (prefix),"w")
        f.write( str(dyn_params) );  f.close()
    
        f = open("%s/_model_params.txt" % (prefix),"w")
        f.write( str(model_params) );  f.close()    

    # Ensemble-resolved info (so 1D)
    if file_output_level >= 1:
        f = open("%s/energies.txt" % (prefix), "w"); f.close()   # average kinetic, potential, and total energies and their fluctuations                                                                
        f = open("%s/D_adi.txt" % (prefix), "w"); f.close()      # average SE-based density matrix in adiabatic basis
        f = open("%s/D_dia.txt" % (prefix), "w"); f.close()      # average SE-based density matrix in diabatic basis
        f = open("%s/SH_pop.txt" % (prefix), "w"); f.close()     # average SH-based populations adiabatic basis
    
    # Trajectory-resolved 1D info (so 2D)
    if file_output_level >= 2:
        f = open("%s/q.txt" % (prefix), "w"); f.close()          # coordinates of all trajectories
        f = open("%s/p.txt" % (prefix), "w"); f.close()          # momenta of all trajectories    
        f = open("%s/C_adi.txt" % (prefix), "w"); f.close()      # TD-SE amplitudes in the adiabatic basis
        f = open("%s/C_dia.txt" % (prefix), "w"); f.close()      # TD-SE amplitudes in the diabatic basis    
        f = open("%s/states.txt" % (prefix), "w"); f.close()     # indices of the quantum states of each trajectory

    if file_output_level >= 3:    
        # Trajectory-resolved 2D info (so 3D)
        for tr in range(ntraj):
            f = open("%s/Hvib_adi_%i.txt" % (prefix, tr), "w"); f.close()   # vibronic Hamiltonians for each trajectory in adiabatic rep. 
            f = open("%s/Hvib_dia_%i.txt" % (prefix, tr), "w"); f.close()   # vibronic Hamiltonians for each trajectory in diabatic rep. 
            f = open("%s/St_%i.txt" % (prefix, tr), "w"); f.close()   # time-overlaps along each trajectory 
            f = open("%s/basis_transform_%i.txt" % (prefix, tr), "w"); f.close() # dia-to-adi transformationfor each trajectory. 
        

    # ======= Hierarchy of Hamiltonians =======
    ham = nHamiltonian(ndia, nadi, nnucl)
    ham.add_new_children(ndia, nadi, nnucl, ntraj)
    ham.init_all(2,1)
    model_params.update({"timestep":0})
    
    update_Hamiltonian_q(dyn_params, q, projectors, ham, compute_model, model_params)
    update_Hamiltonian_p(dyn_params, ham, p, iM)  


    U = []
    for tr in range(ntraj):
        U.append(ham.get_basis_transform(Py2Cpp_int([0, tr]) ))


                
    # Do the propagation
    for i in range(nsteps):
    
        #============ Compute and output properties ===========        
        # Amplitudes, Density matrix, and Populations
        if rep_tdse==0:
            # Diabatic to raw adiabatic
            ham.ampl_dia2adi(Cdia, Cadi, 0, 1) 
            # Raw adiabatic to dynamically-consistent adiabatic
            Cadi = dynconsyst_to_raw(Cadi, projectors)

        elif rep_tdse==1:
            ham.ampl_adi2dia(Cdia, Cadi, 0, 1)

        dm_dia, dm_adi = tsh_stat.compute_dm(ham, Cdia, Cadi, projectors, rep_tdse, 1)        
        pops = tsh_stat.compute_sh_statistics(nadi, states)


        # Energies 
        if force_method in [0, 1]:
            Ekin, Epot, Etot, dEkin, dEpot, dEtot = tsh_stat.compute_etot_tsh(ham, p, Cdia, Cadi, projectors, states, iM, rep_tdse)
        elif force_method in [2]:
            Ekin, Epot, Etot, dEkin, dEpot, dEtot = tsh_stat.compute_etot(ham, p, Cdia, Cadi, projectors, iM, rep_tdse)


        # ============== Debug ====================
        nrm = (Cadi.H() * Cadi ).tr()             

        x = 0.0
        for traj in range(ntraj):
            x += (projectors[traj].H() * projectors[traj]).tr()

        #print(F"step = {i}  nrm = {nrm} projectors measure = {x}")
        #============ End of Debug =================
            
        # Memory output
        if output_level >= 1:
            obs_T.append(i*dt) 
            obs_Ekin.append(Ekin)
            obs_Epot.append(Epot)
            obs_Etot.append(Etot)
            obs_dEkin.append(dEkin)
            obs_dEpot.append(dEpot)
            obs_dEtot.append(dEtot)
            obs_dm_adi.append(CMATRIX(dm_adi))
            obs_dm_dia.append(CMATRIX(dm_dia))
            obs_pop.append(MATRIX(pops))

        # Memory output
        if output_level >= 2:
            obs_q.append(MATRIX(q))
            obs_p.append(MATRIX(p))        
            obs_Cadi.append(CMATRIX(Cadi))
            obs_Cdia.append(CMATRIX(Cdia))
            obs_states.append(list(states))

        # File output
        res12 = q, p, Ekin, Epot, Etot, dEkin, dEpot, dEtot, Cadi, Cdia, dm_adi, dm_dia, pops, states
        dynamics_io.print_results12(i, dt, res12, prefix, file_output_level)

    
        for tr in range(ntraj):
            x = ham.get_basis_transform(Py2Cpp_int([0, tr]) )            
            St = U[tr] * x.H()             
            U[tr] = CMATRIX(x)

            # Memory output            
            if output_level >= 3:
                obs_hvib_adi[tr].append( CMATRIX(ham.get_hvib_adi(Py2Cpp_int([0, tr])) ) )
                obs_hvib_dia[tr].append( CMATRIX(ham.get_hvib_dia(Py2Cpp_int([0, tr])) ) )
                obs_St[tr].append( CMATRIX(St) )
                obs_U[tr].append( CMATRIX(x) )

            # File output   
            res3 = ham.get_hvib_adi(Py2Cpp_int([0, tr])), ham.get_hvib_dia(Py2Cpp_int([0, tr])), St, U[tr]
            dynamics_io.print_results3(i, dt, res3, prefix, file_output_level, tr)


        #============ Propagate ===========        
        model_params.update({"timestep":i})        
        if rep_tdse==0:
            compute_dynamics(q, p, iM, Cdia, projectors, states, ham, compute_model, model_params, dyn_params, rnd)
        elif rep_tdse==1:
            compute_dynamics(q, p, iM, Cadi, projectors, states, ham, compute_model, model_params, dyn_params, rnd)

 

    return obs_T, obs_q, obs_p, obs_Ekin, obs_Epot, obs_Etot, obs_dEkin, obs_dEpot, obs_dEtot, \
           obs_Cadi, obs_Cdia, obs_dm_adi, obs_dm_dia, obs_pop, obs_states, obs_hvib_adi, obs_hvib_dia, obs_St, obs_U




def generic_recipe(q, p, iM, _dyn_params, compute_model, _model_params, _init_elec, rnd):
    """
    This function initializes electronic variables and Hamiltonians for a given set of 
    nuclear variables (those are needed, not initialized here) and then runs the calculations 
    on multiple trajectories

    Args:
        q ( MATRIX(ndof, ntraj) ): the coordinates for all DOFs for all trajectories [ units: a.u.]
        p ( MATRIX(ndof, ntraj) ): the momenta for all DOFs for all trajectories [ units: a.u. ]
        iM ( MATRIX(ndof, 1) ): masses of all nuclear DOFs - same for all trajectories [ units: a.u. ]

        _dyn_params ( dictionary ): control parameters for running the dynamics
            see the documentation of the :func:`libra_py.dynamics.run_dynamics` function

        compute_model ( PyObject ): the pointer to the Python function that performs the Hamiltonian calculations

        _model_params ( dictionary ): contains the selection of a model and the parameters 
            for that model Hamiltonian. In addition to all the parameters, should contain the 
            key

            * **_model_params["model0"]** ( int ): the selection of the model to be used to compute diabatic-to-adiabatic
                transformation. The particular value chosen here depends on the way the `compute_model` functions is defined

               :Note: the function selected should be able to generate the transformation matrix.

        _init_elec ( dictionary ): control parameters for initialization of electronic variables
            see the documentation of the :func:`libra_py.dynamics.init_electronic_dyn_var` function

        rnd ( Random ): random numbers generator object

    Returns:
        tuple : tuple returned by the :func:`libra_py.dynamics.run_dynamics` function



    """


    critical_params = [ "model0" ] 
    default_params = {  }
    comn.check_input(_model_params, default_params, critical_params)


    # Internal parameters
    nnucl, ntraj = q.num_of_rows, q.num_of_cols

    # Initialize electronic variables - either diabatic or adiabatic
    Cdia, Cadi, projectors, states = init_electronic_dyn_var(_init_elec, rnd)
    ndia, nadi = Cdia.num_of_rows, Cadi.num_of_rows


    # In case the initial conditions rep and the propagation rep are different,
    # compute the diabatic-to-adiabatic transformation matrices and
    # transform the amplitudes accordingly
    ham = nHamiltonian(ndia, nadi, nnucl)
    ham.add_new_children(ndia, nadi, nnucl, ntraj)
    ham.init_all(2,1)     

    if _init_elec["rep"]==0:
        if _dyn_params["rep_tdse"]==1:
            model_params = dict(_model_params)
            model_params.update({"model":_model_params["model0"]})
            update_Hamiltonian_q({"rep_tdse":1, "rep_ham":0}, q, projectors, ham, compute_model, model_params )

            Cadi = transform_amplitudes(0, 1, Cdia, ham) 

    elif _init_elec["rep"]==1:    
        if _dyn_params["rep_tdse"]==0:
            model_params = dict(_model_params)
            model_params.update({"model":_model_params["model0"]})
            update_Hamiltonian_q({"rep_tdse":1, "rep_ham":0}, q, projectors, ham, compute_model, model_params )

            Cdia = transform_amplitudes(1, 0, Cdia, ham) 

    res = run_dynamics(q, p, iM, Cdia, Cadi, projectors, states, _dyn_params, compute_model, _model_params, rnd)


    return res



def run_multiple_sets(init_cond, _dyn_params, compute_model, _model_params, _init_nucl, _init_elec, rnd):
    """
    This function runs a series of NA-MD calculations for (potentially) different initial 
    conditions and returns an object that contains results for each of the 
    runs, and/or just stores the results into the corresponding folders.
    
    Args:

      init_cond ( list of 3-element lists ): mean coordinates and momenta as well as the masses of all DOFs

          The format is this:
          init_cond = 
            init_cond[0] = [q0, p0, M0 ]  <- 0-th set of initial variables
            init_cond[1] = [q1, p1, M1 ]  <- 1-st set of initial variables
            ...                           etc.
            init_cond[N] = [qN, pN, MN ]
          
          Here, q0 = [q00, q01, ... q0<ndof>], same for p0 and M0

      _dyn_params ( list of dictionaries ): control parameters for running the dynamics.
          Each item of the list is a dictionary for a particular set of initial conditions.
          see the documentation of the :func:`libra_py.dynamics.run_dynamics` function

      compute_model ( list of PyObject objects ): the pointers to the Python functions 
          that performs the Hamiltonian calculations, one per each initial conditions set

      _model_params ( list of dictionaries ): contains the selection of a model and the parameters 
          for that model Hamiltonian. Each item corresponds to a particular simulation set.
          In addition to all the parameters, should contain the key:

          * **_model_params["model0"]** ( int ): the selection of the model to be used to compute diabatic-to-adiabatic
          transformation. The particular value chosen here depends on the way the `compute_model` functions is defined

              :Note: the function selected should be able to generate the transformation matrix.

      _init_nucl ( list of dictionary ): control parameters for initialization of nuclear variables
          One item per each initial conditions set.
          see the documentation of the :func:`libra_py.dynamics.init_nuclear_dyn_var` function

      _init_elec ( dictionary ): control parameters for initialization of electronic variables.
          One item per each initial conditions set.
          see the documentation of the :func:`libra_py.dynamics.init_electronic_dyn_var` function

      rnd ( Random ): random numbers generator object
    
    Returns:
        list : each element of the list is a tuple returned by the
            :func:`libra_py.dynamics.run_dynamics` function
            Each such element represents the results for a particular set of calculations

    """

    
            
    res = []
    for icond_indx, icond in enumerate(init_cond):

        output_prefix = _dyn_params[icond_indx]["prefix"]
        dyn_params = dict(_dyn_params[icond_indx])
        
        q0, p0, M0 = icond[0], icond[1], icond[2]  # each one is a ndof-elements list
        q, p, iM = init_nuclear_dyn_var(q0, p0, M0, _init_nucl[icond_indx], rnd)

        dyn_params.update({"prefix":F"{output_prefix}_{icond_indx}"})

        res_i = generic_recipe(q, p, iM, dyn_params, compute_model[icond_indx], _model_params[icond_indx], _init_elec[icond_indx], rnd)        
        res.append(res_i)
    
    return res



