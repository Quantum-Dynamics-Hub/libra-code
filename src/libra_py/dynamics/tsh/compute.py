#*********************************************************************************                     
#* Copyright (C) 2019-2022 Alexey V. Akimov                                                   
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
           * init_electronic_dyn_var(params, isNBRA, rnd)
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
import time
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

import util.libutil as comn
import libra_py.units as units
import libra_py.data_outs as data_outs
import libra_py.data_savers as data_savers
import libra_py.tsh_stat as tsh_stat
#import libra_py.dynamics as dynamics_io

from . import save


#def run_dynamics(_q, _p, _iM, _Cdia, _Cadi, _projectors, _states, _dyn_params, compute_model, _model_params, rnd):
def run_dynamics(dyn_var, _dyn_params, ham, compute_model, _model_params, rnd):
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

            ///===============================================================================
            ///================= Computing Hamiltonian-related properties ====================
            ///===============================================================================
      
            * **dyn_params["rep_tdse"]** ( int ): selects the representation in which 
                nuclear/electronic (Ehrenfest core) dynamics is executed. The representation 
                used to integrate TD-SE

                - 0: diabatic representation
                - 1: adiabatic representation [ default ]


            * **dyn_params["rep_ham"]** (int): The representation of the Hamiltonian update: 
                This is the representation in which the computed properties are assumed to be
                For instance, we may have set it to 1, to read the adiabatic energies and couplings,
                to bypass the diabatic-to-adiabatic transformation, which may be useful in some atomistic
                calculations, or with the NBRA

                - 0: diabatic   [ default ]
                - 1: adiabatic                              


            * **dyn_params["rep_sh"]** ( int ): selects the representation which is 
                used to perform surface hopping

                - 0: diabatic representation
                - 1: adiabatic representation [ default ]


            * **dyn_params["rep_lz"]** ( int ): The representation to compute LZ probabilitieis:
 
                - 0: diabatic [ default ] 
                - 1: adiabatic 


            * **dyn_params["rep_force"]** ( int ): In which representation to compute forces.
                To clarify - the forces in both representations shall be equivalent, so this flag actually
                selects the type of the properties needed to compute such forces.
                For instance, if it is set to 0, we may be using the derivatives
                of the diabatic Hamiltonians, diabatic states overlaps, etc.

                - 0: diabatic
                - 1: adiabatic [ default ]


            * **dyn_params["force_method"]** ( int ): How to compute forces in the dynamics: 

                - 0: don't compute forces at all - e.g. we do not really need them
                - 1: state-specific  as in the TSH or adiabatic (including adiabatic excited states) [ default ]
                - 2: Ehrenfest


            * **dyn_params["time_overlap_method"]** ( int ): the flag to select how the time-overlaps are computed:

                - 0: on-the-fly, based on the wavefunction ("basis_transform") [ default ]
                - 1: read externally ( via "time_overlap_adi" ), use this option with the NBRA, but don't forget
                    to set up its update in the Hamiltonian update functions (aka "compute_model")


            * **dyn_params["nac_update_method"]** ( int ): How to update NACs and vibronic Hamiltonian before
               electronic TD-SE propagation:

                - 0: don't update them (e.g. for simplest NAC)
                - 1: update according to changed momentum and existing derivative couplings [ default ]


            * **dyn_params["do_phase_correction"]** ( int ): the algorithm to correct phases on adiabatic states
 
                - 0: no phase correction
                - 1: according to our phase correction algorithm [ default ]


            * **dyn_params["phase_correction_tol"]** ( double ): The minimal magnutude of the matrix element 
               for which we'll be computing the phase correction. If the overlap is zero, then we don't really
               care about the phase, but if it is not, then this parameter sets out threshold for when we do. [ default: 1e-3 ]
 

            * **dyn_params["state_tracking_algo"]** ( int ): the algorithm to keep track of the states' identities
 
                - 0: no state tracking
                - 1: Sato
                - 2: using the mincost, Munkres-Kuhn [ default ]
                - 3: experimental stochastic algorithm, the original version with elimination (known problems)
                - 32: experimental stochastic algorithms with all permutations (too expensive)
                - 33: the improved stochastic algorithm with good scaling and performance, on par with the mincost


            * **dyn_params["MK_alpha"]** ( double ): Munkres-Kuhn alpha 
                (selects the range of orbitals included in reordering) [default: 0.0]


            * **dyn_params["MK_verbosity"]** ( double ): Munkres-Kuhn verbosity: 

                - 0: no extra output [ default ]
                - 1: print details for debugging

            * **dyn_params["convergence"]** ( int ):  A swtich for stochastic reordering algorithm 3 (and variants) 
                to choose what happens when an acceptable permutation isn't generated in the set number of attempts:

                - 0: returns the identity permutation (does not require convergence) [ default ]
                - 1: exits and prints an error (requires convergence)

            * **dyn_params["max_number_attempts"]** ( int ): The maximum number of hops that an be 
                attempted before either choosing the identity or exiting in stochastic reordering algorithm 3 (and variants). 
                Default: 100

            * **dyn_params["min_probability_reordering"]** ( double ): The probability threshold for stochastic 
                state reordering algorithm. If a probability for a multi-state stransition is below this value,
                it will be disregarded and set to 0
                The rest of the probabilities will be renormalized
                Default: 0.0 


            ///===============================================================================
            ///================= Surface hopping: proposal, acceptance =======================
            ///===============================================================================


            * **dyn_params["tsh_method"]** ( int ): Formula for computing SH probabilities:
   
                - -1: adiabatic - no hops [ default ]
                - 0: FSSH
                - 1: GFSH 
                - 2: MSSH
                - 3: DISH

            * **dyn_params["hop_acceptance_algo"]** ( int ): Options to control the acceptance of the proposed hops:

                - 0: accept all proposed hops  [ default ]

                - 10: based on adiabatic energy - accept only those hops that can obey the energy conservation with 
                    adiabatic potential energies
                - 11: based on diabatic energy - same as 10, but we use diabatic potential energies

                - 20: based on derivative coupling vectors - accept only those hops that can obey the energy conservation
                    by rescaling nuclear velocities along the directions of derivative couplings for the quantum nuclear DOFs                   
                - 21: based on difference of state-specific forces - same as 20, but the rescaling is done along the vector
                    parallel to the difference of adiabatic forces on initial and target states

                - 31: accept hops with the probability taken from the quantum Boltzmann distribution
                - 32: accept hops with the probability taken from the classical Maxwell-Boltzmann distribution
                - 33: accept hops with the probability taken from the updated quantum Boltzmann distribution (experimental)


            * **dyn_params["momenta_rescaling_algo"]** ( int ): Options to control nuclear momenta changes upon successful 
                or frustrated hops.

                - 0: don't rescale [ default ]

                - 100: based on adiabatic energy, don't reverse on frustrated hops
                - 101: based on adiabatic energy, reverse on frustrated hops
                - 110: based on diabatic energy, don't reverse on frustrated hops
                - 111: based on diabatic energy, reverse on frustrated hops

                - 200: along derivative coupling vectors, don't reverse on frustrated hops
                - 201: along derivative coupling vectors, reverse on frustrated hops
                - 210: along difference of state-specific forces, don't reverse on frustrated hops
                - 211: along difference of state-specific forces, reverse on frustrated hops


            * **dyn_params["use_boltz_factor"]** ( int ): A flag to scale the proposed hopping probabilities by the
                Boltzmann factor. This is needed for the libra_py/workflows/nbra calculations, where the hop proposal 
                probability also includes the factor to account for the hop acceptance probabilities

                - 0: don't scale [ default ]
                - 1: do scale


            ///===============================================================================
            ///================= Decoherence options =========================================
            ///===============================================================================


            * **dyn_params["decoherence_algo"]** ( int ): selector of the method to incorporate decoherence:

                - [-1]: no decoherence [ default ]
                - 0: SDM and alike
                - 1: instantaneous decoherence options (ID-S, ID-A, ID-C)


            * **dyn_params["sdm_norm_tolerance"]** ( double ): Corresponds to the "tol" parameter in the sdm function. It controls 
                how much the norm of the old state can be larger than 1.0  before the code stops with the error message [ default : 0.0 ]

                Note: only matters if decoherence_algo == 0


            * **dyn_params["dish_decoherence_event_option"]** ( int ):  Selects the how to sample decoherence events in the DISH:

                - 0: compare the coherence time counter with the decoherence time (simplified DISH) 
                - 1: compare the coherence time counter with the time drawn from the exponential distribution
                    with the parameter lambda = 1/decoherence time - this distribution corresponds to 
                    the statistics of wait times between the Poisson-distributed events (decoherence)
                    This is what the original DISH meant to do  [ default ]

                Note: only matters if dyn_params["tsh_method"] == 3


            * **dyn_params["decoherence_times_type"]** ( int ): Type of dephasing times/rates calculation:

                - 0: use the rates read out from the input, need `decoherence_rates` input  [default]
                - 1: use the energy-based decoherence method (EDC)    


            * **dyn_params["decoherence_C_param"]** ( double ): An empirical parameter used in the EDC method: [ default: 1.0 Ha]


            * **dyn_params["decoherence_eps_param"]** ( double ): An empirical parameter used in the EDC method: [ default: 0.1 Ha]

            * **dyn_params["dephasing_informed"]** ( int ): A flag to apply the dephasing-informed approach of Sifain et al. 
                to correct dephasing times: 

                - 0: don't apply [ default ]
                - 1: use it, will need also the `ave_gaps` input


            * **dyn_params["instantaneous_decoherence_variant"]** ( int ): Option to control the instantaneous decoherence methodology,
                only used with `decoherence_algo == 1`

                - 0: ID-S
                - 1: ID-A [ default ]
                - 2: ID-C - consistent ID (experimental)


            * **dyn_params["collapse_option"]** ( int ): How to collapse wavefunction amplitudes in the decoherence schemes:

                - 0: by rescaling the magnitude of the amplitude vector elements, but preserving "phase" [ default ]
                - 1: by resetting the amplitudes to 1.0+0.0j. This option changes phase 


            * **dyn_params["decoherence_rates"]** ( MATRIX(ntates, nstates) ): the matrix of dephasing rates [ units: a.u. of time ^-1 ]


            * **dyn_params["ave_gaps"]** ( MATRIX(ntates, nstates) ): A matrix that contains the averaged moduli of the energy gaps:
                E_ij = <|E_i - E_j|>.  It is needed when dephasing_informed option is used




            ///===============================================================================
            ///================= Entanglement of trajectories ================================
            ///===============================================================================


            * **dyn_params["entanglement_opt"]** ( int ): A selector of a method to couple the trajectories in this ensemble:

                - 0: no coupling across trajectories [ default ]
                - 1: ETHD
                - 2: ETHD3 (experimental)
                - 22: another flavor of ETHD3 (experimental)


            * **dyn_params["ETHD3_alpha"]** ( double ): Gaussian exponents that dresses up the trajectories in the ETHD3 method
                in the coordinate space, that is   ~exp(-alpha*(R-R0)^2 ) [ default: 0.0 ]


            * **dyn_params["ETHD3_beta"]** ( double ): Gaussian exponents that dresses up the trajectories in the ETHD3 method
                in the coordinate space, that is   ~exp(-alpha*(P-P0)^2 ) [ default: 0.0 ]


            ///===============================================================================
            ///================= Bath, Constraints, and Dynamical controls ===================
            ///===============================================================================


            * **dyn_params["Temperature"]** ( double ): Temperature of the system. This parameter
                could be used even in the NVE simulations e.g. as a parameters to compute hop 
                acceptance probabilities based on Boltzmann factors [ default: 300 K]


            * **dyn_params["ensemble"]** ( int ): Which ensemble to use in the dynamics. 

                - 0: NVE [ default ]
                - 1: NVT


            * **dyn_params["thermostat_params"]** ( dict ): Parameters controlling the thermostat,
                only relevant for `ensemble = 1`  [ default: {} ]


            * **dyn_params["thermostat_dofs"]** ( list of ints ): Thermostat DOFs
                This list contains the indices of nuclear DOFs which shall be coupled to a thermostat directly.
                [ default: [] ]


            * **dyn_params["quantum_dofs"]** ( list of ints ): Quantum-classical partitioning
                This list of integers contains the indices of nuclear DOFs which chall be treated "quantum-mechanically", well
                including with TSH that is. These DOFs will determine the velocity-rescaling-based acceptance of the hops,
                and these DOFs will be rescaled when the transition is accepted 
                [ default: all DOFs ]

                Note: this default is different from the C++ default, where we just use [0]


            * **dyn_params["constrained_dofs"]** ( list of ints ):  Constrained DOFs
                This list of integers contains the indices of the nuclear DOFs to be constrained - their momenta will be constantly 
                reset to zero, so the corresponding coordinates will stay fixed
                [ default: [] ]


            * **dyn_params["dt"]** ( double ): the nuclear and electronic integration
                timesteps [ units: a.u. of time, default: 41.0 a.u. = 1 fs ]


            ///===============================================================================
            ///================= Variables specific to Python version: saving ================
            ///===============================================================================


            * **dyn_params["nsteps"]** ( int ): the number of the dynamical steps to do [ default: 1 ]


            * **dyn_params["prefix"]** ( string ): the name of the folder to be created by this function. 
                All the data files will be created in that folder [ default: "out" ]


            * **dyn_params["hdf5_output_level"]** ( int ): controls what info to save into HDF5 files (as we run)

                - -1: all returned values are None [ default ]
                - 0: also, all return values are none  
                - 1: 1-D info - Relevant variables to request to save: ["timestep", "time", "Ekin_ave", "Epot_ave", "Etot_ave", 
                    "Etherm", "E_NHC", "dEkin_ave", "dEpot_ave", "dEtot_ave"]
                - 2: 2-D info - Relevant variables to request to save: ["states"]
                - 3: 3-D info - Relevant variables to request to save: [ "SH_pop", "SH_pop_raw",
                    "D_adi", "D_adi_raw", "D_dia", "D_dia_raw", "q", "p", "Cadi", "Cdia" ]
                - 4: 4-D info - Relevant variables to request to save: ["hvib_adi", "hvib_dia", "St", "basis_transform", "projector"]


            * **dyn_params["mem_output_level"]** ( int ): controls what info to save into HDF5 files (all at the end)

                Same meaning and output as with hdf5_output_level, except all the variables are first stored in memory
                (while the calcs are running) and then they are written into the HDF5 file at the end of the calculations.
                This is a much faster version of hdf5 saver.


            * **dyn_params["txt_output_level"]** ( int ): controls what info to save into TXT files 

                Same meaning and output as with hdf5_output_level, except all the variables are written as text files. [ default: -1 ] 


            * **dyn_params["use_compression"]** ( int ): whether to use the compression when writing the info into HDF5 files
                 
                - 0: don't use it [ default ]
                - 1: use the compression, with the level defined by the `compression_level` parameter. However,
                     we found that the compression slows down the processing significantly and sometimes
                     creates even larger files, so it is not recommended to use it

            * **dyn_params["compression_level"]** ( list of 3 ints ): the compression level for ints, doubles, and complex, respectively.
                The larger the value, the higher the compression. [ default: [0, 0, 0] ]


            * **dyn_params["progress_frequency"]** ( double ):  at what intervals print out some "progress" messages. For
                instance, if you have `nsteps = 100` and `progress_frequency = 0.1`, the code will notify you every 10 steps. [ default : 0.1 ]

                Note: this doesn't affect the frequency of the data saving to the file - at this point, we save the information for every
                timestep

                                                                                                                          
            * **dyn_params["properties_to_save"]** ( list of string ): describes what properties to save to the HDF5 files. Note that
                if some properties are not listed in this variable, then they are not saved, even if `mem_output_level` suggests they may be
                saved. You need to BOTH set the appropriate `mem_output_level` AND `properties_to_save`

                [ default:  [ "timestep", "time", "Ekin_ave", "Epot_ave", "Etot_ave", 
                           "dEkin_ave", "dEpot_ave", "dEtot_ave", "states", "SH_pop", "SH_pop_raw",
                           "D_adi", "D_adi_raw", "D_dia", "D_dia_raw", "q", "p", "Cadi", "Cdia", 
                           "hvib_adi", "hvib_dia", "St", "basis_transform", "projector"
                ] 


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

    nstates = model_params["nstates"]

    # Parameters and dimensions
    critical_params = [  ] 
    default_params = {}
    #================= Computing Hamiltonian-related properties ====================
    default_params.update( { "rep_tdse":1, "rep_ham":0, "ham_update_method":1, "ham_transform_method":1,
                             "rep_sh":1, "rep_lz":0, "rep_force":1,
                             "force_method":1, "enforce_state_following":0, "enforced_state_index":0, 
                             "time_overlap_method":0, "nac_update_method":1, "nac_algo":0,
                             "hvib_update_method":1, "do_ssy":0, 
                             "do_phase_correction":1, "phase_correction_tol":1e-3,
                             "state_tracking_algo":2, "MK_alpha":0.0, "MK_verbosity":0,
                             "convergence":0,  "max_number_attempts":100, "min_probability_reordering":0.0, 
                             "is_nbra":0, "icond":0, "nfiles":-1
                           } )

    #================= Surface hopping: proposal, acceptance =======================
    default_params.update( { "tsh_method":-1, "hop_acceptance_algo":0, "momenta_rescaling_algo":0,
                             "use_boltz_factor":0
                           } )

    #================= Decoherence options =========================================
    default_params.update( { "decoherence_algo":-1, "sdm_norm_tolerance":0.0,
                             "dish_decoherence_event_option":1, "decoherence_times_type":-1, 
                             "decoherence_C_param":1.0, "decoherence_eps_param":0.1, 
                             "dephasing_informed":0, "instantaneous_decoherence_variant":1, 
                             "collapse_option":0,  
                             "decoherence_rates":MATRIX(nstates, nstates),
                             "ave_gaps":MATRIX(nstates,nstates),
                             "schwartz_decoherence_inv_alpha": MATRIX(nstates, 1)
                           } )

    #================= Entanglement of trajectories ================================
    default_params.update( { "entanglement_opt":0, "ETHD3_alpha":0.0, "ETHD3_beta":0.0   } )

    #================= Bath, Constraints, and Dynamical controls ===================
    default_params.update( { "Temperature":300.0, "ensemble":0, "thermostat_params":{},
                             "quantum_dofs":None, "thermostat_dofs":[], "constrained_dofs":[],
                             "dt":1.0*units.fs2au, "num_electronic_substeps":1,
                             "electronic_integrator":0
                           } )

    #================= Variables specific to Python version: saving ================
    default_params.update( { "nsteps":1, "prefix":"out", "prefix2":"out2",
                             "hdf5_output_level":-1, "mem_output_level":-1, "txt_output_level":-1, "txt2_output_level":-1,
                             "use_compression":0, "compression_level":[0,0,0], 
                             "progress_frequency":0.1,
                             "properties_to_save":[ "timestep", "time", "Ekin_ave", "Epot_ave", "Etot_ave", 
                                   "dEkin_ave", "dEpot_ave", "dEtot_ave", "states", "SH_pop", "SH_pop_raw",
                                   "se_pop_adi", "se_pop_dia", "sh_pop_adi",
                                   "D_adi", "D_adi_raw", "D_dia", "D_dia_raw", "q", "p", "f", "Cadi", "Cdia", 
                                   "hvib_adi", "hvib_dia", "St", "basis_transform", "projector" ] 
                           } )

    comn.check_input(dyn_params, default_params, critical_params)
               
    prefix = dyn_params["prefix"] 
    prefix2 = dyn_params["prefix2"] 
    rep_tdse = dyn_params["rep_tdse"]
    rep_ham = dyn_params["rep_ham"]
    nsteps = dyn_params["nsteps"]
    dt = dyn_params["dt"]
    phase_correction_tol = dyn_params["phase_correction_tol"]
    hdf5_output_level = dyn_params["hdf5_output_level"]
    mem_output_level = dyn_params["mem_output_level"]
    txt_output_level = dyn_params["txt_output_level"]
    txt2_output_level = dyn_params["txt2_output_level"]
    do_phase_correction = dyn_params["do_phase_correction"]
    state_tracking_algo = dyn_params["state_tracking_algo"]
    force_method = dyn_params["force_method"]
    properties_to_save = dyn_params["properties_to_save"]
    use_compression = dyn_params["use_compression"]
    compression_level = dyn_params["compression_level"]
    ensemble = dyn_params["ensemble"]
    time_overlap_method = dyn_params["time_overlap_method"]
    decoherence_algo = dyn_params["decoherence_algo"]
    is_nbra = dyn_params["is_nbra"]    
    icond = dyn_params["icond"]    # Setting the initial geoemtry in the dynamics    
    nfiles = dyn_params["nfiles"]  # The number of loaded Ham files
    tsh_method = dyn_params["tsh_method"] 

    #q = MATRIX(_q)
    #p = MATRIX(_p)
    #iM = MATRIX(_iM)
    #Cdia = CMATRIX(_Cdia)
    #Cadi = CMATRIX(_Cadi)
    states = intList()
    #projectors = CMATRIXList() 

    ndia = dyn_var.ndia #Cdia.num_of_rows
    nadi = dyn_var.nadi #Cadi.num_of_rows
    nnucl = dyn_var.ndof #q.num_of_rows
    ntraj = dyn_var.ntraj #q.num_of_cols


    #sys.exit(0)

    if(dyn_params["quantum_dofs"]==None):
        dyn_params["quantum_dofs"] = list(range(nnucl))

      
    if is_nbra == 1:
        for i in range(nstates):
            pass #states.append(_states[i])
        #projectors.append(CMATRIX(_projectors[0]))
    else:
        for i in range(nstates):
            pass #states.append(_states[i])
            #projectors.append(CMATRIX(_projectors[i]))


    # Initialize savers
    _savers = save.init_tsh_savers(dyn_params, model_params, nsteps, ntraj, nnucl, nadi, ndia)

    # Open and close the output files for further writing
    if _savers["txt_saver"]!=None:
        _savers["txt_saver"].save_data_txt( F"{prefix}", properties_to_save, "w", 0)

    if _savers["txt2_saver"]!=None:
        _savers["txt2_saver"].save_data_txt( F"{prefix2}", properties_to_save, "w", 0)


    #sys.exit(0)
    model_params.update({"timestep":icond})    
    #update_Hamiltonian_q(dyn_params, dyn_var, ham, compute_model, model_params)
    #update_Hamiltonian_p(dyn_params, dyn_var, ham)  
    update_Hamiltonian_variables( dyn_params, dyn_var, ham, ham, compute_model, model_params, 0)
    #sys.exit(0)
    update_Hamiltonian_variables( dyn_params, dyn_var, ham, ham, compute_model, model_params, 1)



    #sys.exit(0)

    #U = []
    #if is_nbra == 1:
    #    U.append(ham.get_basis_transform(Py2Cpp_int([0, 0]) ))
    #else:
    #    for tr in range(ntraj):
    #        U.append(ham.get_basis_transform(Py2Cpp_int([0, tr]) ))

    therm = ThermostatList();
    if ensemble==1:
        for traj in range(ntraj):
            therm.append( Thermostat( dyn_params["thermostat_params"] ) )
            therm[traj].set_Nf_t( len(dyn_params["thermostat_dofs"]) )
            therm[traj].init_nhc()

    #sys.exit(0)
    if decoherence_algo==2:
        dyn_var.allocate_afssh()
    elif decoherence_algo==3: # BCSH
        dyn_var.allocate_bcsh()
    if tsh_method==5: # DISH
        dyn_var.allocate_dish()
    if tsh_method==7: #  FSSH2
        dyn_var.allocate_fssh2()
        dyn_var.save_curr_dm_into_prev()


    #sys.exit(0)
    ham_aux = nHamiltonian(ham)


    #sys.exit(0)

    # Do the propagation
    for i in range(nsteps):

        #print(F"========== step {i} ===============")
        #========= Update variables, compute properties, and save ============    
        #dyn_var.update_amplitudes(dyn_params, ham);
        #dyn_var.update_density_matrix(dyn_params, ham, 1);
         
        # Energies 
        Ekin, Epot, Etot, dEkin, dEpot, dEtot = 0.0, 0.0, 0.0,  0.0, 0.0, 0.0
        Etherm, E_NHC = 0.0, 0.0

        """
        if is_nbra != 1:
            if force_method in [0, 1]:
                Ekin, Epot, Etot, dEkin, dEpot, dEtot = tsh_stat.compute_etot_tsh(ham, p, Cdia, Cadi, projectors, states, iM, rep_tdse)
            elif force_method in [2]:
                Ekin, Epot, Etot, dEkin, dEpot, dEtot = tsh_stat.compute_etot(ham, p, Cdia, Cadi, projectors, iM, rep_tdse)
    
            for bath in therm:
                Etherm += bath.energy()
            Etherm = Etherm / float(ntraj)
            E_NHC = Etot + Etherm

        """

        save.save_tsh_data_1234_new(_savers, dyn_params, i, dyn_var, ham)


        #============ Propagate ===========        
        index = i + icond
        if nfiles > 0:
            index = index % nfiles
        model_params.update({"timestep":index})

        compute_dynamics(dyn_var, dyn_params, ham, ham_aux, compute_model, model_params, rnd, therm);

        #dyn_var.update_amplitudes(dyn_params, ham);
        #dyn_var.update_density_matrix(dyn_params, ham, 1);

            
        if _savers["txt_saver"]!=None:
            _savers["txt_saver"].save_data_txt( F"{prefix}", properties_to_save, "a", i)

        if _savers["txt2_saver"]!=None:
            _savers["txt2_saver"].save_data_txt( F"{prefix2}", properties_to_save, "a", 0)


    if _savers["mem_saver"]!=None:
        _savers["mem_saver"].save_data( F"{prefix}/mem_data.hdf", properties_to_save, "w")
        return _savers["mem_saver"]




#def generic_recipe(q, p, iM, _dyn_params, compute_model, _model_params, _init_elec, rnd):
def generic_recipe(_dyn_params, compute_model, _model_params,_init_elec, _init_nucl, rnd):
    """
    This function initializes electronic variables and Hamiltonians for a given set of 
    nuclear variables (those are needed, not initialized here) and then runs the calculations 
    on multiple trajectories

    Args:
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

        _init_nucl (dict): controls how to sample nuclear DOFs from what the user provides. Can contain:
          * **_init_nucl["init_type"]** (int) : the type of the sampling:
            - 0 : Coords and momenta are set exactly to the given value
            - 1 : Coords are set, momenta are sampled
            - 2 : Coords are sampled, momenta are set
            - 3 : Both coords and momenta are samples [ default ]
          * **_init_nucl["ndof"]** (int) : the number of nuclear DOFs [default: 1]
          * **_init_nucl["q"]** ( list of `ndof` doubles ): average coordinate of all trajectories, for each nuclear DOF, 
                  a.u. of lenth [ default: [0.0] ]
          * **_init_nucl["p"]** ( list of `ndof` doubles ): average momentum of all trajectories, for each nuclear DOF,
                  a.u. of lenth [ default: [0.0] ]
          * **_init_nucl["mass"]** (list of `ndof` doubles): the mass for each nuclear DOF, in a.u. of mass [ default: [2000.0] ]
          * **_init_nucl["force_constant"]** (list of `ndof` doubles) : the force constant of the harmonic oscillator in each 
                  nuclear DOF, determins the width of the coordinate and momentum samplings [ default: [0.001] ]
 
        rnd ( Random ): random numbers generator object - controls the initialization of nuclear DOFs

    Returns:
        tuple : tuple returned by the :func:`libra_py.dynamics.run_dynamics` function



    """

    model_params = dict(_model_params)
    dyn_params = dict(_dyn_params)
    init_elec = dict(_init_elec)    
    init_nucl = dict(_init_nucl)

    comn.check_input( model_params, {  }, [ "model0" ] )
    comn.check_input( dyn_params, { "rep_tdse":1, "is_nbra":0, "ntraj":1 }, [  ] )
    comn.check_input( init_nucl, {"init_type":3, "ndof":1, "q":[0.0], "p":[0.0], "mass":[2000.0], "force_constant":[0.01] }, [] )


    is_nbra =  dyn_params["isNBRA"]
    init_elec.update({"is_nbra":is_nbra})
    comn.check_input( init_elec, {  "ndia":1, "nadi":1 }, [  ] )


    # Internal parameters
    ndia = init_elec["ndia"]
    nadi = init_elec["nadi"]
    ndof = len(init_nucl["mass"])
    ntraj = dyn_params["ntraj"]


    # Setup the dynamical variables object
    dyn_var = dyn_variables(ndia, nadi, ndof, ntraj)

    #sys.exit(0)

    # Initialize nuclear variables 
    dyn_var.init_nuclear_dyn_var(init_nucl, rnd)

    #sys.exit(0)

    # Initialize electronic variables
    dyn_var.init_amplitudes(init_elec, rnd)
    dyn_var.init_density_matrix(init_elec)
    dyn_var.init_active_states(init_elec, rnd)

    #sys.exit(0)

    # Setup the hierarchy of Hamiltonians 
    ham = nHamiltonian(ndia, nadi, ndof)
    if is_nbra == 1:
        ham.add_new_children(ndia, nadi, ndof, 1)
    else:
        ham.add_new_children(ndia, nadi, ndof, ntraj)
    ham.init_all(2,1)     


    #sys.exit(0)

    # Compute internals of the Hamiltonian objects
    model_params1 = dict(model_params)
    model_params1.update({"model":model_params["model0"], "timestep":0})

    # We set up this temporary dict such that we we request diabatic Ham calculations
    # followed by the transformation to the adiabatic one - this is what we need to get
    # the transformation matrices to convert amplitudes between the representations
    dyn_params1 = dict(dyn_params)        

    if(dyn_params["ham_update_method"]==2):        
        pass
        #update_Hamiltonian_variables( dyn_params1, dyn_var, ham, ham, compute_model, model_params1, 0)
        #update_Hamiltonian_variables( dyn_params1, dyn_var, ham, ham, compute_model, model_params1, 1)
    else:
        dyn_params1.update({ "ham_update_method":1, "ham_transform_method":1 })
        #update_Hamiltonian_q( dyn_params1, dyn_var, ham, compute_model, model_params1)
        update_Hamiltonian_variables( dyn_params1, dyn_var, ham, ham, compute_model, model_params1, 0)
        update_Hamiltonian_variables( dyn_params1, dyn_var, ham, ham, compute_model, model_params1, 1)

    #sys.exit(0)

    # Update internal dynamical variables using the computed properties of the Hamiltonian objects
    # Set up the "rep_tdse" variable here to the representation that coinsides with the initial representation
    # of electronic variables - this will convert the amplitudes to the proper representation
    dyn_var.update_amplitudes( {"rep_tdse":init_elec["rep"] }, ham)
    dyn_var.update_density_matrix( dyn_params, ham, 1)

    #print("Initial adiabatic amplitudes")
    #dyn_var.get_ampl_adi().show_matrix()

    #print("Initial diabatic amplitudes")
    #dyn_var.get_ampl_dia().show_matrix()

    #print("Initial adiabatic DM")
    #dyn_var.get_dm_adi(0).show_matrix()

    #print("Initial diabatic DM")
    #dyn_var.get_dm_dia(0).show_matrix()

    #print("Active states")
    #print(Cpp2Py(dyn_var.act_states))


    # Finally, start the dynamics calculations
    res = run_dynamics(dyn_var, dyn_params, ham, compute_model, model_params, rnd)
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



