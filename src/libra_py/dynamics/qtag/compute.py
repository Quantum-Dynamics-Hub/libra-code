#*********************************************************************************                     
#* Copyright (C) 2022 Matthew Dutra and Alexey V. Akimov                                                   
#*                                                                                                     
#* This file is distributed under the terms of the GNU General Public License                          
#* as published by the Free Software Foundation, either version 3 of                                   
#* the License, or (at your option) any later version.                                                 
#* See the file LICENSE in the root directory of this distribution   
#* or <http://www.gnu.org/licenses/>.          
#***********************************************************************************
"""
.. module:: compute
   :platform: Unix, Windows
   :synopsis: This module implements a wrapper function for doing QTAG dynamics
       List of functions:  
           * init_nuclear_dyn_var(Q, P, M, params, rnd)
           * init_electronic_dyn_var(params, isNBRA, rnd)
           * init_amplitudes(q, Cdia, Cadi, dyn_params, compute_model, model_params, transform_direction=0)

.. moduleauthor:: Matthew Dutra and Alexey V. Akimov

"""

__author__ = "Matthew Dutra, Alexey V. Akimov"
__copyright__ = "Copyright 2022 Matthew Dutra, Alexey V. Akimov"
__credits__ = ["Matthew Dutra, Alexey V. Akimov"]
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

from . import save
from . import qtag_prop

def qtag_pops(surf_ids, coeff, S, target_states):
    """Returns a single-surface population *n*, calculated from the single-surface basis coefficients *c* 
       and the appropriate overlap matrix *ov* as <c^T|ov|c>. Note that this is equivalent to the norm for 
       a single-surface system.

    Args:
        surf_ids (list): List containing the trajectory indices on various states.

        coeff ( CMATRIX(ntraj, 1) ): The GBF expansion coefficients (in the super-basis) 

        S (CMATRIX(ntraj, ntraj) ): Basis functions super-overlap

        target_states (list): List of states for which the norm should be calculated.

    Returns:
        pops (list of floats): Surface population. The imaginary part should be zero.

    """

    pops = []
    for n in target_states:
        # Extract indices of the trajectories that sit on a given state n
        traj_on_surf = [index for index, traj_id in enumerate(surf_ids) if traj_id == n]
        ntraj_on_surf = len(traj_on_surf)

        # Extract the overlap matrix of the GBFs on that surface
        ov_surf = CMATRIX(ntraj_on_surf,ntraj_on_surf)
        pop_submatrix(S,     ov_surf, traj_on_surf, traj_on_surf)

        # Extract the coefficients of the GBFs on that surface
        c_surf = CMATRIX(ntraj_on_surf,1)
        pop_submatrix(coeff, c_surf,  traj_on_surf, [0])

        pops.append((c_surf.H()*ov_surf*c_surf).get(0).real)

    return pops



def qtag_energy(coeff, H):
    """Returns the system energy *e*, calculated from the total basis coefficients *coeff* 
       and the total Hamiltonian *H* as <G|H|G>.
       E = C^+ * H * C

    Args:
        c (CMATRIX(ntraj, 1) ): The ntraj-by-1 complex matrix of basis coefficients.

        H (CMATRIX(ntraj, ntraj) ): The full system super-Hamiltonian (all trajectories and surfaces).

    Returns:
        float: Total energy of the system. The imaginary part should be zero.
    """


    e = (coeff.H() * H * coeff).get(0).real

    return e



def run_qtag(_q, _p, _alp, _s, _states, _coeff, _iM, _dyn_params, _compute_model, _model_params):
    """
    Args: 
        _q ( MATRIX(nnucl, ntraj) ): coordinates of the "classical" particles [units: Bohr]
        _p ( MATRIX(nnucl, ntraj) ): momenta of the "classical" particles [units: a.u. of momenta]
        _alp ( MATRIX(nnucl, ntraj) ): widths of the GBFs [ units: Bohr^-1 ]
        _s ( MATRIX(1, ntraj) ): phases of all GBFs [ units: no ]
        _states ( intList, or list of ntraj ints ): the quantum state of each trajectory
        _coeff ( CMATRIX(nstates, ntraj) ): amplitudes of the GBFs
        _iM ( MATRIX(nnucl, 1) ): masses of classical particles [units: a.u.^-1]
        _dyn_params ( dictionary ): parameters controlling the execution of the dynamics
            Can contain:

            * **q_update_method** (int): the parameter specifying whether positions should be updated (adaptable) or frozen  

                - 0: frozen
                - 1: adaptable [ default ]


            * **p_update_method** (int): the parameter specifying whether momenta should be updated (adaptable) or frozen  

                - 0: frozen
                - 1: adaptable [ default ]


            * **a_update_method** (int): the parameter specifying whether widths should be updated (adaptable) or frozen  

                - 0: frozen [ default ]
                - 1: adaptable


            * **s_update_method** (int): the parameter specifying whether phases should be updated (adaptable) or frozen 

                - 0: frozen [ default ]
                - 1: adaptable (not available yet)


            * **q_sync_method** (int): the parameter specifying whether positions should be synchronized on low pop surfaces  

               - 0: unsynchronized (leave at initial value)
               - 1: synchronized to most populated state [ default ]


            * **p_sync_method** (int): the parameter specifying whether momenta should be synchronized on low pop surfaces  

              - 0: unsynchronized (leave at initial value)
              - 1: synchronized to most populated state [ default ]


            * **a_sync_method** (int): the parameter specifying whether widths should be synchronized on low pop surfaces  

              - 0: unsynchronized (leave at initial value) [ default ]
              - 1: synchronized to most populated state 


            * **s_sync_method** (int): the parameter specifying whether phases should be synchronized on low pop surfaces  

              - 0: unsynchronized (leave at initial value) [ default ]


            * **decpl_den** (float): a critical population above which a surface can evolve trajectories independently  


            * **qtag_pot_approx_method** (int): the method for approximating potential matrix elements in a basis  

              - 0: BAT (Bra-ket Averaged Taylor expansion)
              - 1: LHA (Local Harmonic Approximation)
              - 2: LHAe (Local Harmonic Approximation w/ exact Gaussian coupling elements)
              - 3: BATe (Bra-ket Averaged Taylor expansion w/ exact Gaussian coupling elements)


            * **mom_calc_type** (int): how to modify the raw basis momenta for numerical stability  

              - 0: unmodified (note that these are frequently unstable!)
              - 1: linear fitting of single-surface momenta


            * **linfit_beta** (float): a threshold for the linear fitting algorithm when needed for momenta modification  

 
            * **prefix** (string): name of the directory where all the results will be stored [ default: "out" ]


            * **hdf5_output_level** (int): the level of output for HDF5 printing of data [ default: -1 ] 


            * **txt2_output_level** (int): the level of output for txt printing of data; 3 is the current maximum  [ default: 0 ]


            * **progress_frequency** (int): how often to print out the results [ default: 1 - every timestep] 


            * **properties_to_save** (list of strings): a list containing the desired output quantities 


            * **target_states** (list of ints): the indices of quantum states for which to compute populations


            * **dt** (float): the timestep between each iteration [ in a.u. of time ]


            * **nsteps** (int): the number of total iterations; this together with `dt` determines the simulated time 


        _compute_model ( PyObject ): the function that computes the Hamiltonian object
        _model_params ( dict ): parameters that are passed to the Hamiltonian-computing function



    """
    

    dyn_params = dict(_dyn_params)
    Q = MATRIX(_q)
    P = MATRIX(_p)
    A = MATRIX(_alp)
    S = MATRIX(_s)
    C = CMATRIX(_coeff)

    default_params = {
        "hdf5_output_level":-1, "prefix":"out", "use_compression":0, "compression_level":[0,0,0], 
        "mem_output_level":4, "txt2_output_level":0, "properties_to_save": [], "progress_frequency": 1,
        "target_states":[0]
    }
    critical_params = []
    comn.check_input(dyn_params, default_params, critical_params)


    ndof = Q.num_of_rows
    ntraj = Q.num_of_cols  # total number of trajectories 
    nstates = len( set(_states) )
    active_states = list(_states)

    #Rename variables locally for convenience...
    #states = sorted(_states)  # <---- TENTATIVELY comment
    #nstates = len(states)
    #active_state = dyn_params["active_state"]
    target_states = dyn_params["target_states"]
    nsteps = dyn_params["nsteps"]
    dt = dyn_params["dt"]
    iM = dyn_params["iM"]
    progress_frequency = dyn_params["progress_frequency"]
    prefix = dyn_params["prefix"]

    #if active_state not in states:
    #    sys.exit("ERROR: Initial state ('active_state') must be in state list ('states') in dyn_params!")
    
    #Initialize the basis parameters {q,p,a,s}...
    #ntraj, Q, P, A, S, active_states = qtag_init.initialize(dyn_params)    
    #print("Initialized "+str(ntraj)+" trajectories across "+str(nstates)+" states")

    #Create initial projection vector b...
    #bt, c0 = qtag_init.coeffs(dyn_params, qpas, active_state)

    #Initialize savers...
    properties_to_save = dyn_params['properties_to_save']
    _savers = save.init_qtag_savers(dyn_params, _model_params, nsteps, ntraj, ndof, nstates)

            
    #Start simulation and walltime variables...
    walltime_start = time.time()
    t=0.0
    
    #Initialize the Hamiltonian objects in the C++ component of Libra...
    ham = nHamiltonian(nstates, nstates, ndof)
    ham.add_new_children(nstates, nstates, ndof, ntraj)    
    ham.init_all(_model_params["deriv_lvl"],1)
    _model_params.update({"timestep":0})
    
    ovlp = CMATRIX(ntraj, ntraj)
    hmat = CMATRIX(ntraj, ntraj)
    
    Coeff = CMATRIX(nstates,ntraj)

    # AVA: originally we use sorted states for some reason
    qtag_hamiltonian_and_overlap(Q, P, A, S, Coeff, Py2Cpp_int(active_states), iM, ham, 
                                 _compute_model, _model_params, dyn_params, ovlp, hmat)
    
    #Run the dynamics...
    for step in range(nsteps):
        # Built-in function for propagation in the non-othogonal basis
        propagate_electronic(0.5*dt, C, hmat, ovlp)


        #Compute the new coefficient vector c_new...
        #ct_new = qtag_calc.basis_diag(ntraj, dt, hmat, ovlp, bt)
        
        #Calculate the total energy and surface populations...
        etot = qtag_energy(C, hmat)
        pops = qtag_pops(_states, C, ovlp, target_states)


        #Update the basis parameters according to the new wavefunction (ct_new)...
        Q, P, A, S, active_states, C = qtag_prop.propagate_basis(Q, P, A, S, active_states, C, dyn_params, pops)

        #Compute the Hamiltonian and overlap matrix elements using the C++ routine 'qtag_ham_and_ovlp'...
        qtag_hamiltonian_and_overlap(Q, P, A, S, Coeff, Py2Cpp_int(active_states), iM, ham, 
                                     _compute_model, _model_params, dyn_params, ovlp, hmat)
    
        # Built-in function for propagation in the non-othogonal basis
        propagate_electronic(0.5*dt, C, hmat, ovlp)


        #Output the energy and populations to the notebook for the user to see...
        if step % progress_frequency == 0:
            print(etot, pops)
        
        #Save the specified data...
        save.save_qtag_data(_savers, dyn_params, step+1, etot, 0, pops, C, Q, P, A, S)        
        if _savers["txt2_saver"]!=None:
            _savers["txt2_saver"].save_data_txt( F"{prefix}", properties_to_save, "a", 0)


    #Print the total simulation time...
    walltime_end = time.time()
    print("Total wall time: ",walltime_end-walltime_start)

    return Q, P, A, S, active_states, C
    