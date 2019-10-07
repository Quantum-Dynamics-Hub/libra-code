#*********************************************************************************                     
#* Copyright (C) 2019 Alexey V. Akimov                                                   
#*                                                                                                     
#* This file is distributed under the terms of the GNU General Public License                          
#* as published by the Free Software Foundation, either version 2 of                                   
#* the License, or (at your option) any later version.                                                 
#* See the file LICENSE in the root directory of this distribution   
#* or <http://www.gnu.org/licenses/>.          
#***********************************************************************************
"""
.. module:: tsh_algo1
   :platform: Unix, Windows
   :synopsis: This module implements the function for many-trajectories TSH hopping
.. moduleauthor:: Alexey V. Akimov

"""

__author__ = "Alexey V. Akimov"
__copyright__ = "Copyright 2019 Alexey V. Akimov"
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

import util.libutil as comn
from . import tsh
from . import tsh_stat
from . import units



def run_tsh(_q, _p, _iM, _Cdia, _Cadi, states, model_params, dyn_params, compute_model, rnd):
    """
  
    This function is similar to the one above, but it stores the computed properties in files
   
    Args: 
        _q ( MATRIX(nnucl, ntraj) ): coordinates of the "classical" particles [units: Bohr]
        _p ( MATRIX(nnucl, ntraj) ): momenta of the "classical" particles [units: a.u. of momenta]
        _iM ( MATRIX(nnucl, 1) ): masses of classical particles [units: a.u.^-1]
        _Cdia ( CMATRIX(ndia, ntraj) ): amplitudes of the diabatic basis states
        _Cadi ( CMATRIX(nadi, ntraj) ): amplitudes of the adiabatic basis states
        states ( intList, or list of ntraj ints ): the quantum state of each trajectory
        model_params ( dictionary ): contains the selection of a model and the parameters 
            for that model Hamiltonian


        dyn_params ( dictionary ): parameters controlling the execution of the dynamics
            Can contain:
      
            * **dyn_params["rep"]** ( int ): selects the representation in which nuclear/electronic (Ehrenfest core)
                dynamics is executed

                - 0: diabatic representation
                - 1: adiabatic representation [default: 1]

            * **dyn_params["ham_rep"]** (int): The representation of the Hamiltonian update: 

                - 0: diabatic
                - 1: adiabatic
                              
                This is the representation in which the computed properties are assumed to be
                For instance, we may have set it to 1, to read the adiabatic energies and couplings,
                to bypass the diabatic-to-adiabatic transformation, which may be useful in some atomistic
                calculations, or with the NBRA

            * **dyn_params["rep_sh"]** ( int ): selects the representation which is 
                used to perform surface hopping

                - 0: diabatic representation
                - 1: adiabatic representation [default: 1]

            * **dyn_params["rep_lz"]** ( int ): The representation to compute LZ probabilitieis:
 
                - 0: diabatic
                - 1: adiabatic 

            * **dyn_params["tsh_method"]** ( int ): Formula for computing SH probabilities: 
   
                - 0: FSSH
                - 1: GFSH 
                - 2: MSSH

            * **dyn_params["use_boltz_factor"]** ( int ): Whether to scale the SH probabilities
                by the Boltzmann factor: 
 
                - 0: do not scale
                - 1: scale

            * **dyn_params["Temperature"]** ( double ): Temperature of the system

            * **dyn_params["do_reverse"]** ( int ): what to do with velocities on frustrated hops:

                - 0: do not revert momenta at the frustrated hops
                - 1: do revert the momenta

            * **dyn_params["vel_rescale_opt"]** ( int ): How to rescale momenta if the hops are successful:

                - 0: rescale along the directions of derivative couplings
                - 1: rescale in the diabatic basis - don't care about the
                    velocity directions, just a uniform rescaling,
                - 2: do not rescale, as in the NBRA.

            * **dyn_params["dt"]** ( double ): the nuclear and electronic integration
                timestep [ units: a.u. of time, default: 41.0 ]

            * **dyn_params["do_phase_correction"]** ( int ): the algorithm to correct phases on adiabatic states
 
                - 0: no phase correction
                - 1: according to our phase correction algorithm [ default: 1 ]

            * **dyn_params["state_tracking_algo"]** ( int ): the algorithm to keep track of the states' identities
 
                - 0: no state tracking
                - 1: Sato
                - 2: using the mincost, Munkres-Kuhn [ default: 2 ]

            * **dyn_params["MK_alpha"]** ( double ): Munkres-Kuhn alpha 
                (selects the range of orbitals included in reordering) [default: 0.0]

            * **dyn_params["MK_verbosity"]** ( double ): Munkres-Kuhn verbosity: 

                - 0: no extra output [ default ]
                - 1: print details


            * **dyn_params["nsteps"]** ( int ): the number of NA-MD steps to do [ default: 1 ]

            * **dyn_params["tsh_version"]** ( int ): 
                - 1:  tsh1(q, p, iM,  Cdia, states, ham, compute_model, model_params, dyn_params, rnd)
                - 2:  tsh1b(q, p, iM,  Cdia, states, ham, compute_model, model_params, dyn_params, rnd) [ default ]

            * **dyn_params["output_level"]** ( int ): controls what info to return as the result of this function

                - -1: all return values are None
                - 0: also, all return values are none 
                - 1: 1-D info - energies, SE and SH populations averaged over ensemble vs. time [ default ]
                - 2: 2-D info - coordinates, momenta, SE amplitudes, and SH states populations 
                    for each individual trajectory vs. time 
                - 3: 3-D info - St, Hvib_adi, Hvib_dia matrices (energies, couplings, etc.) for each
                    individual trajectory vs. time 


            * **dyn_params["file_output_level"]** ( int ): controls what info to print out into files

                - -1: print nothing at all 
                - 0: 0-D info - just the parameters of the simulation [ default ]
                - 1: 1-D info - energies, SE and SH populations averaged over ensemble vs. time 
                - 2: 2-D info - coordinates, momenta, SE amplitudes, and SH states populations 
                    for each individual trajectory vs. time 
                - 3: 3-D info - St, Hvib_adi, Hvib_dia matrices (energies, couplings, etc.) for each
                    individual trajectory vs. time 
            

        compute_model ( PyObject ): the pointer to the Python function that performs the Hamiltonian calculations
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


        Note: the elements are None, if they are excluded by the input optinos

    """

        
    # Create copies of the input dynamical variables, so we could run several such
    # functions with the same input variables without worries that they will be altered
    # inside of each other

    q = MATRIX(_q)
    p = MATRIX(_p)
    iM = MATRIX(_iM)
    Cdia = CMATRIX(_Cdia)
    Cadi = CMATRIX(_Cadi)


    # Parameters and dimensions
    critical_params = [ "prefix" ] 
    default_params = { "rep":1, "nsteps":1, "dt":1.0*units.fs2au, "do_phase_correction":1, "state_tracking_algo":2,
                       "MK_alpha":0.0, "MK_verbosity":0, "file_output_level":1, "tsh_version":2, "tol":1e-3 }
    comn.check_input(dyn_params, default_params, critical_params)
        
    prefix = dyn_params["prefix"]    
    rep = dyn_params["rep"]
    nsteps = dyn_params["nsteps"]
    tsh_version = dyn_params["tsh_version"]
    dt = dyn_params["dt"]
    tol = dyn_params["tol"]
    output_level = dyn_params["output_level"]
    file_output_level = dyn_params["file_output_level"]
    
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
    obs_hvib_adi, obs_hvib_dia, obs_St = None, None, None
         
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
        for tr in xrange(ntraj):
            f = open("%s/Hvib_adi_%i.txt" % (prefix, tr), "w"); f.close()   # vibronic Hamiltonians for each trajectory in adiabatic rep. 
            f = open("%s/Hvib_dia_%i.txt" % (prefix, tr), "w"); f.close()   # vibronic Hamiltonians for each trajectory in diabatic rep. 
            f = open("%s/St_%i.txt" % (prefix, tr), "w"); f.close()   # time-overlaps along each trajectory 
        

    # ======= Hierarchy of Hamiltonians =======
    ham = nHamiltonian(ndia, nadi, nnucl)
    ham.init_all(2)
    ham.phase_corr_ovlp_tol = tol

    ham1 = [] 
    for tr in range(ntraj):
        ham1.append( nHamiltonian(ndia, nadi, nnucl) )        
        ham1[tr].init_all(2)
        ham1[tr].phase_corr_ovlp_tol = tol
        ham.add_child(ham1[tr])
        

    # Initial calculations
    ham.compute_diabatic(compute_model, q, model_params, 1)
    ham.compute_adiabatic(1, 1); 
    ham.ampl_adi2dia(Cdia, Cadi, 0, 1)
    
    U = []
    for tr in range(ntraj):
        U.append(ham1[tr].get_basis_transform())

    if rep==0:
        ham.compute_nac_dia(p, iM, 0, 1);  
        ham.compute_hvib_dia(1); 
    elif rep==1:
        ham.compute_nac_adi(p, iM, 0, 1);
        ham.compute_hvib_adi(1); 


                
    # Do the propagation
    for i in range(nsteps):
    
        #============ Compute and output properties ===========        
        if rep==0:
            ham.ampl_dia2adi(Cdia, Cadi, 0, 1)
        elif rep==1:
            ham.ampl_adi2dia(Cdia, Cadi, 0, 1)

        dm_dia, dm_adi = tsh_stat.compute_dm(ham, Cdia, Cadi, rep, 1)        
        Ekin, Epot, Etot, dEkin, dEpot, dEtot = tsh_stat.compute_etot_tsh(ham, p, Cdia, Cadi, states, iM, rep)
        pops = tsh_stat.compute_sh_statistics(nadi, states)
        
            
        # Memory output
        if output_level>=1:
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

        # File output
        if file_output_level >= 1:
            add_doublelist2file("%s/energies.txt" % (prefix), i*dt, [Ekin, Epot, Etot, dEkin, dEpot, dEtot] )
            add_cmatrix2file("%s/D_adi.txt" % (prefix), i*dt, dm_adi )
            add_cmatrix2file("%s/D_dia.txt" % (prefix), i*dt, dm_dia )
            add_matrix2file("%s/SH_pop.txt" % (prefix), i*dt, pops )


        # Memory output
        if output_level>=2:
            obs_q.append(MATRIX(q))
            obs_p.append(MATRIX(p))        
            obs_Cadi.append(CMATRIX(Cadi))
            obs_Cdia.append(CMATRIX(Cdia))
            obs_states.append(list(states))

        # File output
        if file_output_level >= 2:        
            add_matrix2file("%s/q.txt" % (prefix), i*dt, q )
            add_matrix2file("%s/p.txt" % (prefix), i*dt, p )
            add_cmatrix2file("%s/C_adi.txt" % (prefix), i*dt, Cadi )
            add_cmatrix2file("%s/C_dia.txt" % (prefix), i*dt, Cdia )
            add_intlist2file("%s/states.txt" % (prefix), i*dt, states )        
    
        for tr in range(ntraj):            
            x = ham1[tr].get_basis_transform()
            St = U[tr] * x.H()             
            U[tr] = CMATRIX(x)

            # Memory output            
            if output_level >= 3:
                obs_hvib_adi[tr].append( CMATRIX(ham1[tr].get_hvib_adi()) )
                obs_hvib_dia[tr].append( CMATRIX(ham1[tr].get_hvib_dia()) )
                obs_St[tr].append( CMATRIX(St) )

            # File output
            if file_output_level >= 3:
                add_cmatrix2file("%s/St_%i.txt" % (prefix, tr), i*dt, St)
                add_cmatrix2file("%s/Hvib_adi_%i.txt" % (prefix, tr), i*dt, ham1[tr].get_hvib_adi())
                add_cmatrix2file("%s/Hvib_dia_%i.txt" % (prefix, tr), i*dt, ham1[tr].get_hvib_dia())


        #============ Propagate ===========        
        if rep==0:
            if tsh_version==1:
                tsh1(q, p, iM,  Cdia, states, ham, compute_model, model_params, dyn_params, rnd)
            elif tsh_version==2:
                tsh1b(q, p, iM,  Cdia, states, ham, compute_model, model_params, dyn_params, rnd)
        elif rep==1:
            if tsh_version==1:
                tsh1(q, p, iM,  Cadi, states, ham, compute_model, model_params, dyn_params, rnd)
            elif tsh_version==2:
                tsh1b(q, p, iM,  Cadi, states, ham, compute_model, model_params, dyn_params, rnd)            




    return obs_T, obs_q, obs_p, obs_Ekin, obs_Epot, obs_Etot, obs_dEkin, obs_dEpot, obs_dEtot, \
           obs_Cadi, obs_Cdia, obs_dm_adi, obs_dm_dia, obs_pop, obs_states, obs_hvib_adi, obs_hvib_dia, obs_St



def probabilities_1D_scattering(q, states, nst, params):
    """Computes the scattering probabilities in 1D

    Args:
        _q ( MATRIX(nnucl, ntraj) ): coordinates of the "classical" particles [units: Bohr]
        states ( intList, or list of ntraj ints ): the quantum state of each trajectory
        nst ( int ): the number of possible quantum states in the problem
        params ( dictionary ): parameters of the simulation, should contain
 
            * **params["act_dof"]** ( int ): index of the nuclear DOF that is considered active (scattering coord)
            * **params["left_boundary"] ( double ): the beginning of the reflected particles counter [units: Bohr]
            * **params["right_boundary"] ( double ): the beginning of the transmitted particles counter [units: Bohr]

    Returns:
        tuple: ( pop_refl, pop_transm ): where

            * pop_refl ( MATRIX(nst, 1) ): probabilities of reflection on each state
            * pop_transm ( MATRIX(nst, 1) ): probabilities of transmission on each state

    """

    critical_params = [  ] 
    default_params = {"act_dof":0, "left_boundary":-10.0, "right_boundary":10.0 }
    comn.check_input(params, default_params, critical_params)


    act_dof = params["act_dof"]
    left_boundary = params["left_boundary"]
    right_boundary = params["right_boundary"]


    ntraj = len(states)

    pop_transm = MATRIX(nst, 1)  # transmitted
    pop_refl = MATRIX(nst, 1)    # reflected

    ntransm, nrefl = 0.0, 0.0
    for traj in range(0,ntraj):

        if q.get(act_dof, traj) < left_boundary:
            pop_refl.add(states[traj], 0, 1.0)
            nrefl += 1.0

        if q.get(act_dof, traj) > right_boundary:
            pop_transm.add(states[traj], 0, 1.0)
            ntransm += 1.0         

    ntot = ntransm + nrefl 
    if ntot > 0.0:
        pop_transm = pop_transm / ntot
        pop_refl = pop_refl / ntot

    return pop_refl, pop_transm
