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
.. module:: tsh_algo2
   :platform: Unix, Windows
   :synopsis: This module implements the function for many-trajectories TSH hopping
       and also applied a thermostat
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
from . import tsh
from . import tsh_stat
from . import units


def run_tsh(_q, _p, _iM, _Cdia, _Cadi, states, model_params, dyn_params, compute_model, rnd):
    """

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

            * **dyn_params["rep_sh"]** ( int ): selects the representation which is 
                used to perform surface hopping

                - 0: diabatic representation
                - 1: adiabatic representation [default: 1]

            * **dyn_params["nsteps"]** ( int ): the number of NA-MD steps to do [ default: 1 ]

            * **dyn_params["dt"]** ( double ): the nuclear and electronic integration
                timestep [ units: a.u. of time, default: 1.0 ]
 
            * **dyn_params["BATH_params"]** ( Dictionary ): 

        compute_model ( PyObject ): the pointer to the Python function that performs the Hamiltonian calculations
        rnd ( Random ): random numbers generator object

    Returns:
        tuple: ( obs_T, obs_q, obs_p, obs_Ekin, obs_Epot, obs_Etot, obs_dEkin, obs_dEpot, obs_dEtot, obs_Cadi, obs_Cdia, obs_dm_adi, obs_dm_dia, obs_pops ), where

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
              
    """

    
    obs_T = [] # time
    obs_q = [] # coordinates of all trajectories
    obs_p = [] # momenta of all trajectories
    obs_Ekin = []  # average kinetic energy 
    obs_Epot = []  # average potential energy 
    obs_Etot = []  # average total energy 
    obs_dEkin = []  # kinetic energy fluctuation
    obs_dEpot = []  # potential energy fluctuation
    obs_dEtot = []  # total energy fluctuation
    obs_Cadi = []  # average TD-SE amplitudes in the adiabatic basis
    obs_Cdia = []  # average TD-SE amplitudes in the diabatic basis
    obs_dm_adi = []  # average SE-based density matrix in adiabatic basis
    obs_dm_dia = []  # average SE-based density matrix in diabatic basis
    obs_pop = []  # average SH-based populations adiabatic basis
    obs_states = []  # indices of the quantum states of each trajectory
    #obs_ind = []  # ??
    
    
    # Create copies of the input dynamical variables, so we could run several run_test 
    # functions with the same input variables without worries that they will be altered
    # inside of the run_test

    q = MATRIX(_q)
    p = MATRIX(_p)
    iM = MATRIX(_iM)
    Cdia = CMATRIX(_Cdia)
    Cadi = CMATRIX(_Cadi)

    # Parameters and dimensions
    critical_params = [ ] 
    default_params = { "rep":1, "nsteps":1, "dt":1.0*units.fs2au, "do_phase_correction":1, "state_tracking_algo":2,
                       "MK_alpha":0.0, "MK_verbosity":0, "tsh_version":1 }
    comn.check_input(dyn_params, default_params, critical_params)
    

    rep = dyn_params["rep"]
    nsteps = dyn_params["nsteps"]
    dt = dyn_params["dt"]
    tsh_version = dyn_params["tsh_version"]
    BATH_params = dyn_params["BATH_params"]


    ndia = Cdia.num_of_rows
    nadi = Cadi.num_of_rows
    nnucl= q.num_of_rows
    ntraj= q.num_of_cols

    # ======= Hierarchy of Hamiltonians =======
    ham = nHamiltonian(ndia, nadi, nnucl)
    ham.init_all(2)

    ham1 = [] 
    for tr in range(0,ntraj):
        ham1.append( nHamiltonian(ndia, nadi, nnucl) )        
        ham1[tr].init_all(2)
        ham.add_child(ham1[tr])
        

    # Initial calculations
    ham.compute_diabatic(compute_model, q, model_params, 1)
    ham.compute_adiabatic(1, 1); 
    ham.ampl_adi2dia(Cdia, Cadi, 0, 1)


    if rep==0:
        ham.compute_nac_dia(p, iM, 0, 1);  
        ham.compute_hvib_dia(1); 
    elif rep==1:
        ham.compute_nac_adi(p, iM, 0, 1);
        ham.compute_hvib_adi(1); 


    Ekin, Epot, Etot, dEkin, dEpot, dEtot = tsh_stat.compute_etot_tsh(ham, p, Cdia, Cadi, states, iM, rep) 

    # Thermostats
    therms = []

    for tr in range(0,ntraj):
        therms.append( Thermostat( BATH_params ) )
        therms[tr].set_Nf_t(nnucl)
        therms[tr].set_Nf_r(0)
        therms[tr].set_Nf_b(0)
        therms[tr].init_nhc()   

        
    # Do the propagation
    for i in range(0,nsteps):

        for tr in range(0,ntraj):
            # Rescale momenta
            scl = therms[tr].vel_scale(0.25*dt)
            p.scale(-1, tr, scl)

            # Update thermostat variables
            ekin = 0.0
            for dof in range(0,nnucl):
                ekin = ekin + 0.5 * iM.get(dof, 0) * p.get(dof, tr)**2

            therms[tr].propagate_nhc(0.5*dt, ekin, 0.0, 0.0)

            # Rescale momenta
            scl = therms[tr].vel_scale(0.25*dt)
            p.scale(-1, tr, scl)


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


        for tr in range(0,ntraj):
            # Rescale momenta
            scl = therms[tr].vel_scale(0.25*dt)
            p.scale(-1, tr, scl)

            # Update thermostat variables
            ekin = 0.0
            for dof in range(0,nnucl):
                ekin = ekin + 0.5 * iM.get(dof, 0) * p.get(dof, tr)**2

            therms[tr].propagate_nhc(0.5*dt, ekin, 0.0, 0.0)

            # Rescale momenta
            scl = therms[tr].vel_scale(0.25*dt)
            p.scale(-1, tr, scl)




        #=========== Properties ==========
        if rep==0:
            ham.ampl_dia2adi(Cdia, Cadi, 0, 1)
        elif rep==1:
            ham.ampl_adi2dia(Cdia, Cadi, 0, 1)

        dm_dia, dm_adi = tsh_stat.compute_dm(ham, Cdia, Cadi, rep, 1)
        Ekin, Epot, Etot, dEkin, dEpot, dEtot = tsh_stat.compute_etot_tsh(ham, p, Cdia, Cadi, states, iM, rep)
        pops = tsh_stat.compute_sh_statistics(nadi, states)
    

        ind = 0.0
        for tr in range(0,ntraj):
            ind = ind + ham1[tr].get_ordering_adi()[0]
        ind = ind/float(ntraj)
        
        
        obs_T.append(i*dt) 
        obs_q.append(MATRIX(q))
        obs_p.append(MATRIX(p))
        obs_Ekin.append(Ekin)
        obs_Epot.append(Epot)
        obs_Etot.append(Etot)
        obs_dEkin.append(dEkin)
        obs_dEpot.append(dEpot)
        obs_dEtot.append(dEtot)
        obs_Cadi.append(CMATRIX(Cadi))
        obs_Cdia.append(CMATRIX(Cdia))
        obs_dm_adi.append(CMATRIX(dm_adi))
        obs_dm_dia.append(CMATRIX(dm_dia))
        obs_pop.append(MATRIX(pops))
        obs_states.append(list(states))
        #obs_ind.append(ind)
        
    return obs_T, obs_q, obs_p, obs_Ekin, obs_Epot, obs_Etot, obs_dEkin, obs_dEpot, obs_dEtot, obs_Cadi, obs_Cdia, obs_dm_adi, obs_dm_dia, obs_pop, obs_states


