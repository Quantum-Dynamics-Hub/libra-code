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

import tsh
import tsh_stat


def run_tsh(_q, _p, _iM, _Cdia, _Cadi, states, model_params, dyn_params, rnd):
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

        rnd ( Random ): random numbers generator object

    Returns:
        tuple: ()
              
    """

    
    obs_T = [] # time
    obs_q = [] # coordinate of the first DOF
    obs_p = [] # momentum of the first DOF
    obs_Ekin = []  # average kinetic energy 
    obs_Epot = []  # average potential energy 
    obs_Etot = []  # average total energy 
    obs_dEkin = []  # kinetic energy fluctuation
    obs_dEpot = []  # potential energy fluctuation
    obs_dEtot = []  # total energy fluctuation
    obs_dm_adi00 = []  # average SE-based population of the state 0 in adiabatic basis
    obs_dm_adi11 = []  # average SE-based population of the state 1 in adiabatic basis
    obs_dm_dia00 = []  # average SE-based population of the state 0 in diabatic basis
    obs_dm_dia11 = []  # average SE-based population of the state 1 in diabatic basis
    obs_pop00 = []  # average SH-based population of the state 0 in diabatic basis
    obs_pop01 = []  # average ??
    obs_ind = []  # ??
    
    
    # Create copies of the input dynamical variables, so we could run several run_test 
    # functions with the same input variables without worries that they will be altered
    # inside of the run_test

    q = MATRIX(_q)
    p = MATRIX(_p)
    iM = MATRIX(_iM)
    Cdia = CMATRIX(_Cdia)
    Cadi = CMATRIX(_Cadi)
    
    rep = dyn_params["rep"]
    dt = dyn_params["dt"]
    nsteps = dyn_params["nsteps"]


    ndia = Cdia.num_of_rows
    nadi = Cadi.num_of_rows
    nnucl= q.num_of_rows
    ntraj= q.num_of_cols

    # ======= Hierarchy of Hamiltonians =======
    ham = nHamiltonian(ndia, nadi, nnucl)
    ham.init_all(2)
    #print "id=", ham.id, " level=", ham.level

    ham1 = [] 
    for tr in xrange(ntraj):
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


#    sys.exit(0)
    Ekin, Epot, Etot, dEkin, dEpot, dEtot = tsh_stat.compute_etot_tsh(ham, p, Cdia, Cadi, states, iM, rep) 
    #print Ekin, Epot, Etot, dEkin, dEpot, dEtot


    
    # Do the propagation
    for i in xrange(nsteps):

        if rep==0:
            tsh1(dt, q, p, iM,  Cdia, states, ham, compute_model, model_params, dyn_params, rnd)
        elif rep==1:
            tsh1(dt, q, p, iM,  Cadi, states, ham, compute_model, model_params, dyn_params, rnd, 1, 1)

        #=========== Properties ==========
        if rep==0:
            ham.ampl_dia2adi(Cdia, Cadi, 0, 1)
        elif rep==1:
            ham.ampl_adi2dia(Cdia, Cadi, 0, 1)

        dm_dia, dm_adi = tsh_stat.compute_dm(ham, Cdia, Cadi, rep, 1)
        Ekin, Epot, Etot, dEkin, dEpot, dEtot = tsh_stat.compute_etot_tsh(ham, p, Cdia, Cadi, states, iM, rep)
        pops = tsh_stat.compute_sh_statistics(nadi, states)
    

        ind = 0.0
        for tr in xrange(ntraj):
            ind = ind + ham1[tr].get_ordering_adi()[0]
        ind = ind/float(ntraj)
        
        
        obs_T.append(i*dt) 
        obs_q.append(q)
        obs_p.append(p)
        obs_Ekin.append(Ekin)
        obs_Epot.append(Epot)
        obs_Etot.append(Etot)
        obs_dEkin.append(dEkin)
        obs_dEpot.append(dEpot)
        obs_dEtot.append(dEtot)
        obs_Cadi.append(Cadi)
        obs_Cdia.append(Cdia)
        obs_dm_adi.append(dm_adi)
        obs_dm_dia.append(dm_dia)
        obs_pops.append(pops)
        obs_ind.append(ind)
        
    return obs_T, obs_q, obs_p, obs_dm_adi00, obs_dm_adi11, obs_dm_dia00, obs_dm_dia11, obs_pop00, obs_pop01, obs_ind
