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
.. module:: dynamics_recipes
   :platform: Unix, Windows
   :synopsis: This module implements a number of predefined execution types
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
from . import units
from . import data_outs
from . import dynamics


def generic_recipe(q, p, iM, compute_model, init_elec, dyn_params, model_params, rnd):

    # Internal parameters
    nnucl, ntraj = q.num_of_rows, q.num_of_cols

    # Initialize electronic variables - either diabatic or adiabatic
    Cdia, Cadi, states = dynamics.init_electronic_dyn_var(init_elec, rnd)
    ndia, nadi = Cdia.num_of_rows, Cadi.num_of_rows

    # Compute the diabatic-to-adiabatic transformation matrices
    ham = nHamiltonian(ndia, nadi, nnucl)
    ham.add_new_children(ndia, nadi, nnucl, ntraj)
    ham.init_all(2,1)     

    mdl_prms = dict(model_params)
    mdl_prms.update({"model":1})
    update_Hamiltonian_q({"rep_tdse":1, "rep_ham":0}, q, ham, compute_model, mdl_prms )

    if init_elec["rep"]==0:
        Cadi = transform_amplitudes(0, 1, Cdia, ham) 
    elif init_elec["rep"]==1:    
        Cdia = transform_amplitudes(1, 0, Cdia, ham) 

    res = dynamics.run_dynamics(q, p, iM, Cdia, Cadi, states, model_params, dyn_params, compute_model, rnd)

    return res




def Ehrenfest_dia0_dia_diah(q, p, iM, compute_model, _init_elec, _dyn_params, _model_params, rnd):
    """
    dynamics - Ehrenfest
    initial conditions - diabatic
    propagation: TDSE/Nuclear - diabatic
    compute_model - diabatic
    """

    init_elec = dict(_init_elec)
    dyn_params = dict(_dyn_params)
    model_params = dict(_model_params)


    init_elec.update({"rep":0})
    dyn_params.update({"force_method":2, "rep_ham":0, "rep_force":0, "rep_tdse":0, "nac_update_method":0, "tsh_method":-1})
    dyn_params.update({"do_phase_correction": 0, "state_tracking_algo":0})

    return generic_recipe(q, p, iM, compute_model, init_elec, dyn_params, model_params, rnd)



def Ehrenfest_adi0_dia_diah(q, p, iM, compute_model, _init_elec, _dyn_params, _model_params, rnd):
    """
    dynamics - Ehrenfest
    initial conditions - adiabatic
    propagation: TDSE/Nuclear - diabatic
    compute_model - diabatic
    """

    init_elec = dict(_init_elec)
    dyn_params = dict(_dyn_params)
    model_params = dict(_model_params)


    init_elec.update({"rep":1})
    dyn_params.update({"force_method":2, "rep_ham":0, "rep_force":0, "rep_tdse":0, "nac_update_method":0, "tsh_method":-1})
    dyn_params.update({"do_phase_correction": 0, "state_tracking_algo":0})

    return generic_recipe(q, p, iM, compute_model, init_elec, dyn_params, model_params, rnd)




def Ehrenfest_dia0_adi_adih(q, p, iM, compute_model, _init_elec, _dyn_params, _model_params, rnd):
    """
    dynamics - Ehrenfest
    initial conditions - diabatic
    propagation: TDSE/Nuclear - adiabatic
    compute_model - adiabatic
    """

    init_elec = dict(_init_elec)
    dyn_params = dict(_dyn_params)
    model_params = dict(_model_params)


    init_elec.update({"rep":0})
    dyn_params.update({"force_method":2, "rep_ham":1, "rep_force":1, "rep_tdse":1, "nac_update_method":1, "tsh_method":-1})

    return generic_recipe(q, p, iM, compute_model, init_elec, dyn_params, model_params, rnd)



def Ehrenfest_adi0_adi_adih(q, p, iM, compute_model, _init_elec, _dyn_params, _model_params, rnd):
    """
    dynamics - Ehrenfest
    initial conditions - adiabatic
    propagation: TDSE/Nuclear - adiabatic
    compute_model - adiabatic
    """

    init_elec = dict(_init_elec)
    dyn_params = dict(_dyn_params)
    model_params = dict(_model_params)


    init_elec.update({"rep":1})
    dyn_params.update({"force_method":2, "rep_ham":1, "rep_force":1, "rep_tdse":1, "nac_update_method":1, "tsh_method":-1})

    return generic_recipe(q, p, iM, compute_model, init_elec, dyn_params, model_params, rnd)






def Ehrenfest_dia0_adi_diah(q, p, iM, compute_model, _init_elec, _dyn_params, _model_params, rnd):
    """
    dynamics - Ehrenfest
    initial conditions - diabatic
    propagation: TDSE/Nuclear - adiabatic
    compute_model - diabatic
    """

    init_elec = dict(_init_elec)
    dyn_params = dict(_dyn_params)
    model_params = dict(_model_params)


    init_elec.update({"rep":0})
    dyn_params.update({"force_method":2, "rep_ham":0, "rep_force":1, "rep_tdse":1, "nac_update_method":1, "tsh_method":-1})

    return generic_recipe(q, p, iM, compute_model, init_elec, dyn_params, model_params, rnd)



def Ehrenfest_adi0_adi_diah(q, p, iM, compute_model, _init_elec, _dyn_params, _model_params, rnd):
    """
    dynamics - Ehrenfest
    initial conditions - adiabatic
    propagation: TDSE/Nuclear - adiabatic
    compute_model - diabatic
    """

    init_elec = dict(_init_elec)
    dyn_params = dict(_dyn_params)
    model_params = dict(_model_params)


    init_elec.update({"rep":1})
    dyn_params.update({"force_method":2, "rep_ham":0, "rep_force":1, "rep_tdse":1, "nac_update_method":1, "tsh_method":-1})

    return generic_recipe(q, p, iM, compute_model, init_elec, dyn_params, model_params, rnd)


