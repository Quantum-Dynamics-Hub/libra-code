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
.. module:: dynamics_recipes_Ehrenfest
   :platform: Unix, Windows
   :synopsis: This module implements a number of predefined simulations runs of the Ehrenfest dynamics
       Contains:
           * Ehrenfest_dia0_dia_diah(q, p, iM, _dyn_params, compute_model, _model_params, _init_elec, rnd)
           * Ehrenfest_adi0_dia_diah(q, p, iM, _dyn_params, compute_model, _model_params, _init_elec, rnd)
           * Ehrenfest_dia0_adi_adih(q, p, iM, _dyn_params, compute_model, _model_params, _init_elec, rnd)
           * Ehrenfest_adi0_adi_adih(q, p, iM, _dyn_params, compute_model, _model_params, _init_elec, rnd)
           * Ehrenfest_dia0_adi_diah(q, p, iM, _dyn_params, compute_model, _model_params, _init_elec, rnd)
           * Ehrenfest_adi0_adi_diah(q, p, iM, _dyn_params, compute_model, _model_params, _init_elec, rnd)

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
from . import compute



def Ehrenfest_dia0_dia_diah(q, p, iM, _dyn_params, compute_model, _model_params, _init_elec, rnd):
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


    return compute.generic_recipe(q, p, iM, dyn_params, compute_model, model_params, init_elec, rnd)        



def Ehrenfest_adi0_dia_diah(q, p, iM, _dyn_params, compute_model, _model_params, _init_elec, rnd):
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


    return compute.generic_recipe(q, p, iM, dyn_params, compute_model, model_params, init_elec, rnd)        



def Ehrenfest_dia0_adi_adih(q, p, iM, _dyn_params, compute_model, _model_params, _init_elec, rnd):
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

    return compute.generic_recipe(q, p, iM, dyn_params, compute_model, model_params, init_elec, rnd)        



def Ehrenfest_adi0_adi_adih(q, p, iM, _dyn_params, compute_model, _model_params, _init_elec, rnd):
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


    return compute.generic_recipe(q, p, iM, dyn_params, compute_model, model_params, init_elec, rnd)        






def Ehrenfest_dia0_adi_diah(q, p, iM, _dyn_params, compute_model, _model_params, _init_elec, rnd):
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

    return compute.generic_recipe(q, p, iM, dyn_params, compute_model, model_params, init_elec, rnd)        



def Ehrenfest_adi0_adi_diah(q, p, iM, _dyn_params, compute_model, _model_params, _init_elec, rnd):
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

    return compute.generic_recipe(q, p, iM, dyn_params, compute_model, model_params, init_elec, rnd)        



