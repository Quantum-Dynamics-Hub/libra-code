#***********************************************************
# * Copyright (C) 2017-2018 Brendan Smith, Wei Li and Alexey V. Akimov
# * This file is distributed under the terms of the
# * GNU General Public License as published by the
# * Free Software Foundation; either version 3 of the
# * License, or (at your option) any later version.
# * http://www.gnu.org/copyleft/gpl.txt
#***********************************************************/

import os
import sys

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *

from utils import*

def compute_properties_dia_gamma(params, es_curr, es_next, curr_index):
    """
    This fucntions computes the properties needed to construct the
    vibrionic Hamiltonian, and then computes it.
    \param[in] params A dictionary containing important simulation parameters
    \param[in] es_curr A dictionary containing the data for the g-vectors and pw coefficients
                     for the current timestep
    \param[in] es_next A dictionary containing the data for the g-vectors and pw coefficients 
                     for the next timestep
    \param[in] curr_index This is index represents the current time step
  
    Returns: The vibrionic Hamiltonian in Ha = a.u. of energy    
    """

    Ca_curr = es_curr["Coeff_dia"]
    Ca_next = es_next["Coeff_dia"]
    Ea_curr = es_curr["E_dia"]
    Ea_next = es_next["E_dia"]

    rd = get_value(params,"rd",os.getcwd()+"../../res","s")   # of where the files will be printed out
    dt = get_value(params,"dt","41.34145","f") # time step in a.u - rescale NAC if actual dt is different

    S  = Ca_curr.H() * Ca_curr
    St = Ca_curr.H() * Ca_next
    nac = (0.5/dt)*(St - St.H())
    edia = 0.5*(Ea_curr + Ea_next)
    hvib = edia - 1.0j * nac

    os.system("mkdir %s" % rd)
    #========== Print out ================
    #S.real().show_matrix("%s/S_dia_ks_%d_re" % (rd, curr_index))

    St.real().show_matrix("%s/St_dia_ks_%d_re" % (rd, curr_index))
    St.imag().show_matrix("%s/St_dia_ks_%d_im" % (rd, curr_index))

    #nac.real().show_matrix("%s/nac_dia_ks_%d_re" % (rd, curr_index))
    #nac.imag().show_matrix("%s/nac_dia_ks_%d_im" % (rd, curr_index))

    edia.real().show_matrix("%s/E_dia_ks_%d_re" % (rd, curr_index))
    #edia.imag().show_matrix("%s/E_dia_ks_%d_im" % (rd, curr_index))

    hvib.real().show_matrix("%s/hvib_dia_%d_re" % (rd, curr_index))
    hvib.imag().show_matrix("%s/hvib_dia_%d_im" % (rd, curr_index))
    return hvib

