#***********************************************************
# * Copyright (C) 2017-2018 Brendan A. Smith, Wei Li, and Alexey V. Akimov
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

from utils import * 


def compute_properties_gamma(params, es_curr, es_next, curr_index):

    Ca_curr = es_curr["Coeff_adi"][0]
    Ca_next = es_next["Coeff_adi"][0]
    Ea_curr = es_curr["E_adi"]
    Ea_next = es_next["E_adi"]       

    rd = get_value(params,"rd",os.getcwd()+"../../res","s")   # of where the files will be printed out
    dt = get_value(params,"dt","1.0","f") # time step in fs - rescale NAC if actual dt is different
    dt = 41.34145 * dt  # convert to a.u., so the NACs are in a.u.

    is_st   = get_value(params,"compute_st",  "1","i")  # <current|next>
    is_nac  = get_value(params,"compute_nac", "1","i")  
    is_eadi = get_value(params,"compute_eadi","1","i")  
    is_hvib = get_value(params,"compute_hvib","1","i")  

    St, nac, eadi, hvib = None, None, None, None

    #========== Computations ================
    if is_st==1 or is_nac==1 or is_hvib:    
        St = Ca_curr.H() * Ca_next

    if is_nac==1 or is_hvib==1:
        nac = (0.5/dt)*(St - St.H())

    if is_eadi==1 or is_hvib==1:
        eadi = 0.5*(Ea_curr + Ea_next)

    if is_hvib==1:
        hvib = eadi - 1.0j * nac

    S = Ca_curr.H() * Ca_curr
    S_adi_pw = Ca_curr

    ##### Pauli matrices ###
    #
    #        | 0  1 |         | 0  -i |         | 1   0 |
    #  sig1 =|      |  sig2 = |       |  sig3 = |       | 
    #        | 1  0 |         | i   0 |         | 0  -1 |
    #
    ######

    sx = S.num_of_cols
    sig1 = CMATRIX(sx/2, sx/2)
    sig2 = CMATRIX(sx/2, sx/2)
    sig3 = CMATRIX(sx/2, sx/2)
    for n in xrange(sx/2):
        for k in xrange(sx/2):
            sig1.set(n,k, S.get(2*n,2*k+1) + S.get(2*n+1,2*k) )
            sig2.set(n,k, (-S.get(2*n,2*k+1) + S.get(2*n+1,2*k))*(1.0j+0.0) )
            sig3.set(n,k, S.get(2*n,2*k) - S.get(2*n+1,2*k+1) )

    #sig1.real().show_matrix("%s/0_sig1_%d_re" % (rd, curr_index) )
    #sig1.imag().show_matrix("%s/0_sig1_%d_im" % (rd, curr_index) )
    #sig2.real().show_matrix("%s/0_sig2_%d_re" % (rd, curr_index) )
    #sig2.imag().show_matrix("%s/0_sig2_%d_im" % (rd, curr_index) )
    #sig3.real().show_matrix("%s/0_sig3_%d_re" % (rd, curr_index) )
    #sig3.imag().show_matrix("%s/0_sig3_%d_im" % (rd, curr_index) )
    S.real().show_matrix("%s/S_adi_ks_%d_re" % (rd, curr_index) )
    S.imag().show_matrix("%s/S_adi_ks_%d_im" % (rd, curr_index) )
    S_adi_pw.real().show_matrix("%s/S_adi_pw_%d_re" % (rd, curr_index) )
    S_adi_pw.imag().show_matrix("%s/S_adi_pw_%d_im" % (rd, curr_index) )

  
    os.system("mkdir %s" % rd)
    #========== Print out ================
    if is_st==1:    
        St.real().show_matrix("%s/St_adi_ks_%d_re" % (rd, curr_index))
        St.imag().show_matrix("%s/St_adi_ks_%d_im" % (rd, curr_index))

    if is_nac==1:
        nac.real().show_matrix("%s/nac_adi_ks_%d_re" % (rd, curr_index))
        nac.imag().show_matrix("%s/nac_adi_ks_%d_im" % (rd, curr_index))

    if is_eadi==1:
        eadi.real().show_matrix("%s/E_adi_ks_%d_re" % (rd, curr_index))
        eadi.imag().show_matrix("%s/E_adi_ks_%d_im" % (rd, curr_index))

    if is_hvib==1:
        hvib.real().show_matrix("%s/hvib_adi_%d_re" % (rd, curr_index))
        hvib.imag().show_matrix("%s/hvib_adi_%d_im" % (rd, curr_index))

    return hvib


