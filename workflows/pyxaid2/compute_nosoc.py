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
import math

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *

from utils import *
from trajectory import *


def compute_properties_gamma(params, es_curr, es_next, curr_index):

    Ca_curr = es_curr["Coeff_dia"][0]
    Ca_next = es_next["Coeff_dia"][0]
    Ea_curr = es_curr["E_dia"]
    Ea_next = es_next["E_dia"]

    rd = get_value(params,"rd",os.getcwd()+"../../res","s")   # of where the files will be printed out
    dt = get_value(params,"dt","1.0","f") # time step in fs - rescale NAC if actual dt is different
    dt = 41.34145 * dt  # convert to a.u., so the NACs are in a.u.
    hbar = 1.0 #0.658218 # units eV * fs
    is_st   = get_value(params,"compute_st",  "1","i")  # <current|next>
    is_nac  = get_value(params,"compute_nac", "1","i")  
    is_edia = get_value(params,"compute_edia","1","i")  
    is_hvib = get_value(params,"compute_hvib","1","i")  

    St, nac, edia, hvib = None, None, None, None
    #========== Computations ================
    if is_st==1 or is_nac==1 or is_hvib:    
        St = Ca_curr.H() * Ca_next

    if is_nac==1 or is_hvib==1:                                   
        nac = (0.5/dt)*(St - St.H())

    if is_edia==1 or is_hvib==1:
        edia = 0.5*(Ea_curr + Ea_next)

    if is_hvib==1:
        hvib = edia - 1.0j * hbar * nac


    S = Ca_curr.H() * Ca_curr
  
    os.system("mkdir %s" % rd)
    #========== Print out ================
    if is_st==1:    
        St.real().show_matrix("%s/St_dia_ks_%d_re" % (rd, curr_index))
        St.imag().show_matrix("%s/St_dia_ks_%d_im" % (rd, curr_index))

    if is_nac==1:
        nac.real().show_matrix("%s/nac_dia_ks_%d_re" % (rd, curr_index))
        nac.imag().show_matrix("%s/nac_dia_ks_%d_im" % (rd, curr_index))

    if is_edia==1:
        edia.real().show_matrix("%s/E_dia_ks_%d_re" % (rd, curr_index))
        edia.imag().show_matrix("%s/E_dia_ks_%d_im" % (rd, curr_index))

    if is_hvib==1:
        hvib.real().show_matrix("%s/hvib_dia_%d_re" % (rd, curr_index))
        hvib.imag().show_matrix("%s/hvib_dia_%d_im" % (rd, curr_index))


    return hvib




def compute_properties_general(params, es_curr, es_next, curr_index):

    print "You are dealing with multiple kpoints"

    as_sz = len(act_sp1)
    H = CMATRIX(info0["nk"]*as_sz, info0["nk"]*as_sz )
    S = CMATRIX(info0["nk"]*as_sz, info0["nk"]*as_sz )

    # optional orthogonalization - to mitigate the round off errors
    orthogonalize=0

    if orthogonalize==1:
        print "Do internal orbital orthogonalization"
        for ik1 in xrange(info0["nk"]):
            ovlp_cc = pw_overlap(info0["k"][ik1], info0["k"][ik1], coeff_curr0[ik1], coeff_curr0[ik1], grid_curr0[ik1], grid_curr0[ik1])
            ovlp_nn = pw_overlap(info0["k"][ik1], info0["k"][ik1], coeff_next0[ik1], coeff_next0[ik1], grid_next0[ik1], grid_next0[ik1])

            ovlp_cc_half = CMATRIX(ovlp_cc.num_of_rows, ovlp_cc.num_of_cols)
            ovlp_cc_i_half = CMATRIX(ovlp_cc.num_of_rows, ovlp_cc.num_of_cols)
            sqrt_matrix(ovlp_cc, ovlp_cc_half, ovlp_cc_i_half)
            coeff_curr0[ik1] = coeff_curr0[ik1] * ovlp_cc_i_half

            ovlp_nn_half = CMATRIX(ovlp_nn.num_of_rows, ovlp_nn.num_of_cols)
            ovlp_nn_i_half = CMATRIX(ovlp_nn.num_of_rows, ovlp_nn.num_of_cols)
            sqrt_matrix(ovlp_nn, ovlp_nn_half, ovlp_nn_i_half)
            coeff_next0[ik1] = coeff_next0[ik1] * ovlp_nn_i_half




    """
    The convention for the matrices for multiple k-points is:

          |  x_11  x_12 ... |
    X =   |  x_21  x_22 ... |
          |  ...      ...   |

    here, each x_ij block is a  as_sz x as_sz matrix describing
    the interactions of the as_sz orbitals of a k-point i and 
    as_sz orbitals of a k-point j

    So X is has a k-point-first block-structure

    """                     

    for ik1 in xrange(info0["nk"]):
        for ik2 in range(ik1, info0["nk"]):
            tim.start()
            ovlp_cc = pw_overlap(info0["k"][ik1], info0["k"][ik2], coeff_curr0[ik1], coeff_curr0[ik2], grid_curr0[ik1], grid_curr0[ik2])
            ovlp_nn = pw_overlap(info0["k"][ik1], info0["k"][ik2], coeff_next0[ik1], coeff_next0[ik2], grid_next0[ik1], grid_next0[ik2])
            #ovlp_nc = pw_overlap(info["k"][ik1], info["k"][ik2], coeff_next[ik1], coeff_curr[ik2], grid_next[ik1], grid_curr[ik2])
            ovlp_cn = pw_overlap(info0["k"][ik1], info0["k"][ik2], coeff_curr0[ik1], coeff_next0[ik2], grid_curr0[ik1], grid_next0[ik2])
    
            print "Time to compute 3 overlaps for the pair of k-points ", ik1, " ", ik2," is ", tim.stop()

            
            h_cc = CMATRIX(as_sz, as_sz)
            h_nn = CMATRIX(as_sz, as_sz)

            tim.start()
            for i1 in xrange(as_sz):
                for j1 in xrange(as_sz):
                    h_cc.set(i1, j1, 0.5*(e_curr0[ik1].get(i1,i1) + e_curr0[ik2].get(j1,j1))*ovlp_cc.get(i1,j1)) 
                    h_nn.set(i1, j1, 0.5*(e_next0[ik1].get(i1,i1) + e_next0[ik2].get(j1,j1))*ovlp_nn.get(i1,j1))

            h = 0.5*(h_cc + h_nn)  - (0.5j/dt)*(ovlp_cn - ovlp_cn.H()) 
            s = 0.5*(ovlp_cc + ovlp_nn)
    
            if ik2!=ik1:
                push_submatrix(S, s, range(ik1*as_sz, (ik1+1)*as_sz), range(ik2*as_sz, (ik2+1)*as_sz))
                push_submatrix(S, s.H(), range(ik2*as_sz, (ik2+1)*as_sz), range(ik1*as_sz, (ik1+1)*as_sz))

                push_submatrix(H, h, range(ik1*as_sz, (ik1+1)*as_sz), range(ik2*as_sz, (ik2+1)*as_sz))
                push_submatrix(H, h.H(), range(ik2*as_sz, (ik2+1)*as_sz), range(ik1*as_sz, (ik1+1)*as_sz))


            else:
                push_submatrix(S, s, range(ik1*as_sz, (ik1+1)*as_sz), range(ik2*as_sz, (ik2+1)*as_sz))
                push_submatrix(H, h, range(ik1*as_sz, (ik1+1)*as_sz), range(ik2*as_sz, (ik2+1)*as_sz))
            print "Time to push matrices is ", tim.stop()


