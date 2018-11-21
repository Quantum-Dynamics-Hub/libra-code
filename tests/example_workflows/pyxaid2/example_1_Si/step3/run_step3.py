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
import time
import math

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *
from libra_py.workflows.pyxaid2 import *

#############################################################################################
# This is an example of an input file to run namd calculations with SOC
# For test casesfor using NA-MD w/ slater determinants, we will only consider
# the diabatic cases for now. 
#############################################################################################

params = {}

rt="/mnt/c/cygwin/home/Alexey-user/Programming/Project_libra/tests/example_workflows/pyxaid2/example_1_Si/step2"
##### Extract Spin-diabatic Information #####
params["E_dia_ks_re_prefix"] = rt + "/res/E_dia_ks_"
params["E_dia_ks_re_suffix"] = "_re"
params["E_dia_ks_im_prefix"] = rt + "/res/E_dia_ks_"
params["E_dia_ks_im_suffix"] = "_im"

params["S_dia_ks_re_prefix"] = rt + "/res/S_dia_ks_"
params["S_dia_ks_re_suffix"] = "_re"
params["S_dia_ks_im_prefix"] = rt + "/res/S_dia_ks_"
params["S_dia_ks_im_suffix"] = "_im"

params["St_dia_ks_re_prefix"] = rt + "/res/St_dia_ks_"
params["St_dia_ks_re_suffix"] = "_re"
params["St_dia_ks_im_prefix"] = rt + "/res/St_dia_ks_"
params["St_dia_ks_im_suffix"] = "_im"

### Set up simulation specific parameters ###
# Initialize the SD basis
params["psi_dia_ks"] = range(1,7)  # 2 pairs of KS orbitals, indexing starts from 1

Phi_basis = []
Phi_0 = [ 3,-3 ]
Phi_1 = [ 3,-4 ]
Phi_2 = [ 4,-3 ]
Phi_3 = [ 3,-5 ]
Phi_4 = [ 5,-3 ]
Phi_5 = [ 3,-6 ]
Phi_6 = [ 6,-3 ]

Phi_basis.append( Phi_0 )
Phi_basis.append( Phi_1 ); Phi_basis.append( Phi_2 )
Phi_basis.append( Phi_3 ); Phi_basis.append( Phi_4 )
Phi_basis.append( Phi_5 ); Phi_basis.append( Phi_6 )
params["Phi_basis"] = Phi_basis

### Initilize the spin-adapated basis ###
coeff = []
coeff = [ [1,0,0,0,0,0,0], [0,1,-1,0,0,0,0], [0,0,0,1,-1,0,0], [0,0,0,0,0,1,-1] ]
P2C = CMATRIX(len(Phi_basis),len(coeff))
for i in xrange(len(Phi_basis)):
    for j in xrange(len(coeff)):
        P2C.set(i,j,coeff[j][i])

#P2C.show_matrix()
#sys.exit(0)

# Account for normalization constants #
N = []
for i in xrange(P2C.num_of_cols):
    count = 0.0
    N.append(0.0)
    for j in xrange(P2C.num_of_rows):
        if P2C.get(j,i) != 0:
            count += 1.0 
    N[i] = 1.0 / math.sqrt(count)
    for j in xrange(P2C.num_of_rows):
        P2C.set(j,i, N[i] * P2C.get(j,i))        



### ###
params["P2C"] = P2C
params["init_Chi"] = 3
params["Phi_dE"] = [];
for i in xrange(len(params["Phi_basis"])):
    params["Phi_dE"].append(0.0)
print params["Phi_dE"]
print params["Phi_basis"]

# Actual simulation paramters
params["init_time"] = 0  # starting from the first file 
params["len_traj"] = 30
params["do_state_reordering"] = 1
params["state_reordering_alpha"] = 0.01
params["do_phase_correction"] = 1
params["do_rescale"]  = 1 # boltzmann rescalling
params["sh_method"]   = 1 # 0 - MSSH, 1 - FSSH
params["do_collapse"] = 1 # 0 - no decoherence, 1 - decoherence (ID-A)
params["num_sh_traj"] = 1000
params["dt"] = 1.0
params["T"] = 300.0

print "\nPrinting params from run_step3.py"
print params

namd.run_namd(params)

