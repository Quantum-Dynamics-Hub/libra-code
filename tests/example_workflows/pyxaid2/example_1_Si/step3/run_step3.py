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
import namd

#############################################################################################
# This is an example of an input file to run namd calculations with SOC
# For test casesfor using NA-MD w/ slater determinants, we will only consider
# the diabatic cases for now. 
#############################################################################################

params = {}

rt="/projects/academic/alexeyak/brendan/libra-code/tests/example_workflows/pyxaid2/example_1_Si/step2"

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

##### Set up simulation specific parameters #####
# Spin-diabatic Parameters
params["psi_dia_ks"] = range(1,12)  # 2 pairs of KS orbitals, indexing starts from 1

Chi_0 = [ [ 1.0 , [1,-1,2,-2,3,-3,4,-4,5,-5,6,-6] ] ]
Chi_1 = [ [ 1.0 , [1,-1,2,-2,3,-3,4,-4,5,-5,6,-7] ] , [ -1.0 , [1,-1,2,-2,3,-3,4,-4,5,-5,-6,7] ] ]

params["Chi_basis"] = []
params["Chi_basis"].append( Chi_0 )
params["Chi_basis"].append( Chi_1 )

params["Phi_dE"] = [];
for i in xrange(len(params["psi_dia_ks"])):
    params["Phi_dE"].append(0.0)

print params["Phi_dE"]
print params["Chi_basis"]

params["init_Chi"] = 1  # Spin-Adapated Singlet Wavefunction 

# Actual simulation paramters
params["init_time"] = 0  # starting from the first file 
params["len_traj"] = 1000
params["sh_method"] = 1   # 0, 1
params["do_collapse"] = 0 # 0 - no decoherence, 1 - decoherence
params["num_sh_traj"] = 100
params["dt"] = 1
params["T"] = 400

print "\nPrinting params from run_step3.py"
print params

namd.run_namd(params)

