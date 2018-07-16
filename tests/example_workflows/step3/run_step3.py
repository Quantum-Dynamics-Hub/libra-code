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

#cwd = os.getcwd()
#print "Current working directory", cwd

# Set path to "step3" folder, and please make sure run_namd.py is in this folder.
sys.path.insert(1,"/projects/academic/alexeyak/brendan/libra/libra-code/tests/pyxaid2_libra_test/step3")
import run_namd

#############################################################################################
# This is an example of an input file to run namd calculations with SOC
# For test casesfor using NA-MD w/ slater determinants, we will only consider
# the diabatic cases for now. 
#############################################################################################

params = {}

rt="/gpfs/scratch/brendan/perov/CsPbI3/pyxaid"

##### Extract diabatic Information #####
#"""
params["E_dia_ks_re_prefix"] = rt + "/res/E_dia_ks_"
params["E_dia_ks_re_suffix"] = "_re"
params["E_dia_ks_im_prefix"] = rt + "/res/E_dia_ks_"
params["E_dia_ks_im_suffix"] = "_im"

params["S_dia_ks_re_prefix"] = rt + "/res/0_S_dia_ks_"
params["S_dia_ks_re_suffix"] = "_re"
params["S_dia_ks_im_prefix"] = rt + "/res/0_S_dia_ks_"
params["S_dia_ks_im_suffix"] = "_im"

params["St_dia_ks_re_prefix"] = rt + "/res/St_dia_ks_"
params["St_dia_ks_re_suffix"] = "_re"
params["St_dia_ks_im_prefix"] = rt + "/res/St_aia_ks_"
params["St_dia_ks_im_suffix"] = "_im"
#"""

##### Extract adibatic Information #####
"""
params["H_vib_re_prefix"] = rt + "/res/0_Ham_"
params["H_vib_re_suffix"] = "_re"
params["H_vib_im_prefix"] = rt + "/res/0_Ham_"
params["H_vib_im_suffix"] = "_im"

params["E_adi_ks_re_prefix"] = rt + "/res/E_adi_ks_"
params["E_adi_ks_re_suffix"] = "_re"
params["E_adi_ks_im_prefix"] = rt + "/res/E_adi_ks_"
params["E_adi_ks_im_suffix"] = "_im"

params["S_adi_ks_re_prefix"] = rt + "/res/0_S_adi_ks_"
params["S_adi_ks_re_suffix"] = "_re"                    
params["S_adi_ks_im_prefix"] = rt + "/res/0_S_adi_ks_"
params["S_adi_ks_im_suffix"] = "_im"                    

params["St_adi_ks_re_prefix"] = rt + "/res/St_adi_ks_"
params["St_adi_ks_re_suffix"] = "_re"
params["St_adi_ks_im_prefix"] = rt + "/res/St_adi_ks_"
params["St_adi_ks_im_suffix"] = "_im"
"""

##### Extract diabatic, adiabatic pw for dia2adi projections #####

# NOTE - This section is under construction, but is not needed for diabatic
# dynamics. 

"""
params["S_dia_pw_re_prefix"] = rt+"/res/S_dia_pw_"  
params["S_dia_pw_re_suffix"] = "_re"
params["S_dia_pw_im_prefix"] = rt+"/res/S_dia_pw_"
params["S_dia_pw_im_suffix"] = "_im"

params["S_adi_pw_re_prefix"] = rt+"/res/S_adi_pw_"          
params["S_adi_pw_re_suffix"] = "_re"
params["S_adi_pw_im_prefix"] = rt+"/res/S_adi_pw_"
params["S_adi_pw_im_suffix"] = "_im"
"""


##### Set up simulation specific parameters #####
## Adiabatic Parameters
#params["psi_adi_ks"] = [1,2,3,4]  #2-component spinor KS orbitals
#params["Psi_basis"] = [];                    params["Psi_dE"] = [];
#params["Psi_basis"].append([1,-1,2,-2]);     params["Psi_dE"].append(0.0);
#params["Psi_basis"].append([1,-1,2,-3]);     params["Psi_dE"].append(0.0);
#params["Psi_basis"].append([1,2,3,4,9,10]);  params["Psi_dE"].append(0.0);

## Diabatic Parameters
params["psi_dia_ks"] = range(1,22)  # 2 pairs of KS orbitals, indexing starts from 1
params["Chi_basis"] = [];                                               params["Phi_dE"] = [];
params["Chi_basis"].append( [ [1.0, [11,-11]] ] );                     params["Phi_dE"].append(0.0);
params["Chi_basis"].append( [ [1.0, [15,-15]] ] );                     params["Phi_dE"].append(0.0);
#params["Chi_basis"].append( [ [1.0, [ 10, 11]] ] );                     params["Phi_dE"].append(0.0);
#params["Chi_basis"].append( [ [1.0, [ 11,-11]] ] );                     params["Phi_dE"].append(0.0);
#params["Chi_basis"].append( [ [1.0, [ 10,-11]], [-1.0, [-10, 11]] ] );  params["Phi_dE"].append(0.0);
params["init_Chi"] = 1  # Spin-Adapated Singlet Wavefunction 

# Actual simulation paramters
params["init_time"] = 0  # starting from the first file 
params["len_traj"] = 75
params["sh_method"] = 2 # 0, 1, 2 # 2 = boltzman scalling
params["do_collapse"] = 0 # 0 - no decoherence, 1 - decoherence
params["num_sh_traj"] = 100
params["dt"] = 1
params["T"] = 300

print "\nPrinting params from run_step3.py"
print params

run_namd.run_namd(params)

