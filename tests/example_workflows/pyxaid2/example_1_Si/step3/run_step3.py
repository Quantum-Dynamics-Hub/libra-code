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

import libra_py.workflows.nbra.step3 as step3

#############################################################################################
# This is an example of an input file to run namd calculations with SOC
# For test casesfor using NA-MD w/ slater determinants, we will only consider
# the diabatic cases for now. 
#############################################################################################

setup = 1  # both setups should give the same results


params = {}

params["nsteps"] = 30
params["dt"] = 41.3413  # in a.u.
params["is_pyxaid_format"] = False
params["do_state_reordering"] = 2
params["state_reordering_alpha"] = 0.0
params["do_phase_correction"] = 1

params["data_set_paths"] = ["/mnt/c/cygwin/home/Alexey-user/Programming/Project_libra/tests/example_workflows/pyxaid2/example_1_Si/step2/res/"]
params["E_re_prefix"] = "E_dia_ks_"; params["E_re_suffix"] = "_re"
params["E_im_prefix"] = "E_dia_ks_"; params["E_im_suffix"] = "_im"
params["S_re_prefix"] = "S_dia_ks_"; params["S_re_suffix"] = "_re"
params["S_im_prefix"] = "S_dia_ks_"; params["S_im_suffix"] = "_im"
params["St_re_prefix"] = "St_dia_ks_"; params["St_re_suffix"] = "_re"
params["St_im_prefix"] = "St_dia_ks_"; params["St_im_suffix"] = "_im"

if setup==1:
    params["output_set_paths"] = ["/mnt/c/cygwin/home/Alexey-user/Programming/Project_libra/tests/example_workflows/pyxaid2/example_1_Si/step3/res_setup1/"]
    os.system("mkdir res_setup1")
elif setup==2:
    params["output_set_paths"] = ["/mnt/c/cygwin/home/Alexey-user/Programming/Project_libra/tests/example_workflows/pyxaid2/example_1_Si/step3/res_setup2/"]
    os.system("mkdir res_setup2")

params["Hvib_re_prefix"] = "Hvib_"; params["Hvib_re_suffix"] = "_re"
params["Hvib_im_prefix"] = "Hvib_"; params["Hvib_im_suffix"] = "_im"


### Set up basis and basis transformations ###
"""
    Example: 

    Setup 1

    What is in file   -->   What we actually need             --->  How states are indexed in the SDs

    norbitals = 12      active_space = range(2,6)+range(5,12)     numbers used in params["Phi_basis"]

    alpha  beta
            11 ----------------------------> 11 --------------------------> -4
            10 ----------------------------> 10 --------------------------> -3
         (L) 9 ---------------------------->  9 --------------------------> -2
         (H) 8 ---------------------------->  8 --------------------------> -1
             7 
             6
     5 ---------------------------->  5  -------------------------->    4
     4 ---------------------------->  4  -------------------------->    3
 (L) 3 ---------------------------->  3  -------------------------->    2
 (H) 2 ---------------------------->  2  -------------------------->    1
     1
     0



    Setup 2

    What is in file   -->   What we actually need             --->  How states are indexed in the SDs

    norbitals = 12      active_space = range(0,12)     numbers used in params["Phi_basis"]

    alpha  beta
            11 ----------------------------> 11 --------------------------> -6
            10 ----------------------------> 10 --------------------------> -5
         (L) 9 ---------------------------->  9 --------------------------> -4
         (H) 8 ---------------------------->  8 --------------------------> -3
             7 ---------------------------->  7 --------------------------> -2
             6 ---------------------------->  6 --------------------------> -1
     5 ---------------------------->  5  -------------------------->    6
     4 ---------------------------->  4  -------------------------->    5
 (L) 3 ---------------------------->  3  -------------------------->    4
 (H) 2 ---------------------------->  2  -------------------------->    3
     1 ---------------------------->  1  -------------------------->    2
     0 ---------------------------->  0  -------------------------->    1


"""


if setup==1:
    # ===  Basis of KS orbitals ====
    params["norbitals"] = 12  # 6 pairs of KS orbitals, 12 spin-orbitals
    params["active_space"] = range(2,6)+range(8,12)

    # ===  Basis of SDs (Phi) ====
    params["Phi_basis"] = [ [ 1,-1 ], [ 1,-2 ], [ 2,-1 ], [ 1,-3 ], [ 3,-1 ], [ 1,-4 ], [ 4,-1 ] ]

elif setup==2:
    # ===  Basis of KS orbitals ====
    params["norbitals"] = 12  # 6 pairs of KS orbitals, 12 spin-orbitals
    params["active_space"] = range(0,12)

    # ===  Basis of SDs (Phi) ====
    params["Phi_basis"] = [ [ 3,-3 ], [ 3,-4 ], [ 4,-3 ], [ 3,-5 ], [ 5,-3 ], [ 3,-6 ], [ 6,-3 ] ]

# ===  Basis of spin-adapted SDs (Chi) ====
params["P2C"] = [ [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ], 
                  [0.0, 1.0,-1.0, 0.0, 0.0, 0.0, 0.0 ],
                  [0.0, 0.0, 0.0, 1.0,-1.0, 0.0, 0.0 ],
                  [0.0, 0.0, 0.0, 0.0, 0.0, 1.0,-1.0 ] 
                ]

params["Phi_dE"] = [];
for i in xrange(len(params["Phi_basis"])):
    params["Phi_dE"].append(0.0)


step3.run(params)



"""
I keep it here just to have it handy - need to find out what to do with this

### Initialize parameters to scale energies and NACs ###
en_gap = 3.21*units.ev2Ha # experimental band gap, units in Ha
params["Chi_en_gap"] = en_gap
# Scale NAC < i | d/dt | j >. Format = [ [i,j], scaling factor ]
params["do_scale"] = 1
params["sc_nac_method"] = 0
params["NAC_dE"] = []
params["NAC_dE"].append( [ [0,1], 1.0 ] )
print params["Phi_dE"]
print params["NAC_dE"]
"""

