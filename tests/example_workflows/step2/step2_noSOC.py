#***********************************************************
# * Copyright (C) 2017-2018 Brendan A. Smith and Alexey V. Akimov
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
sys.path.insert(1,"/projects/academic/alexeyak/brendan/official_workflows_libra/libra-code/workflows/pyxaid2")

import trajectory
import utils

"""
    this file is a reformulation of Pyxaid2 in terms of Libra's workflow module.
    Steps 2&3 in Pyxaid2 will each have their own "recipe", based on the functions
    defined in Libra's workflow module. Here, we use the functions defined in 
    Libra's workflow module to rcreate step2 in Pyxaid2.
"""

"""
    First, we need to initilize step2 - this means to properly create and distribute 
    the input files for the electronic structure calculations.
"""

### ========== From py-scr2.py ========== ###
nsteps_per_job = 2
tot_nsteps = 2

os.system("mkdir res")
#QE_methods.out2inp("x.md.out","x0.scf.in","wd","x0.scf",0,tot_nsteps,1)
utils.xyz2inp("CsPbI3.xyz","x0.scf.in","wd","x0.scf",0,tot_nsteps,1)

os.system("cp submit_templ.slm wd"); os.system("cp x0.exp.in wd"); os.chdir("wd")
QE_methods.distribute(0,tot_nsteps,nsteps_per_job,"submit_templ.slm",["x0.exp.in"],["x0.scf"],2)


