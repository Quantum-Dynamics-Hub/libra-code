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
from workflows.pyxaid2 import *


# Example of usage
opt = 1
scl1 = 2.0*13.60569253 # Ha to eV
scl2 = 1000.0 # some scaling to get plottable data
tmin = 0
tmax = 100

HOMO = 1
minE = 0.0
maxE = 10.0
dE = 0.05

# Energy, couplings and H' in space of orbital indexes
excitation_spectrum.ham_map("res/0_Ham_",   tmin,tmax,"_re" ,opt,scl1,"spectr/ave_Ham_re.dat")
excitation_spectrum.ham_map("res/0_Ham_",   tmin,tmax,"_im" ,opt,scl1,"spectr/ave_Ham_im.dat")

