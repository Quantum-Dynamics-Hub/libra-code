#*********************************************************************************
#* Copyright (C) 2017 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
###################################################################
#
# This example demonstrates how to build a quantum dot (QD) 
# of a given radius
# 
###################################################################


import os
import sys
import math
import copy

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *


## Si QD

a = [     5.4307100000000000,    0.0000000000000000,    0.0000000000000000 ]
b = [     0.0000000000000000,    5.4307100000000000,    0.0000000000000000 ]
c = [     0.0000000000000000,    0.0000000000000000,    5.4307100000000000 ]
build.generate_replicas_xyz(a, b, c, 5, 5, 5 , "Si.xyz", "si-5.xyz")
build.crop_sphere_xyz("si-5.xyz", "si_5.xyz", 5.0)


## CdSe QD
a = [5.833, 0.0, 0.0]
b = [0.0, 5.833, 0.0]
c = [0.0, 0.0, 5.833]
build.generate_replicas_xyz(a, b, c, 20, 20, 20 , "CdSe.xyz", "cdse-20.xyz")
label, R = build.crop_sphere_xyz("cdse-20.xyz", "cdse_15.xyz", 15.0)
g = [VECTOR(0.0, 0.0, 0.0)] * len(R)


## And this is how you can create a System object containing this QD
rnd = Random()
sz = len(R)

# Note, it is typical that .xyz file would contain the coordinates in Angstrom units
# so this is what numbers in R are. This means we need to convert these data to Bohrs
#
for i in xrange(sz):
    R[i] = R[i] * 1.889725989  # convert Angstrom -> Bohr

syst = init_system.init_system(label, R, g, rnd, 300.0, 0.01, 0, "elements.dat")
syst.print_ent("cdse_15.ent")

