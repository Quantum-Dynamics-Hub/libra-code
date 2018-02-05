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
# Tutorial: Create interactions and compute them
###################################################################

import os
import sys
import math

cwd = os.getcwd()
print "Current working directory", cwd
sys.path.insert(1,cwd+"/../../_build/src/hamiltonian/Hamiltonian_Atomistic/Hamiltonian_MM")
sys.path.insert(1,cwd+"/../../_build/src/converters")
sys.path.insert(1,cwd+"/../../_build/src/math_linalg")

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    #from cyglibra_core import *
    from cygconverters import *
    from cyghamiltonian_mm import *
    from cyglinalg import *

elif sys.platform=="linux" or sys.platform=="linux2":
    #from liblibra_core import *
    from libconverters import *
    from libhamiltonian_mm import *
    from liblinalg import *


Rall = Py2Cpp_VECTOR(
       [VECTOR(0.0, 0.0, 0.0), 
        VECTOR(1.0, 0.0, 0.0),
        VECTOR(2.0, 0.0, 0.0),
        VECTOR(3.0, 0.0, 0.0),
        VECTOR(4.0, 0.0, 0.0),
        VECTOR(5.0, 0.0, 0.0)
       ] )

Tall = Py2Cpp_VECTOR(
       [VECTOR(0.0, 0.0, 0.0), 
        VECTOR(-4.0, 0.0, 0.0),
        VECTOR(0.0, 0.0, 0.0),
        VECTOR(0.0, 0.0, 0.0),
        VECTOR(0.0, 0.0, 0.0),
        VECTOR(0.0, 0.0, 0.0)
       ] )


Fall = Py2Cpp_VECTOR(
       [VECTOR(0.0, 0.0, 0.0), 
        VECTOR(0.0, 0.0, 0.0),
        VECTOR(0.0, 0.0, 0.0),
        VECTOR(0.0, 0.0, 0.0),
        VECTOR(0.0, 0.0, 0.0),
        VECTOR(0.0, 0.0, 0.0)
       ] )

#i1 = Interaction_N_Body(2)
#i1 = Interaction_2_Body()
i1 = Bond_Interaction()
r1 = VECTOR(0.0, 0.0, 0.0)
r2 = VECTOR(1.0, 0.0, 0.0)
#i1.set_coords(r1,r2)
#i1.set_coords(Rall[0], Rall[1])
i1.set_coords(Rall, Py2Cpp_int([0,1]))
i1.set_transl(Tall, Py2Cpp_int([0,1]))
i1.set_forces(Fall, Py2Cpp_int([0,1]))

i1.set_functional("Harmonic")
i1.set_params({"K":10.0, "r0":2.0})
i1.compute()
i1.compute()

print i1.int_type, i1.functional, i1.K, i1.r0, i1.energy
print Fall[0].x, Fall[0].y, Fall[0].z
print Fall[1].x, Fall[1].y, Fall[1].z







