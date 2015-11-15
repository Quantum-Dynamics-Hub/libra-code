#*********************************************************************************
#* Copyright (C) 2015 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/

import os
import sys
import math

# Fisrt, we add the location of the library to test to the PYTHON path
cwd = os.getcwd()
print "Current working directory", cwd
sys.path.insert(1,cwd+"/../../_build/src/mmath")
sys.path.insert(1,cwd+"/../../_build/src/mmath/meigen")
sys.path.insert(1,cwd+"/../../_build/src/mmath/linalg")


print "\nTest 1: Importing the library and its content"
from cygmmath import *
from cygmeigen import *
from cyglinalg import *


print "\n Test2: Setting up "
H = MATRIX(2,2)
H.set(0,0, -0.001);  H.set(0,1, 0.001)
H.set(1,0, 0.001);   H.set(1,1,  0.001)


S = MATRIX(2,2)
S.set(0,0, 1.0);   S.set(0,1, 0.0)
S.set(1,0, 0.0);   S.set(1,1, 1.0)

E = CMATRIX(2,2)
C = CMATRIX(2,2)

solve_eigen(2, H, S, E, C)

C.show()


print "H = \n"; H.show_matrix()

X = CMATRIX(H)
print "X = \n"; X.show()

