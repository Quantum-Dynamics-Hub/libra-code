#*********************************************************************************
#* Copyright (C) 2015-2016 Alexey V. Akimov
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
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *



print "\n Test2: Setting up "
H = MATRIX(2,2)
H.set(0,0, -0.001);  H.set(0,1, 0.001)
H.set(1,0, 0.001);   H.set(1,1,  0.001)
print "H = \n"; H.show_matrix()

cH = CMATRIX(H); cH.show_matrix()

S = MATRIX(2,2)
S.set(0,0, 1.0);   S.set(0,1, 0.5)
S.set(1,0, 0.5);   S.set(1,1, 1.0)
print "S = \n"; S.show_matrix()

cS = CMATRIX(S); cS.show_matrix()


E = CMATRIX(2,2)
C = CMATRIX(2,2)

solve_eigen(2, H, S, E, C)

print "E = \n"; E.show_matrix()
print "C = \n"; C.show_matrix()

print "H*C = \n"; (cH*C).show_matrix()
print "S*C*E = \n"; (cS*C*E).show_matrix()

print "C.H() * S * C\n"; (C.H() * cS * C).show_matrix()

#X = CMATRIX(H)
#print "X = \n"; X.show_matrix()

