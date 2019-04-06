#*********************************************************************************
#* Copyright (C) 2019 Alexey V. Akimov
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


H = CMATRIX(2,2)
H.set(0,0, -0.001+0.0j);         H.set(0,1, 0.001+0.001j)
H.set(1,0, 0.001+0.001j);        H.set(1,1,  0.001+0.0j)
print "H = \n"; H.show_matrix()

U = CMATRIX(2,2)
S = CMATRIX(2,2)
V = CMATRIX(2,2)

print "============ Jacobi SVD ============="
JacobiSVD_decomposition(H, U, S, V)

print "U = \n"; U.show_matrix()
print "S = \n"; S.show_matrix()
print "V = \n"; V.show_matrix()
print "H = \n"; H.show_matrix()
print "U*S*V^* = \n"; (U*S*V.conj()).show_matrix()
print "U*S*V = \n"; (U*S*V).show_matrix()
print "U.H() * U = \n"; (U.H() * U).show_matrix()
print "V.H() * V = \n"; (V.H() * V).show_matrix()



print "============ BDC SVD ============="
BDCSVD_decomposition(H, U, S, V)

print "U = \n"; U.show_matrix()
print "S = \n"; S.show_matrix()
print "V = \n"; V.show_matrix()
print "H = \n"; H.show_matrix()
print "U*S*V^* = \n"; (U*S*V.conj()).show_matrix()
print "U*S*V = \n"; (U*S*V).show_matrix()
print "U.H() * U = \n"; (U.H() * U).show_matrix()
print "V.H() * V = \n"; (V.H() * V).show_matrix()



