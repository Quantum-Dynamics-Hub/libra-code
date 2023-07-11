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


print "============ Matrices =============="
H = CMATRIX(2,2)
H.set(0,0, -0.001+0.0j);         H.set(0,1, 0.001+0.001j)
H.set(1,0, 0.001+0.001j);        H.set(1,1,  0.002+0.0j)
print "H = \n"; H.show_matrix()

I = CMATRIX(2,2)
I.set(0,0, 1.0, 0.0);   I.set(0,1, 0.0, 0.0)
I.set(1,0, 0.0, 0.0);   I.set(1,1, 1.0, 0.0)
print "I = \n"; I.show_matrix()

H2 = CMATRIX(2,2)
H2 = H.H() * H
H2p = H * H.H()
print "H2 = \n"; H2.show_matrix()
print "H2p = \n"; H2p.show_matrix()


print "============ Eigenvectors =============="
E2 = CMATRIX(2,2)
C2 = CMATRIX(2,2)

E2p = CMATRIX(2,2)
C2p = CMATRIX(2,2)

"""
   H2 * C2 = I * C2 * E2

   H2p * C2p = I * C2p * E2p

"""

solve_eigen(H2, I, E2, C2, 0)
solve_eigen(H2p, I, E2p, C2p, 0)

print "E2 = \n"; E2.show_matrix()
print "C2 = \n"; C2.show_matrix()

print "E2p = \n"; E2p.show_matrix()
print "C2p = \n"; C2p.show_matrix()


"""
Assume that:  H  = U * S * V.H()
  
             H.H() * H = V * S.H() * U.H() * U * S * V.H()

             H.H() * H = V * S.H() * S * V.H()

             H.H() * H * V = I * V * S.H() * S

That is: V = C2 
         E2 = S.H() * S

Likewise:

           H * H.H() = U * S * V.H() * V * S.H() * U.H()

           H * H.H() * U = U * S * S.H()

"""


U = CMATRIX(2,2)
S = CMATRIX(2,2)
V = CMATRIX(2,2)

print "============ Jacobi SVD ============="
JacobiSVD_decomposition(H, U, S, V)


print "S = \n"; S.show_matrix()
print "U*S*V.H() (should be same as H) = \n"; (U*S*(V.H())).show_matrix()
print "U.H() * U = \n"; (U.H() * U).show_matrix()
print "V.H() * V = \n"; (V.H() * V).show_matrix()
print "V (should be same as C2) = \n"; V.show_matrix()
print "S.H() * S (should be same as E2) = \n"; (S.H()*S).show_matrix()
#print "U (should be same as C2p) = \n"; U.show_matrix()
#print "S * S.H() (should be same as E2p) = \n"; (S*S.H()).show_matrix()




print "============ BDC SVD ============="
BDCSVD_decomposition(H, U, S, V)


print "S = \n"; S.show_matrix()
print "U*S*V.H() (should be same as H) = \n"; (U*S*(V.H())).show_matrix()
print "U.H() * U = \n"; (U.H() * U).show_matrix()
print "V.H() * V = \n"; (V.H() * V).show_matrix()
print "V (should be same as C2) = \n"; V.show_matrix()
print "S.H() * S (should be same as E2) = \n"; (S.H()*S).show_matrix()
#print "U (should be same as C2p) = \n"; U.show_matrix()
#print "S * S.H() (should be same as E2p) = \n"; (S*S.H()).show_matrix()


