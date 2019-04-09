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

print "============ Real Matrix =============="
H = MATRIX(2,2)
H.set(0,0, 1.0);        H.set(0,1, 2.0)
H.set(1,0, 2.0);        H.set(1,1, 4.0)
print "H = \n"; H.show_matrix()

dt = 1.0

print "Method 1"
eH1 = exp_(H, dt)
eH1.show_matrix()

print "Method 1a"
eH1 = exp_(H, 0.5*dt)
eH1 = eH1 * eH1
eH1.show_matrix()

print "Method 2"
eH2 = exp_2(H, dt)
eH2.show_matrix()

print "Method 2a"
eH2 = exp_2(H, 0.5*dt)
eH2 = eH2 * eH2
eH2.show_matrix()

#print "Method 2, only 2 terms"
#eH2 = exp_2(H, dt, 10000, 1e-10)
#eH2.show_matrix()



print "============ Hermitian Matrix =============="
H = CMATRIX(2,2)
H.set(0,0, -0.001+0.0j);         H.set(0,1, 0.001-0.001j)
H.set(1,0, 0.001+0.001j);        H.set(1,1,  0.002+0.0j)
print "H = \n"; H.show_matrix()

dt = 1.0+0.0j


print "Method 1"
eH1 = exp_(H, dt)
eH1.show_matrix()

print "Method 1a"
eH1 = exp_(H, 0.5*dt)
eH1 = eH1 * eH1
eH1.show_matrix()

print "Method 2"
eH2 = exp_2(H, dt)
eH2.show_matrix()

print "Method 2a"
eH2 = exp_2(H, 0.5*dt)
eH2 = eH2 * eH2
eH2.show_matrix()

#print "Method 2, only 2 terms"
#eH2 = exp_2(H, dt, 10000, 1e-10)
#eH2.show_matrix()



print "============ non-Hermitian Matrix =============="
H = CMATRIX(2,2)
H.set(0,0, -0.001+0.0j);         H.set(0,1, 0.001+0.001j)
H.set(1,0, 0.001+0.001j);        H.set(1,1,  0.002+0.0j)
print "H = \n"; H.show_matrix()

dt = 0.0-1.0j

print "Methods 1 and 1a do not make sense for non-Hermitian matrices"
print "Method 1"
eH1 = exp_(H, dt)
eH1.show_matrix()

print "Method 1a"
eH1 = exp_(H, 0.5*dt)
eH1 = eH1 * eH1
eH1.show_matrix()

print "Method 2"
eH2 = exp_2(H, dt)
eH2.show_matrix()

print "Method 2a"
eH2 = exp_2(H, 0.5*dt)
eH2 = eH2 * eH2
eH2.show_matrix()

#print "Method 2, only 2 terms"
#eH2 = exp_2(H, dt, 10000, 1e-10)
#eH2.show_matrix()


