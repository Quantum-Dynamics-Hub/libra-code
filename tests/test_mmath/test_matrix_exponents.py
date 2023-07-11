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
import cmath

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
H.set(0,0, -0.001+0.0j);         H.set(0,1, 0.002+0.003j)
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

print "Method 2 and 2a in principle work even for non-Hermitian matrices"
eH2 = exp_2(H, dt)
eH2.show_matrix()

print "Method 2a"
eH2 = exp_2(H, 0.5*dt)
eH2 = eH2 * eH2
eH2.show_matrix()


print "Method 3 is based on eigendecomposition"
R = CMATRIX(2,2);  R_dag = CMATRIX(2,2);  iR = CMATRIX(2,2)
L = CMATRIX(2,2);  L_dag = CMATRIX(2,2);  iL_dag = CMATRIX(2,2)
ER = CMATRIX(2,2); EL_dag = CMATRIX(2,2)

# H * R = R * E
solve_eigen_nosort(H, ER, R, 0)
R_dag = R.H()

H_dag = H.H()
solve_eigen_nosort(H_dag, EL_dag, L, 0)
L_dag = L.H()

FullPivLU_inverse(R, iR)
FullPivLU_inverse(L_dag, iL_dag)

print "ER = "; ER.show_matrix()
print "EL_dag = "; EL_dag.show_matrix()

D = L.H() * R;

print "L.H() * R = "; D.show_matrix()
print "|D(0,0)| = ", abs(D.get(0,0))
print "|D(1,1)| = ", abs(D.get(1,1))


eE3 = CMATRIX(2,2)
for i in xrange(2):
    eE3.set(i,i,  cmath.exp(dt*ER.get(i,i)) )

print "exp(dt*ER) = "; eE3.show_matrix()

exH3 = iL_dag * D * eE3 * iR

print "exp(dt*H) = "; exH3.show_matrix()

#sys.exit(0)


print "============ non-Hermitian Matrix =============="
H = CMATRIX(3,3)
H.set(0,0, -0.001+0.0j);         H.set(0,1, 0.002+0.003j);   H.set(0,2, 0.002+0.003j)
H.set(1,0, 0.001+0.001j);        H.set(1,1,  0.002+0.0j);    H.set(1,2,  -0.002+0.1j)
H.set(2,0, 0.001-0.001j);        H.set(2,1,  0.002+0.1j);    H.set(2,2,   0.004+0.0j)
print "H = \n"; H.show_matrix()

R = CMATRIX(3,3);  R_dag = CMATRIX(3,3); iR = CMATRIX(3,3)
L = CMATRIX(3,3);  L_dag = CMATRIX(3,3); iL_dag = CMATRIX(3,3)
ER = CMATRIX(3,3); EL_dag = CMATRIX(3,3)

# H * R = R * E
solve_eigen_nosort(H, ER, R, 0)
R_dag = R.H()

H_dag = H.H()
solve_eigen_nosort(H_dag, EL_dag, L, 0)
L_dag = L.H()

FullPivLU_inverse(R, iR)
FullPivLU_inverse(L_dag, iL_dag)

print "ER = "; ER.show_matrix()
print "EL_dag = "; EL_dag.show_matrix()
D = L.H() * R
print "L.H() * R = "; D.show_matrix()
print "|D(0,0)| = ", abs(D.get(0,0))
print "|D(1,1)| = ", abs(D.get(1,1))
print "|D(2,2)| = ", abs(D.get(2,2))

eE3 = CMATRIX(3,3)
for i in xrange(3):
    eE3.set(i,i,  cmath.exp(dt*ER.get(i,i)) )

print "exp(dt*ER) = "; eE3.show_matrix()

exH3 = iL_dag * D * eE3 * iR
print "exp(dt*H) = "; exH3.show_matrix()


print "Method 2:"
eH2 = exp_2(H, dt)
eH2.show_matrix()

