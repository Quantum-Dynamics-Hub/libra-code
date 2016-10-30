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



print "\n Test2: Setting up complex matrix also test show() function"
print "a = CMATRIX(4,4)"
print "for i in range(0,4):"
print "    for j in range(0,4):"
print "        a.set(i,j,0.5*i+j, 0.1*(i-j))"
a = CMATRIX(4,4)
for i in range(0,4):
    for j in range(0,4):
        a.set(i,j,0.5*i+j, 0.1*(i-j))

print "for i in range(0,4):"
print "    for j in range(0,4):"
print "        print i,j, a.get(i,j)"
for i in range(0,4):
    for j in range(0,4):
        print i,j, a.get(i,j)

print "a.show()"
a.show()


print "\n Test3: Basic operations on a single matrix"
print "Transpose..."
print "x = a.T()"
x = a.T()
print "print a"
print a
print "print x"
print x
print "x.show()"
x.show()

print "Conjugate..."
print "x = a.conj()"
x = a.conj()
print "print a"
print a
print "print x"
print x
print "x.show()"
x.show()

print "Hermitian..."
print "x = a.H()"
x = a.H()
print "print a"
print a
print "print x"
print x
print "x.show()"
x.show()


print "load_identity..."
print "x.load_identity()"
x.load_identity()
print "print x"
print x
print "x.show()"
x.show()

print "t1 = x.col(0)"
t1 = x.col(0)
print "print t1"
print t1
print "t1.show()"
t1.show()

print "t1 = x.row(3)"
t2 = x.row(3)
print "print t2"
print t2
print "t2.show()"
t2.show()


print "\n Test4: Matrix multiplications"
print "(t1*t2).show()"
(t1*t2).show()

print "(t2*t1).show()"
(t2*t1).show()


print "\n Test5: Fun with symmetric (Hermitian) matrices"
print "Here we make a Hermitian matrix"
print "H = a + a.H()"
print "H.show()"
H = a + a.H()
H.show()


print "Tridiagonalization..."
print "T = CMATRIX(4,4)"
T = CMATRIX(4,4)
print "H.tridiagonalize(T)"
H.tridiagonalize(T)
print "T.show()"
T.show()


print "Tridiagonalization with Householder matrices..."
print "T = CMATRIX(4,4)"
hoho = CMATRIX(4,4)  # Householder matrix
print "H.tridiagonalize(T,hoho)"
H.tridiagonalize(T,hoho)

print "\nTridiagonal form"
print "T.show()"
T.show()

print "\nHouseholder transformation matrix"
print "hoho.show()"
hoho.show()

print "\nThe original matrix, H, is related to T and hoho matrices as:  H = hoho * T * hoho"
print "H.show()"
H.show()
print "(hoho*T*hoho).show()"
(hoho*T*hoho).show()





