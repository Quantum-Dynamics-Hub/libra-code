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

import sys
import cmath
import math

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

#from libra_py import *

A = MATRIX(2,2)
A.set(0,0, 1.0);   A.set(0,1, 2.0);
A.set(1,0, -1.0);  A.set(1,1, 0.5);


print "========== Scaling 1 ============"
B = MATRIX(2,2)
B.set(0,0, 1.0);   B.set(0,1, 1.0);
B.set(1,0, 1.0);   B.set(1,1, 1.0);
C = MATRIX(2,2)
C.dot_product(A,B)

print "Scaling matrix"; B.show_matrix()
print "Original matrix"; A.show_matrix()
print "Result"; C.show_matrix()


print "========== Scaling 2 ============"
B.set(0,0, 1.0);   B.set(0,1, 0.1);
B.set(1,0, 10.0);  B.set(1,1, 1.0);
C.dot_product(A,B)

print "Scaling matrix"; B.show_matrix()
print "Original matrix"; A.show_matrix()
print "Result"; C.show_matrix()


print "========== Scaling 3 ============"
# What if one of matrices is also used to store the result
print "Original matrix"; A.show_matrix()
B.set(0,0, 1.0);   B.set(0,1, 0.1);
B.set(1,0, 10.0);  B.set(1,1, 1.0);
A.dot_product(A,B)

print "Scaling matrix"; B.show_matrix()
print "Result"; A.show_matrix()






