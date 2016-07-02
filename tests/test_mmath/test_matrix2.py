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

import sys
import cmath
import math

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from libra_py import *


N = 5
A = MATRIX(N,N)
B = MATRIX(N,N)
for i in range(0,N):
    for j in range(0,N):
        A.set(i,j,-2.0*i+j*j )
        B.set(i,j,0.0 )

print "The original matrix A"
A.show_matrix()

print "The original matrix B"
B.show_matrix()

a2x2 = MATRIX(2,2)
print "Extracting with stencil [0,1]"
pop_submatrix(A,a2x2, [0,1]); a2x2.show_matrix()

print "Extracting with stencil [3,4]"
pop_submatrix(A,a2x2, [3,4]); a2x2.show_matrix()

print "Extracting with stencil [1,3]"
pop_submatrix(A,a2x2, [1,3]); a2x2.show_matrix()

print "Pushing 2x2 matrix with stencil [0,1]"
push_submatrix(B,a2x2, [0,1]); B.show_matrix()

B = B*0.0

a3x3 = MATRIX(3,3)
print "Extracting 3x3 matrix with stencil [0,1,2]"
pop_submatrix(A,a3x3, [0,1,2]); a3x3.show_matrix()
print "Pushing 3x3 matrix into B with stencil [0,3,4]"
push_submatrix(B,a3x3, [0,3,4]); B.show_matrix()
print "Pushing another, a 2x2, matrix into B with stencil [0,1]"
push_submatrix(B,a2x2, [0,1]); B.show_matrix()



print "The original matrix A"
A.show_matrix()

a2x3 = MATRIX(2,3)
print "Extracting 2x3 matrix with stencils [2,4] and [1,3,4]"
pop_submatrix(A,a2x3, [2,4],[1,3,4]); a2x3.show_matrix()

print "Extracting 2x3 matrix with stencils [0,1] and [2,3,4]"
pop_submatrix(A,a2x3, [0,1],[2,3,4]); a2x3.show_matrix()


a3x2 = MATRIX(3,2)
print "Extracting 3x2 matrix with stencils [2,3,4] and [0,1]"
pop_submatrix(A,a3x2, [2,3,4],[0,1]); a3x2.show_matrix()


B = B*0.0
print "Pushing 2x3 matrix with stencils [3,4] and [0,1,2]"
push_submatrix(B,a2x3, [3,4],[0,1,2]); B.show_matrix()

print "Pushing 3x2 matrix with stencils [0,1,2] and [3,4]"
push_submatrix(B,a3x2, [0,1,2],[3,4]); B.show_matrix()










