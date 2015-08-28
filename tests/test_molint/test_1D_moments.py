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
sys.path.insert(1,cwd+"/../../_build/src/qchem")


print "\nTest 1: Importing the library and its content"
from cygmmath import *
from cygqchem import *



print "\nTest 2: 1D moments (normalized)"
print "Normalization coefficients"
s_norm = gaussian_norm(0, 1.3)
p_norm = gaussian_norm(1, 1.3)
d_norm = gaussian_norm(2, 1.3)
f_norm = gaussian_norm(3, 1.3)
print "s_norm = ", s_norm, "|<s|s>|^2 = ", gaussian_overlap(0, 1.3, 0.0,   0, 1.3, 0.0, 0 )*(s_norm**2)
print "p_norm = ", p_norm, "|<p|p>|^2 = ", gaussian_overlap(1, 1.3, 0.0,   1, 1.3, 0.0, 0 )*(p_norm**2)
print "d_norm = ", d_norm, "|<d|d>|^2 = ", gaussian_overlap(2, 1.3, 0.0,   2, 1.3, 0.0, 0 )*(d_norm**2)
print "f_norm = ", f_norm, "|<f|f>|^2 = ", gaussian_overlap(3, 1.3, 0.0,   3, 1.3, 0.0, 0 )*(f_norm**2)

f = open("1D_moments_norm.txt","w")
for i in range(0,50):
    x = 0.1 * i

    ss1 = gaussian_moment(0, 1.3, 0.0,  0, 0.0, 0.0,   0, 1.3, x )  # <s(A)| 1 | s(B)>
    ss2 = gaussian_moment(0, 1.3, 0.0,  1, 0.0, 0.0,   0, 1.3, x )  # <s(A)| (x-x(A))| s(B) >
    ss3 = gaussian_moment(0, 1.3, 0.0,  1, 0.0, 0.5*x, 0, 1.3, x )  # <s(A)| (x-(x(A)+x(B))/2) | s(B)>
    ss4 = gaussian_moment(0, 1.3, 0.0,  2, 0.0, 0.0,   0, 1.3, x )  # <s(A)| (x-X(A))^2 | s(B)>
                                                       
    pp1 = gaussian_moment(1, 1.3, 0.0,  0, 0.0, 0.0,   1, 1.3, x )  # <p(A)| 1 | p(B)>
    pp2 = gaussian_moment(1, 1.3, 0.0,  1, 0.0, 0.0,   1, 1.3, x )  # <p(A)| (x-x(A))| s(B) >
    pp3 = gaussian_moment(1, 1.3, 0.0,  1, 0.0, 0.5*x, 1, 1.3, x )  # <p(A)| (x-(x(A)+x(B))/2) | p(B)>
    pp4 = gaussian_moment(1, 1.3, 0.0,  2, 0.0, 0.0,   1, 1.3, x )  # <p(A)| (x-X(A))^2 | p(B)>
                                                       
    sp1 = gaussian_moment(0, 1.3, 0.0,  0, 0.0, 0.0,   1, 1.3, x )  # <s(A)| 1 | p(B)>
    sp2 = gaussian_moment(0, 1.3, 0.0,  1, 0.0, 0.0,   1, 1.3, x )  # <s(A)| (x-x(A))| s(B) >
    sp3 = gaussian_moment(0, 1.3, 0.0,  1, 0.0, 0.5*x, 1, 1.3, x )  # <s(A)| (x-(x(A)+x(B))/2) | p(B)>
    sp4 = gaussian_moment(0, 1.3, 0.0,  2, 0.0, 0.0,   1, 1.3, x )  # <s(A)| (x-X(A))^2 | p(B)>



    f.write("%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f \n" % (x, ss1, ss2, ss3, ss4, pp1, pp2, pp3, pp4, sp1, sp2, sp3, sp4) )


f.close()


