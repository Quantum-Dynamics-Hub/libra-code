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


is_norm = 0

t = Timer()
t.start()
print "\nTest 2: 1D overlaps: Analytic derivatives"
f = open("1D_overlaps_der_anal.txt","w")
for i in range(0,50):
    x = 0.1 * i
    # Derivatives
    ss = gaussian_overlap(0, 1.3, 0.0, 0, 1.3, x, is_norm, 1 )  # s(H)-s(H)
    sp = gaussian_overlap(0, 1.3, 0.0, 1, 1.3, x, is_norm, 1 )  # s(H)-p(H)
    sd = gaussian_overlap(0, 1.3, 0.0, 2, 1.3, x, is_norm, 1 )  # s(H)-d(H)
    sf = gaussian_overlap(0, 1.3, 0.0, 3, 1.3, x, is_norm, 1 )  # s(H)-f(H)
    pp = gaussian_overlap(1, 1.3, 0.0, 1, 1.3, x, is_norm, 1 )  # p(H)-p(H)
    pd = gaussian_overlap(1, 1.3, 0.0, 2, 1.3, x, is_norm, 1 )  # p(H)-d(H)
    pf = gaussian_overlap(1, 1.3, 0.0, 3, 1.3, x, is_norm, 1 )  # p(H)-f(H)
    dd = gaussian_overlap(2, 1.3, 0.0, 2, 1.3, x, is_norm, 1 )  # d(H)-d(H)
    df = gaussian_overlap(2, 1.3, 0.0, 3, 1.3, x, is_norm, 1 )  # d(H)-f(H)
    ff = gaussian_overlap(3, 1.3, 0.0, 3, 1.3, x, is_norm, 1 )  # f(H)-f(H)

    f.write("%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f \n" % (x, ss[1], sp[1], sd[1], sf[1], pp[1], pd[1], pf[1], dd[1], df[1], ff[1]) )

f.close()

t.stop()
print "Time to compute = ", t.show(), " sec"


dx = 0.0001
t.start()
print "\nTest 3: 1D overlaps: Numerical derivatives"
f = open("1D_overlaps_der_num.txt","w")
for i in range(0,50):
    x = 0.1 * i
    ss = (gaussian_overlap(0, 1.3, 0.5*dx, 0, 1.3, x, is_norm )  - gaussian_overlap(0, 1.3, -0.5*dx, 0, 1.3, x, is_norm ) )/dx
    sp = (gaussian_overlap(0, 1.3, 0.5*dx, 1, 1.3, x, is_norm )  - gaussian_overlap(0, 1.3, -0.5*dx, 1, 1.3, x, is_norm ) )/dx
    sd = (gaussian_overlap(0, 1.3, 0.5*dx, 2, 1.3, x, is_norm )  - gaussian_overlap(0, 1.3, -0.5*dx, 2, 1.3, x, is_norm ) )/dx
    sf = (gaussian_overlap(0, 1.3, 0.5*dx, 3, 1.3, x, is_norm )  - gaussian_overlap(0, 1.3, -0.5*dx, 3, 1.3, x, is_norm ) )/dx
    pp = (gaussian_overlap(1, 1.3, 0.5*dx, 1, 1.3, x, is_norm )  - gaussian_overlap(1, 1.3, -0.5*dx, 1, 1.3, x, is_norm ) )/dx
    pd = (gaussian_overlap(1, 1.3, 0.5*dx, 2, 1.3, x, is_norm )  - gaussian_overlap(1, 1.3, -0.5*dx, 2, 1.3, x, is_norm ) )/dx
    pf = (gaussian_overlap(1, 1.3, 0.5*dx, 3, 1.3, x, is_norm )  - gaussian_overlap(1, 1.3, -0.5*dx, 3, 1.3, x, is_norm ) )/dx
    dd = (gaussian_overlap(2, 1.3, 0.5*dx, 2, 1.3, x, is_norm )  - gaussian_overlap(2, 1.3, -0.5*dx, 2, 1.3, x, is_norm ) )/dx
    df = (gaussian_overlap(2, 1.3, 0.5*dx, 3, 1.3, x, is_norm )  - gaussian_overlap(2, 1.3, -0.5*dx, 3, 1.3, x, is_norm ) )/dx
    ff = (gaussian_overlap(3, 1.3, 0.5*dx, 3, 1.3, x, is_norm )  - gaussian_overlap(3, 1.3, -0.5*dx, 3, 1.3, x, is_norm ) )/dx

    f.write("%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f \n" % (x, ss, sp, sd, sf, pp, pd, pf, dd, df, ff) )

f.close()

t.stop()
print "Time to compute = ", t.show(), " sec"



