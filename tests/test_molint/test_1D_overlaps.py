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


t = Timer()
t.start()
print "\nTest 2: 1D overlaps: Unnormalized gaussians"
f = open("1D_overlaps_ref.txt","w")
for i in range(0,50):
    x = 0.1 * i
    ss = gaussian_overlap_ref(0, 1.3, 0.0, 0, 1.3, x )  # s(H)-s(H)
    sp = gaussian_overlap_ref(0, 1.3, 0.0, 1, 1.3, x )  # s(H)-p(H)
    sd = gaussian_overlap_ref(0, 1.3, 0.0, 2, 1.3, x )  # s(H)-d(H)
    sf = gaussian_overlap_ref(0, 1.3, 0.0, 3, 1.3, x )  # s(H)-f(H)
    pp = gaussian_overlap_ref(1, 1.3, 0.0, 1, 1.3, x )  # p(H)-p(H)
    pd = gaussian_overlap_ref(1, 1.3, 0.0, 2, 1.3, x )  # p(H)-d(H)
    pf = gaussian_overlap_ref(1, 1.3, 0.0, 3, 1.3, x )  # p(H)-f(H)
    dd = gaussian_overlap_ref(2, 1.3, 0.0, 2, 1.3, x )  # d(H)-d(H)
    df = gaussian_overlap_ref(2, 1.3, 0.0, 3, 1.3, x )  # d(H)-f(H)
    ff = gaussian_overlap_ref(3, 1.3, 0.0, 3, 1.3, x )  # f(H)-f(H)

    f.write("%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f \n" % (x, ss, sp, sd, sf, pp, pd, pf, dd, df, ff) )

f.close()

t.stop()
print "Time to compute = ", t.show(), " sec"



t.start()
print "\nTest 3: 1D overlaps: Fast version"
f = open("1D_overlaps.txt","w")
for i in range(0,50):
    x = 0.1 * i
#    print x, binomial_expansion(0, 1, 0.5*x, -0.5*x, 0)
    ss = gaussian_overlap(0, 1.3, 0.0, 0, 1.3, x, 0 )  # s(H)-s(H)
    sp = gaussian_overlap(0, 1.3, 0.0, 1, 1.3, x, 0 )  # s(H)-p(H)
    sd = gaussian_overlap(0, 1.3, 0.0, 2, 1.3, x, 0 )  # s(H)-d(H)
    sf = gaussian_overlap(0, 1.3, 0.0, 3, 1.3, x, 0 )  # s(H)-f(H)
    pp = gaussian_overlap(1, 1.3, 0.0, 1, 1.3, x, 0 )  # p(H)-p(H)
    pd = gaussian_overlap(1, 1.3, 0.0, 2, 1.3, x, 0 )  # p(H)-d(H)
    pf = gaussian_overlap(1, 1.3, 0.0, 3, 1.3, x, 0 )  # p(H)-f(H)
    dd = gaussian_overlap(2, 1.3, 0.0, 2, 1.3, x, 0 )  # d(H)-d(H)
    df = gaussian_overlap(2, 1.3, 0.0, 3, 1.3, x, 0 )  # d(H)-f(H)
    ff = gaussian_overlap(3, 1.3, 0.0, 3, 1.3, x, 0 )  # f(H)-f(H)

    f.write("%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f \n" % (x, ss, sp, sd, sf, pp, pd, pf, dd, df, ff) )

f.close()

t.stop()
print "Time to compute = ", t.show(), " sec"



