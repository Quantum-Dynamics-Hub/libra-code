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



print "\nTest 2: 1D overlaps: normalized STOs"  # similar to Test 3
print "These are actually normalization coefficients"
Rcut = 12.0
s_norm = sto_norm(0, 1.3)
p_norm = sto_norm(1, 1.3)
d_norm = sto_norm(2, 1.3)
f_norm = sto_norm(3, 1.3)
print "sto_overlap  returns overlaps assuming STOs are normalized already!"
print "s_norm = ", s_norm, "|<s|s>|^2 = ", sto_overlap(0, 0, 0, 1.3,  0, 0, 0, 1.3,  0.0, Rcut )#*(s_norm**2)
print "p_norm = ", p_norm, "|<p|p>|^2 = ", sto_overlap(1, 0, 0, 1.3,  1, 0, 0, 1.3,  0.0, Rcut )#*(p_norm**2)
print "d_norm = ", d_norm, "|<d|d>|^2 = ", sto_overlap(2, 0, 0, 1.3,  2, 0, 0, 1.3,  0.0, Rcut )#*(d_norm**2)
print "f_norm = ", f_norm, "|<f|f>|^2 = ", sto_overlap(3, 0, 0, 1.3,  3, 0, 0, 1.3,  0.0, Rcut )#*(f_norm**2)

f = open("1D_sto_overlaps_norm.txt","w")
for i in range(0,50):
    x = 0.1 * i

    ss = sto_overlap(0, 0, 0, 1.3,  0, 0, 0, 1.3,  x, Rcut )
    sp = sto_overlap(0, 0, 0, 1.3,  1, 0, 0, 1.3,  x, Rcut )
    sd = sto_overlap(0, 0, 0, 1.3,  2, 0, 0, 1.3,  x, Rcut )
    sf = sto_overlap(0, 0, 0, 1.3,  3, 0, 0, 1.3,  x, Rcut )
    pp = sto_overlap(1, 0, 0, 1.3,  1, 0, 0, 1.3,  x, Rcut )
    pd = sto_overlap(1, 0, 0, 1.3,  2, 0, 0, 1.3,  x, Rcut )
    pf = sto_overlap(1, 0, 0, 1.3,  3, 0, 0, 1.3,  x, Rcut )
    dd = sto_overlap(2, 0, 0, 1.3,  2, 0, 0, 1.3,  x, Rcut )
    df = sto_overlap(2, 0, 0, 1.3,  3, 0, 0, 1.3,  x, Rcut )
    ff = sto_overlap(3, 0, 0, 1.3,  3, 0, 0, 1.3,  x, Rcut )


    f.write("%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f \n" % (x, ss, sp, sd, sf, pp, pd, pf, dd, df, ff) )


f.close()


f = open("1D_sto_overlaps_norm1.txt","w")
for i in range(0,50):
    x = 0.1 * i

    ss = sto_overlap_fast(0, 0, 0, 1.3,  0, 0, 0, 1.3,  x, Rcut )
    sp = sto_overlap_fast(0, 0, 0, 1.3,  1, 0, 0, 1.3,  x, Rcut )
    sd = sto_overlap_fast(0, 0, 0, 1.3,  2, 0, 0, 1.3,  x, Rcut )
    sf = sto_overlap_fast(0, 0, 0, 1.3,  3, 0, 0, 1.3,  x, Rcut )
    pp = sto_overlap_fast(1, 0, 0, 1.3,  1, 0, 0, 1.3,  x, Rcut )
    pd = sto_overlap_fast(1, 0, 0, 1.3,  2, 0, 0, 1.3,  x, Rcut )
    pf = sto_overlap_fast(1, 0, 0, 1.3,  3, 0, 0, 1.3,  x, Rcut )
    dd = sto_overlap_fast(2, 0, 0, 1.3,  2, 0, 0, 1.3,  x, Rcut )
    df = sto_overlap_fast(2, 0, 0, 1.3,  3, 0, 0, 1.3,  x, Rcut )
    ff = sto_overlap_fast(3, 0, 0, 1.3,  3, 0, 0, 1.3,  x, Rcut )


    f.write("%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f \n" % (x, ss, sp, sd, sf, pp, pd, pf, dd, df, ff) )


f.close()





