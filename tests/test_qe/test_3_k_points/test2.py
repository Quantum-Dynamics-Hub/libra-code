#*********************************************************************************
#* Copyright (C) 2016 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
## \file test2.py
# Here, we will implement the functions for computing overlaps
# These are Python implementation of the functions to compute the overlap integrals with multiple k-points
# This is VERY SLOW!!!


import os
import sys
import math
import cmath

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *


def I_1D(kx, kxp, gx, gxp):
    # All arguments are float
    # kx,  gx  - refer to k-point 1 and its corresponding grid points
    # kxp, gxp - refer to k-point 2 and its corresponding grid points

    zero = 1e-8

    res = 0.0 + 0.0j

    delt = kx + gx - kxp - gxp

    if math.fabs(delt) <= zero:
        res = 1.0 + 0.0j
    else:
        res = -1.0j *( cmath.exp(2.0j*math.pi*delt) - 1.0+0.0j ) / (2.0 * math.pi*delt)

    return res


def I_3D(k, kp, g, gp):
    # All arguments are VECTOR (float)
    # k,  g  - refer to k-point 1 and its corresponding grid points
    # kp, gp - refer to k-point 2 and its corresponding grid points

    ix = I_1D(k.x, kp.x, g.x, gp.x)
    iy = I_1D(k.y, kp.y, g.y, gp.y)
    iz = I_1D(k.z, kp.z, g.z, gp.z)

    res =  ix * iy * iz

    return res




def overlap(k1, k2, coeff1, coeff2, grid1, grid2):
    # all k- and g-points are in units of 2*pi/a
    # k1, k2 - are k-point vectors (VECTOR of float) 
    # coeff1 - is a matrix (complex) of coefficeints for all states for given k-point (1), dimensions: npw1 x nbands1
    # coeff2 - is a matrix (complex) of coefficeints for all states for given k-point (2), dimensions: npw2 x nbands2
    # grid1 - a list of vectors for all G-points for given k-point (1): dimension npw1
    # grid2 - a list of vectors for all G-points for given k-point (2): dimension npw2

    npw1 = coeff1.num_of_rows
    nbands1 = coeff1.num_of_cols

    npw2 = coeff2.num_of_rows
    nbands2 = coeff2.num_of_cols


    S = CMATRIX(nbands1, nbands2)  # all orbitals for given pair of k-points (a block of entire matrix)

    # A double sum over the grid points (may be different for the two k-points)
    for g1 in xrange(npw1):
        for g2 in xrange(npw2):

            s = I_3D(k1, k2, grid1[g1], grid2[g2])

            for i1 in xrange(nbands1):
                for i2 in xrange(nbands2):

                    tmp = coeff1.get(g1,i1).conjugate() * s * coeff2.get(g2,i2)

                    S.set(i1,i2, S.get(i1,i2) + tmp)

    return S        

    


info, all_e = QE_methods.read_qe_index("x.export/index.xml", [1,2], 1)

print "The total # of k-points is: ", info["nk"]

coeff = []
grid = []

for ik in xrange(info["nk"]):
    print ik, info["k"][ik]

    coeff.append(QE_methods.read_qe_wfc("x.export/wfc.%i" % (ik+1), [1,2], 0))

    grid.append( QE_methods.read_qe_wfc_grid("x.export/grid.%i" % (ik+1) , 0) )


for ik1 in xrange(info["nk"]):
    for ik2 in xrange(info["nk"]):

        S = overlap(info["k"][ik1], info["k"][ik2], coeff[ik1], coeff[ik2], grid[ik1], grid[ik2])

        print ik1, ik2;  S.show_matrix()

sys.exit(0)


### The following way of computing overlaps is INCORRECT for multiple k-points!!! ###

c1 = QE_methods.read_qe_wfc("x.export/wfc.1", [1,2], 1)
c2 = QE_methods.read_qe_wfc("x.export/wfc.2", [1,2], 1)
c3 = QE_methods.read_qe_wfc("x.export/wfc.3", [1,2], 1)


s1  = c1.H() * c1
s2  = c2.H() * c2
s3  = c3.H() * c3
print "s1"; s1.show_matrix()
print "s2"; s2.show_matrix()
print "s3"; s3.show_matrix()

g1 = QE_methods.read_qe_wfc_grid("x.export/grid.1", 1)
g2 = QE_methods.read_qe_wfc_grid("x.export/grid.2", 1)
g3 = QE_methods.read_qe_wfc_grid("x.export/grid.3", 1)



### The following operations will not work, because there may be different number of 
### planewaves needed for each k-point. The proper way to compute the overlaps is 
### via additional matrix with integrals over planewaves

#s12  = c1.H() * c2;  print "s12"; s12.show_matrix()
#s13  = c1.H() * c3;  print "s13"; s13.show_matrix()
#s23  = c2.H() * c3;  print "s23"; s23.show_matrix()



print "\n"


