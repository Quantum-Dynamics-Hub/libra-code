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


t = Timer()
t.start()
print "\nTest 2: 3D overlaps: normalized"
f = open("3D_overlaps_norm.txt","w")
Ra = VECTOR(0.0, 0.0, 0.0)
Rb = VECTOR(0.0, 0.0, 0.0)
u = VECTOR(1.0, 2.0, -0.5)
u.normalize()
print "u (normalized) = ", u.x, u.y, u.z


for i in range(0,50):
    x = 0.1 * i
    Rb = Ra + x*u    

    ss  = gaussian_overlap(0,0,0, 1.3, Ra,  0,0,0, 1.3, Rb )  # s(H)-s(H)
    spx = gaussian_overlap(0,0,0, 1.3, Ra,  1,0,0, 1.3, Rb )  # s(H)-px(H)
    spy = gaussian_overlap(0,0,0, 1.3, Ra,  0,1,0, 1.3, Rb )  # s(H)-py(H)
    spz = gaussian_overlap(0,0,0, 1.3, Ra,  0,0,1, 1.3, Rb )  # s(H)-pz(H)
    pxpx = gaussian_overlap(1,0,0, 1.3, Ra,  1,0,0, 1.3, Rb )  # px(H)-px(H) 
    pxpy = gaussian_overlap(1,0,0, 1.3, Ra,  0,1,0, 1.3, Rb )  # px(H)-py(H) = 0
    pypz = gaussian_overlap(0,1,0, 1.3, Ra,  0,0,1, 1.3, Rb )  # py(H)-pz(H) = 0
    sdz2 = gaussian_overlap(0,0,0, 1.3, Ra,  0,0,2, 1.3, Rb )  # s(H)-dz2(H) 
    sdxy = gaussian_overlap(0,0,0, 1.3, Ra,  0,1,1, 1.3, Rb )  # s(H)-dxy(H) 
    pxdxy = gaussian_overlap(1,0,0, 1.3, Ra,  0,1,1, 1.3, Rb )  # s(H)-dxy(H) 

    f.write("%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f \n" % (x, ss, spx, spy, spz, pxpx, pxpy, pypz, sdz2, sdxy, pxdxy) )


f.close()

t.stop()
print "Time to compute = ", t.show(), " sec"



