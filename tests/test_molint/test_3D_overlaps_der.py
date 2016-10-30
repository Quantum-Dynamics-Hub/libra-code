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


is_norm = 1


t = Timer()
t.start()
print "\nTest 2: 3D overlaps: Analytic derivatives"
Ra = VECTOR(0.0, 0.0, 0.0)
Rb = VECTOR(0.0, 0.0, 0.0)
u = VECTOR(1.0, 2.0, -0.5)
u.normalize()
print "u (normalized) = ", u.x, u.y, u.z


f_x = open("3D_overlaps_der_anal_x.txt","w")
f_y = open("3D_overlaps_der_anal_y.txt","w")
f_z = open("3D_overlaps_der_anal_z.txt","w")

for i in range(0,50):
    x = 0.1 * i
    Rb = Ra + x*u    
    ss  = gaussian_overlap(0,0,0, 1.3, Ra,  0,0,0, 1.3, Rb, is_norm, 1)  # s(H)-s(H)
    spx = gaussian_overlap(0,0,0, 1.3, Ra,  1,0,0, 1.3, Rb, is_norm, 1 )  # s(H)-px(H)
    spy = gaussian_overlap(0,0,0, 1.3, Ra,  0,1,0, 1.3, Rb, is_norm, 1 )  # s(H)-py(H)
    spz = gaussian_overlap(0,0,0, 1.3, Ra,  0,0,1, 1.3, Rb, is_norm, 1 )  # s(H)-pz(H)
    pxpx = gaussian_overlap(1,0,0, 1.3, Ra,  1,0,0, 1.3, Rb, is_norm, 1 )  # px(H)-px(H) 
    pxpy = gaussian_overlap(1,0,0, 1.3, Ra,  0,1,0, 1.3, Rb, is_norm, 1 )  # px(H)-py(H) = 0
    pypz = gaussian_overlap(0,1,0, 1.3, Ra,  0,0,1, 1.3, Rb, is_norm, 1 )  # py(H)-pz(H) = 0
    sdz2 = gaussian_overlap(0,0,0, 1.3, Ra,  0,0,2, 1.3, Rb, is_norm, 1 )  # s(H)-dz2(H) 
    sdxy = gaussian_overlap(0,0,0, 1.3, Ra,  0,1,1, 1.3, Rb, is_norm, 1 )  # s(H)-dxy(H) 
    pxdxy = gaussian_overlap(1,0,0, 1.3, Ra,  0,1,1, 1.3, Rb, is_norm, 1 )  # s(H)-dxy(H) 

    f_x.write("%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f \n" % (x, ss[2].x, spx[2].x, spy[2].x, spz[2].x, pxpx[2].x, pxpy[2].x, pypz[2].x, sdz2[2].x, sdxy[2].x, pxdxy[2].x) )
    f_y.write("%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f \n" % (x, ss[2].y, spx[2].y, spy[2].y, spz[2].y, pxpx[2].y, pxpy[2].y, pypz[2].y, sdz2[2].y, sdxy[2].y, pxdxy[2].y) )
    f_z.write("%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f \n" % (x, ss[2].z, spx[2].z, spy[2].z, spz[2].z, pxpx[2].z, pxpy[2].z, pypz[2].z, sdz2[2].z, sdxy[2].z, pxdxy[2].z) )


f_x.close()
f_y.close()
f_z.close()

t.stop()
print "Time to compute = ", t.show(), " sec"


Rbl = VECTOR(0.0, 0.0, 0.0)
Rbr = VECTOR(0.0, 0.0, 0.0)

dx = 0.0001
f_x = open("3D_overlaps_der_num_x.txt","w")
f_y = open("3D_overlaps_der_num_y.txt","w")
f_z = open("3D_overlaps_der_num_z.txt","w")
for i in range(0,50):
    x = 0.1 * i

    # Numerical derivatives w.r.t. coordinate B along x direction
    Rbl = Ra + x*u + VECTOR(-0.5*dx, 0.0, 0.0)
    Rbr = Ra + x*u + VECTOR( 0.5*dx, 0.0, 0.0)  

    ss    = (gaussian_overlap(0,0,0, 1.3, Ra,  0,0,0, 1.3, Rbr, is_norm) - gaussian_overlap(0,0,0, 1.3, Ra,  0,0,0, 1.3, Rbl, is_norm) )/dx
    spx   = (gaussian_overlap(0,0,0, 1.3, Ra,  1,0,0, 1.3, Rbr, is_norm) - gaussian_overlap(0,0,0, 1.3, Ra,  1,0,0, 1.3, Rbl, is_norm) )/dx
    spy   = (gaussian_overlap(0,0,0, 1.3, Ra,  0,1,0, 1.3, Rbr, is_norm) - gaussian_overlap(0,0,0, 1.3, Ra,  0,1,0, 1.3, Rbl, is_norm) )/dx
    spz   = (gaussian_overlap(0,0,0, 1.3, Ra,  0,0,1, 1.3, Rbr, is_norm) - gaussian_overlap(0,0,0, 1.3, Ra,  0,0,1, 1.3, Rbl, is_norm) )/dx
    pxpx  = (gaussian_overlap(1,0,0, 1.3, Ra,  1,0,0, 1.3, Rbr, is_norm) - gaussian_overlap(1,0,0, 1.3, Ra,  1,0,0, 1.3, Rbl, is_norm) )/dx
    pxpy  = (gaussian_overlap(1,0,0, 1.3, Ra,  0,1,0, 1.3, Rbr, is_norm) - gaussian_overlap(1,0,0, 1.3, Ra,  0,1,0, 1.3, Rbl, is_norm) )/dx
    pypz  = (gaussian_overlap(0,1,0, 1.3, Ra,  0,0,1, 1.3, Rbr, is_norm) - gaussian_overlap(0,1,0, 1.3, Ra,  0,0,1, 1.3, Rbl, is_norm) )/dx
    sdz2  = (gaussian_overlap(0,0,0, 1.3, Ra,  0,0,2, 1.3, Rbr, is_norm) - gaussian_overlap(0,0,0, 1.3, Ra,  0,0,2, 1.3, Rbl, is_norm) )/dx
    sdxy  = (gaussian_overlap(0,0,0, 1.3, Ra,  0,1,1, 1.3, Rbr, is_norm) - gaussian_overlap(0,0,0, 1.3, Ra,  0,1,1, 1.3, Rbl, is_norm) )/dx
    pxdxy = (gaussian_overlap(1,0,0, 1.3, Ra,  0,1,1, 1.3, Rbr, is_norm) - gaussian_overlap(1,0,0, 1.3, Ra,  0,1,1, 1.3, Rbl, is_norm) )/dx

    f_x.write("%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f \n" % (x, ss, spx, spy, spz, pxpx, pxpy, pypz, sdz2, sdxy, pxdxy) )



    # Numerical derivatives w.r.t. coordinate B along y direction
    Rbl = Ra + x*u + VECTOR(0.0, -0.5*dx, 0.0)
    Rbr = Ra + x*u + VECTOR(0.0,  0.5*dx, 0.0)  

    ss    = (gaussian_overlap(0,0,0, 1.3, Ra,  0,0,0, 1.3, Rbr, is_norm) - gaussian_overlap(0,0,0, 1.3, Ra,  0,0,0, 1.3, Rbl, is_norm) )/dx
    spx   = (gaussian_overlap(0,0,0, 1.3, Ra,  1,0,0, 1.3, Rbr, is_norm) - gaussian_overlap(0,0,0, 1.3, Ra,  1,0,0, 1.3, Rbl, is_norm) )/dx
    spy   = (gaussian_overlap(0,0,0, 1.3, Ra,  0,1,0, 1.3, Rbr, is_norm) - gaussian_overlap(0,0,0, 1.3, Ra,  0,1,0, 1.3, Rbl, is_norm) )/dx
    spz   = (gaussian_overlap(0,0,0, 1.3, Ra,  0,0,1, 1.3, Rbr, is_norm) - gaussian_overlap(0,0,0, 1.3, Ra,  0,0,1, 1.3, Rbl, is_norm) )/dx
    pxpx  = (gaussian_overlap(1,0,0, 1.3, Ra,  1,0,0, 1.3, Rbr, is_norm) - gaussian_overlap(1,0,0, 1.3, Ra,  1,0,0, 1.3, Rbl, is_norm) )/dx
    pxpy  = (gaussian_overlap(1,0,0, 1.3, Ra,  0,1,0, 1.3, Rbr, is_norm) - gaussian_overlap(1,0,0, 1.3, Ra,  0,1,0, 1.3, Rbl, is_norm) )/dx
    pypz  = (gaussian_overlap(0,1,0, 1.3, Ra,  0,0,1, 1.3, Rbr, is_norm) - gaussian_overlap(0,1,0, 1.3, Ra,  0,0,1, 1.3, Rbl, is_norm) )/dx
    sdz2  = (gaussian_overlap(0,0,0, 1.3, Ra,  0,0,2, 1.3, Rbr, is_norm) - gaussian_overlap(0,0,0, 1.3, Ra,  0,0,2, 1.3, Rbl, is_norm) )/dx
    sdxy  = (gaussian_overlap(0,0,0, 1.3, Ra,  0,1,1, 1.3, Rbr, is_norm) - gaussian_overlap(0,0,0, 1.3, Ra,  0,1,1, 1.3, Rbl, is_norm) )/dx
    pxdxy = (gaussian_overlap(1,0,0, 1.3, Ra,  0,1,1, 1.3, Rbr, is_norm) - gaussian_overlap(1,0,0, 1.3, Ra,  0,1,1, 1.3, Rbl, is_norm) )/dx

    f_y.write("%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f \n" % (x, ss, spx, spy, spz, pxpx, pxpy, pypz, sdz2, sdxy, pxdxy) )



    # Numerical derivatives w.r.t. coordinate B along z direction
    Rbl = Ra + x*u + VECTOR(0.0, 0.0, -0.5*dx)
    Rbr = Ra + x*u + VECTOR(0.0, 0.0,  0.5*dx)  

    ss    = (gaussian_overlap(0,0,0, 1.3, Ra,  0,0,0, 1.3, Rbr, is_norm) - gaussian_overlap(0,0,0, 1.3, Ra,  0,0,0, 1.3, Rbl, is_norm) )/dx
    spx   = (gaussian_overlap(0,0,0, 1.3, Ra,  1,0,0, 1.3, Rbr, is_norm) - gaussian_overlap(0,0,0, 1.3, Ra,  1,0,0, 1.3, Rbl, is_norm) )/dx
    spy   = (gaussian_overlap(0,0,0, 1.3, Ra,  0,1,0, 1.3, Rbr, is_norm) - gaussian_overlap(0,0,0, 1.3, Ra,  0,1,0, 1.3, Rbl, is_norm) )/dx
    spz   = (gaussian_overlap(0,0,0, 1.3, Ra,  0,0,1, 1.3, Rbr, is_norm) - gaussian_overlap(0,0,0, 1.3, Ra,  0,0,1, 1.3, Rbl, is_norm) )/dx
    pxpx  = (gaussian_overlap(1,0,0, 1.3, Ra,  1,0,0, 1.3, Rbr, is_norm) - gaussian_overlap(1,0,0, 1.3, Ra,  1,0,0, 1.3, Rbl, is_norm) )/dx
    pxpy  = (gaussian_overlap(1,0,0, 1.3, Ra,  0,1,0, 1.3, Rbr, is_norm) - gaussian_overlap(1,0,0, 1.3, Ra,  0,1,0, 1.3, Rbl, is_norm) )/dx
    pypz  = (gaussian_overlap(0,1,0, 1.3, Ra,  0,0,1, 1.3, Rbr, is_norm) - gaussian_overlap(0,1,0, 1.3, Ra,  0,0,1, 1.3, Rbl, is_norm) )/dx
    sdz2  = (gaussian_overlap(0,0,0, 1.3, Ra,  0,0,2, 1.3, Rbr, is_norm) - gaussian_overlap(0,0,0, 1.3, Ra,  0,0,2, 1.3, Rbl, is_norm) )/dx
    sdxy  = (gaussian_overlap(0,0,0, 1.3, Ra,  0,1,1, 1.3, Rbr, is_norm) - gaussian_overlap(0,0,0, 1.3, Ra,  0,1,1, 1.3, Rbl, is_norm) )/dx
    pxdxy = (gaussian_overlap(1,0,0, 1.3, Ra,  0,1,1, 1.3, Rbr, is_norm) - gaussian_overlap(1,0,0, 1.3, Ra,  0,1,1, 1.3, Rbl, is_norm) )/dx


    f_z.write("%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f \n" % (x, ss, spx, spy, spz, pxpx, pxpy, pypz, sdz2, sdxy, pxdxy) )
 


f_x.close()
f_y.close()
f_z.close()




