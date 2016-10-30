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



print "\nTest 2: Pseudopotentials with 3D Gaussians (normalized): Analytic derivatives"

Ra = VECTOR(0.0, 0.0, 0.0)
Rb = VECTOR(0.0, 0.0, 0.0)
Rc = VECTOR(0.0, 0.0, 0.0)

u = VECTOR(1.0, 2.0, -0.5)
u.normalize()
print "u (normalized) = ", u.x, u.y, u.z



f_Xx = open("3D_pseudopot02_der_anal_Xx.txt","w")
f_Xy = open("3D_pseudopot02_der_anal_Xy.txt","w")
f_Xz = open("3D_pseudopot02_der_anal_Xz.txt","w")

f_Ax = open("3D_pseudopot02_der_anal_Ax.txt","w")
f_Ay = open("3D_pseudopot02_der_anal_Ay.txt","w")
f_Az = open("3D_pseudopot02_der_anal_Az.txt","w")


for i in range(0,50):
    x = 0.1 * i
    Rb = Ra + x*u    

    ss_pp   = pseudopot02(1.0, 1.0, 2.0,Rc,   0,0,0, 1.3, Ra,  0,0,0, 1.3, Rb, 1, 1 )  # s(H)-s(H)
    spx_pp  = pseudopot02(1.0, 1.0, 2.0,Rc,   0,0,0, 1.3, Ra,  1,0,0, 1.3, Rb, 1, 1 )  # s(H)-px(H)
    spy_pp  = pseudopot02(1.0, 1.0, 2.0,Rc,   0,0,0, 1.3, Ra,  0,1,0, 1.3, Rb, 1, 1 )  # s(H)-py(H)
    pxpx_pp = pseudopot02(1.0, 1.0, 2.0,Rc,   1,0,0, 1.3, Ra,  1,0,0, 1.3, Rb, 1, 1 )  # px(H)-px(H)
    pxpy_pp = pseudopot02(1.0, 1.0, 2.0,Rc,   1,0,0, 1.3, Ra,  0,1,0, 1.3, Rb, 1, 1 )  # px(H)-py(H)


    f_Xx.write("%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f \n" % (x, ss_pp[1].x, spx_pp[1].x, spy_pp[1].x, pxpx_pp[1].x, pxpy_pp[1].x ) )
    f_Xy.write("%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f \n" % (x, ss_pp[1].y, spx_pp[1].y, spy_pp[1].y, pxpx_pp[1].y, pxpy_pp[1].y ) )
    f_Xz.write("%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f \n" % (x, ss_pp[1].z, spx_pp[1].z, spy_pp[1].z, pxpx_pp[1].z, pxpy_pp[1].z ) )

    f_Ax.write("%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f \n" % (x, ss_pp[2].x, spx_pp[2].x, spy_pp[2].x, pxpx_pp[2].x, pxpy_pp[2].x ) )
    f_Ay.write("%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f \n" % (x, ss_pp[2].y, spx_pp[2].y, spy_pp[2].y, pxpx_pp[2].y, pxpy_pp[2].y ) )
    f_Az.write("%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f \n" % (x, ss_pp[2].z, spx_pp[2].z, spy_pp[2].z, pxpx_pp[2].z, pxpy_pp[2].z ) )


f_Xx.close()
f_Xy.close()
f_Xz.close()
f_Ax.close()
f_Ay.close()
f_Az.close()




print "\nTest 3: Pseudopotentials with 3D Gaussians (normalized): Numerical derivatives"

Ra = VECTOR(0.0, 0.0, 0.0)
Rb = VECTOR(0.0, 0.0, 0.0)
Rc = VECTOR(0.0, 0.0, 0.0)

u = VECTOR(1.0, 2.0, -0.5)
u.normalize()
print "u (normalized) = ", u.x, u.y, u.z



f_Xx = open("3D_pseudopot02_der_num_Xx.txt","w")
f_Xy = open("3D_pseudopot02_der_num_Xy.txt","w")
f_Xz = open("3D_pseudopot02_der_num_Xz.txt","w")

f_Ax = open("3D_pseudopot02_der_num_Ax.txt","w")
f_Ay = open("3D_pseudopot02_der_num_Ay.txt","w")
f_Az = open("3D_pseudopot02_der_num_Az.txt","w")

dx = 0.0001
drx = VECTOR(0.5*dx, 0.0, 0.0)
dry = VECTOR(0.0, 0.5*dx, 0.0)
drz = VECTOR(0.0, 0.0, 0.5*dx)

for i in range(0,50):
    x = 0.1 * i
    Rb = Ra + x*u    

    # d/dXx :
    ss_pp   = (pseudopot02(1.0, 1.0, 2.0,Rc + drx,   0,0,0, 1.3, Ra,  0,0,0, 1.3, Rb ) - pseudopot02(1.0, 1.0, 2.0,Rc - drx,   0,0,0, 1.3, Ra,  0,0,0, 1.3, Rb ) )/dx
    spx_pp  = (pseudopot02(1.0, 1.0, 2.0,Rc + drx,   0,0,0, 1.3, Ra,  1,0,0, 1.3, Rb ) - pseudopot02(1.0, 1.0, 2.0,Rc - drx,   0,0,0, 1.3, Ra,  1,0,0, 1.3, Rb ) )/dx
    spy_pp  = (pseudopot02(1.0, 1.0, 2.0,Rc + drx,   0,0,0, 1.3, Ra,  0,1,0, 1.3, Rb ) - pseudopot02(1.0, 1.0, 2.0,Rc - drx,   0,0,0, 1.3, Ra,  0,1,0, 1.3, Rb ) )/dx
    pxpx_pp = (pseudopot02(1.0, 1.0, 2.0,Rc + drx,   1,0,0, 1.3, Ra,  1,0,0, 1.3, Rb ) - pseudopot02(1.0, 1.0, 2.0,Rc - drx,   1,0,0, 1.3, Ra,  1,0,0, 1.3, Rb ) )/dx
    pxpy_pp = (pseudopot02(1.0, 1.0, 2.0,Rc + drx,   1,0,0, 1.3, Ra,  0,1,0, 1.3, Rb ) - pseudopot02(1.0, 1.0, 2.0,Rc - drx,   1,0,0, 1.3, Ra,  0,1,0, 1.3, Rb ) )/dx
                                          
    f_Xx.write("%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f \n" % (x, ss_pp, spx_pp, spy_pp, pxpx_pp, pxpy_pp ) )

    # d/dXy :
    ss_pp   = (pseudopot02(1.0, 1.0, 2.0,Rc + dry,   0,0,0, 1.3, Ra,  0,0,0, 1.3, Rb ) - pseudopot02(1.0, 1.0, 2.0,Rc - dry,   0,0,0, 1.3, Ra,  0,0,0, 1.3, Rb ) )/dx
    spx_pp  = (pseudopot02(1.0, 1.0, 2.0,Rc + dry,   0,0,0, 1.3, Ra,  1,0,0, 1.3, Rb ) - pseudopot02(1.0, 1.0, 2.0,Rc - dry,   0,0,0, 1.3, Ra,  1,0,0, 1.3, Rb ) )/dx
    spy_pp  = (pseudopot02(1.0, 1.0, 2.0,Rc + dry,   0,0,0, 1.3, Ra,  0,1,0, 1.3, Rb ) - pseudopot02(1.0, 1.0, 2.0,Rc - dry,   0,0,0, 1.3, Ra,  0,1,0, 1.3, Rb ) )/dx
    pxpx_pp = (pseudopot02(1.0, 1.0, 2.0,Rc + dry,   1,0,0, 1.3, Ra,  1,0,0, 1.3, Rb ) - pseudopot02(1.0, 1.0, 2.0,Rc - dry,   1,0,0, 1.3, Ra,  1,0,0, 1.3, Rb ) )/dx
    pxpy_pp = (pseudopot02(1.0, 1.0, 2.0,Rc + dry,   1,0,0, 1.3, Ra,  0,1,0, 1.3, Rb ) - pseudopot02(1.0, 1.0, 2.0,Rc - dry,   1,0,0, 1.3, Ra,  0,1,0, 1.3, Rb ) )/dx
                                          
    f_Xy.write("%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f \n" % (x, ss_pp, spx_pp, spy_pp, pxpx_pp, pxpy_pp ) )

    # d/dXz :
    ss_pp   = (pseudopot02(1.0, 1.0, 2.0,Rc + drz,   0,0,0, 1.3, Ra,  0,0,0, 1.3, Rb ) - pseudopot02(1.0, 1.0, 2.0,Rc - drz,   0,0,0, 1.3, Ra,  0,0,0, 1.3, Rb ) )/dx
    spx_pp  = (pseudopot02(1.0, 1.0, 2.0,Rc + drz,   0,0,0, 1.3, Ra,  1,0,0, 1.3, Rb ) - pseudopot02(1.0, 1.0, 2.0,Rc - drz,   0,0,0, 1.3, Ra,  1,0,0, 1.3, Rb ) )/dx
    spy_pp  = (pseudopot02(1.0, 1.0, 2.0,Rc + drz,   0,0,0, 1.3, Ra,  0,1,0, 1.3, Rb ) - pseudopot02(1.0, 1.0, 2.0,Rc - drz,   0,0,0, 1.3, Ra,  0,1,0, 1.3, Rb ) )/dx
    pxpx_pp = (pseudopot02(1.0, 1.0, 2.0,Rc + drz,   1,0,0, 1.3, Ra,  1,0,0, 1.3, Rb ) - pseudopot02(1.0, 1.0, 2.0,Rc - drz,   1,0,0, 1.3, Ra,  1,0,0, 1.3, Rb ) )/dx
    pxpy_pp = (pseudopot02(1.0, 1.0, 2.0,Rc + drz,   1,0,0, 1.3, Ra,  0,1,0, 1.3, Rb ) - pseudopot02(1.0, 1.0, 2.0,Rc - drz,   1,0,0, 1.3, Ra,  0,1,0, 1.3, Rb ) )/dx
                                          
    f_Xz.write("%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f \n" % (x, ss_pp, spx_pp, spy_pp, pxpx_pp, pxpy_pp ) )



    # d/dAx :
    ss_pp   = (pseudopot02(1.0, 1.0, 2.0,Rc,   0,0,0, 1.3, Ra + drx,  0,0,0, 1.3, Rb ) - pseudopot02(1.0, 1.0, 2.0,Rc,   0,0,0, 1.3, Ra - drx,  0,0,0, 1.3, Rb ) )/dx
    spx_pp  = (pseudopot02(1.0, 1.0, 2.0,Rc,   0,0,0, 1.3, Ra + drx,  1,0,0, 1.3, Rb ) - pseudopot02(1.0, 1.0, 2.0,Rc,   0,0,0, 1.3, Ra - drx,  1,0,0, 1.3, Rb ) )/dx
    spy_pp  = (pseudopot02(1.0, 1.0, 2.0,Rc,   0,0,0, 1.3, Ra + drx,  0,1,0, 1.3, Rb ) - pseudopot02(1.0, 1.0, 2.0,Rc,   0,0,0, 1.3, Ra - drx,  0,1,0, 1.3, Rb ) )/dx
    pxpx_pp = (pseudopot02(1.0, 1.0, 2.0,Rc,   1,0,0, 1.3, Ra + drx,  1,0,0, 1.3, Rb ) - pseudopot02(1.0, 1.0, 2.0,Rc,   1,0,0, 1.3, Ra - drx,  1,0,0, 1.3, Rb ) )/dx
    pxpy_pp = (pseudopot02(1.0, 1.0, 2.0,Rc,   1,0,0, 1.3, Ra + drx,  0,1,0, 1.3, Rb ) - pseudopot02(1.0, 1.0, 2.0,Rc,   1,0,0, 1.3, Ra - drx,  0,1,0, 1.3, Rb ) )/dx
                                          
    f_Ax.write("%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f \n" % (x, ss_pp, spx_pp, spy_pp, pxpx_pp, pxpy_pp ) )

    # d/dAy :
    ss_pp   = (pseudopot02(1.0, 1.0, 2.0,Rc,   0,0,0, 1.3, Ra + dry,  0,0,0, 1.3, Rb ) - pseudopot02(1.0, 1.0, 2.0,Rc,   0,0,0, 1.3, Ra - dry,  0,0,0, 1.3, Rb ) )/dx
    spx_pp  = (pseudopot02(1.0, 1.0, 2.0,Rc,   0,0,0, 1.3, Ra + dry,  1,0,0, 1.3, Rb ) - pseudopot02(1.0, 1.0, 2.0,Rc,   0,0,0, 1.3, Ra - dry,  1,0,0, 1.3, Rb ) )/dx
    spy_pp  = (pseudopot02(1.0, 1.0, 2.0,Rc,   0,0,0, 1.3, Ra + dry,  0,1,0, 1.3, Rb ) - pseudopot02(1.0, 1.0, 2.0,Rc,   0,0,0, 1.3, Ra - dry,  0,1,0, 1.3, Rb ) )/dx
    pxpx_pp = (pseudopot02(1.0, 1.0, 2.0,Rc,   1,0,0, 1.3, Ra + dry,  1,0,0, 1.3, Rb ) - pseudopot02(1.0, 1.0, 2.0,Rc,   1,0,0, 1.3, Ra - dry,  1,0,0, 1.3, Rb ) )/dx
    pxpy_pp = (pseudopot02(1.0, 1.0, 2.0,Rc,   1,0,0, 1.3, Ra + dry,  0,1,0, 1.3, Rb ) - pseudopot02(1.0, 1.0, 2.0,Rc,   1,0,0, 1.3, Ra - dry,  0,1,0, 1.3, Rb ) )/dx
                                          
    f_Ay.write("%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f \n" % (x, ss_pp, spx_pp, spy_pp, pxpx_pp, pxpy_pp ) )

    # d/dAz :
    ss_pp   = (pseudopot02(1.0, 1.0, 2.0,Rc,   0,0,0, 1.3, Ra + drz,  0,0,0, 1.3, Rb ) - pseudopot02(1.0, 1.0, 2.0,Rc,   0,0,0, 1.3, Ra - drz,  0,0,0, 1.3, Rb ) )/dx
    spx_pp  = (pseudopot02(1.0, 1.0, 2.0,Rc,   0,0,0, 1.3, Ra + drz,  1,0,0, 1.3, Rb ) - pseudopot02(1.0, 1.0, 2.0,Rc,   0,0,0, 1.3, Ra - drz,  1,0,0, 1.3, Rb ) )/dx
    spy_pp  = (pseudopot02(1.0, 1.0, 2.0,Rc,   0,0,0, 1.3, Ra + drz,  0,1,0, 1.3, Rb ) - pseudopot02(1.0, 1.0, 2.0,Rc,   0,0,0, 1.3, Ra - drz,  0,1,0, 1.3, Rb ) )/dx
    pxpx_pp = (pseudopot02(1.0, 1.0, 2.0,Rc,   1,0,0, 1.3, Ra + drz,  1,0,0, 1.3, Rb ) - pseudopot02(1.0, 1.0, 2.0,Rc,   1,0,0, 1.3, Ra - drz,  1,0,0, 1.3, Rb ) )/dx
    pxpy_pp = (pseudopot02(1.0, 1.0, 2.0,Rc,   1,0,0, 1.3, Ra + drz,  0,1,0, 1.3, Rb ) - pseudopot02(1.0, 1.0, 2.0,Rc,   1,0,0, 1.3, Ra - drz,  0,1,0, 1.3, Rb ) )/dx
                                          
    f_Az.write("%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f \n" % (x, ss_pp, spx_pp, spy_pp, pxpx_pp, pxpy_pp ) )




f_Xx.close()
f_Xy.close()
f_Xz.close()
f_Ax.close()
f_Ay.close()
f_Az.close()





