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
print "\nTest 2: Gaussian primitive objects"
f = open("G_prim.txt","w")
Ra = VECTOR(0.0, 0.0, 0.0)
Rb = VECTOR(0.0, 0.0, 0.0)
u = VECTOR(1.0, 2.0, -0.5)
u.normalize()
print "u (normalized) = ", u.x, u.y, u.z

ga = PrimitiveG()
gb = PrimitiveG()


for i in range(0,50):
    x = 0.1 * i
    Rb = Ra + x*u    

    ss  = gaussian_overlap(0,0,0, 1.3, Ra,  0,0,0, 1.3, Rb )  # s(H)-s(H)
    pxpx = gaussian_overlap(1,0,0, 1.3, Ra,  1,0,0, 1.3, Rb )  # px(H)-px(H) 
    dz2dz2 = gaussian_overlap(0,0,2, 1.3, Ra,  0,0,2, 1.3, Rb )  # dz2(H)-dz2(H) 

    
    # Now do the same with Gaussinas:
    ga.init(0,0,0, 1.3, Ra)
    gb.init(0,0,0, 1.3, Rb)
    g_ss = gaussian_overlap(ga,gb)

    ga.x_exp, ga.y_exp, ga.z_exp = 1, 0, 0 
    gb.x_exp, gb.y_exp, gb.z_exp = 1, 0, 0 
    g_pxpx = gaussian_overlap(ga,gb)

    ga.x_exp, ga.y_exp, ga.z_exp = 0, 0, 2 
    gb.x_exp, gb.y_exp, gb.z_exp = 0, 0, 2
    g_dz2dz2 = gaussian_overlap(ga,gb)

    
    f.write("%8.5f   %8.5f %8.5f %8.5f   %8.5f %8.5f %8.5f \n" % (x, ss, pxpx, dz2dz2, g_ss, g_pxpx, g_dz2dz2) )

f.close()

t.stop()
print "Time to compute = ", t.show(), " sec"


print "\nTest 3: Move Gaussian by reference"
ga = PrimitiveG(0,0,0, 1.3, Ra)
gb = PrimitiveG(0,0,0, 1.3, Rb)

f = open("G_prim1.txt","w")

for i in range(0,50):
    x = 0.1 * i
    Rb = Ra + x*u    
    ss  = gaussian_overlap(0,0,0, 1.3, Ra,  0,0,0, 1.3, Rb )  # s(H)-s(H)

    gb.R = Rb
    g_ss = gaussian_overlap(ga,gb)  
    
    f.write("%8.5f   %8.5f %8.5f \n" % (x, ss, g_ss) )


f.close()



