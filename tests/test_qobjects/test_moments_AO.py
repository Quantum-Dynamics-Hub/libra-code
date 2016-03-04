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

import os
import sys
import math

# Fisrt, we add the location of the library to test to the PYTHON path
cwd = os.getcwd()
print "Current working directory", cwd
sys.path.insert(1,cwd+"/../../_build/src/mmath")
sys.path.insert(1,cwd+"/../../_build/src/qchem")


print "\nImporting the library and its content"
from cygmmath import *
from cygqchem import *


print "\nConstruct AO objects: s-type primitives only"

# Pople values for C 2s orbitals
alp_2s = [0.994203, 0.231031,  0.0751386]
coeff_2s = [-0.09996723,  0.39951283,  0.70011547]

# Pople values for C 2p orbitals
alp_2p = [0.994203, 0.231031,  0.0751386]
coeff_2p = [0.15591627, 0.60768372, 0.39195739]


ksi = 1.3

#======== s-orbitals ====================
aoa_2s_v1 = AO()
aoa_2s_v2 = AO()
for i in range(0,3):
    aoa_2s_v1.add_primitive(coeff_2s[i],PrimitiveG(0,0,0, ksi*ksi*alp_2s[i], VECTOR(0,0,0)))
    aoa_2s_v2.add_primitive(coeff_2s[i],PrimitiveG(0,0,0, (2.0*ksi)*(2.0*ksi)*alp_2s[i], VECTOR(0,0,0)))
aoa_2s_v1.normalize()
aoa_2s_v2.normalize()

#========== p-orbitals ==================
aoa_2p_v1 = AO()
aoa_2p_v2 = AO()
for i in range(0,3):
    aoa_2p_v1.add_primitive(coeff_2p[i],PrimitiveG(1,0,0, ksi*ksi*alp_2p[i], VECTOR(0,0,0)))
    aoa_2p_v2.add_primitive(coeff_2p[i],PrimitiveG(1,0,0, (2.0*ksi)*(2.0*ksi)*alp_2p[i], VECTOR(0,0,0)))
aoa_2p_v1.normalize()
aoa_2p_v2.normalize()



print "\nTest STO-3G moments"
v = VECTOR(0.0,0.0,0.0)
gx = PrimitiveG(1,0,0, 0.00001,v)

f = open("moments.txt","w")
for i in range(0,50):
    x = 0.1*i

    aoa_2s_v2.set_position(VECTOR(x,0.0,0.0))
    aoa_2p_v2.set_position(VECTOR(x,0.0,0.0))

    s_s = gaussian_overlap(aoa_2s_v1, aoa_2s_v2)
    s_p = gaussian_overlap(aoa_2s_v1, aoa_2p_v2)

    s_x_s = gaussian_moment(aoa_2s_v1, gx, aoa_2s_v2)
    s_x_p = gaussian_moment(aoa_2s_v1, gx, aoa_2p_v2)

# Using explicit 3D moment on a couple of primitives
#    nx1, ny1, nz1 = aoa_2s_v1.primitives[0].x_exp, aoa_2s_v1.primitives[0].y_exp, aoa_2s_v1.primitives[0].z_exp
#    gx1, gy1, gz1 = gx.x_exp, gx.y_exp, gx.z_exp
#    nx2, ny2, nz2 = aoa_2s_v2.primitives[0].x_exp, aoa_2s_v2.primitives[0].y_exp, aoa_2s_v2.primitives[0].z_exp
#    s_x_s = gaussian_moment(nx1,ny1,nz1, aoa_2s_v1.primitives[0].alpha, aoa_2s_v1.primitives[0].R,
#                            gx1,gy1,gz1, gx.alpha, gx.R,
#                            nx2,ny2,nz2, aoa_2s_v2.primitives[0].alpha, aoa_2s_v2.primitives[0].R
#                           )

    
    f.write("%8.5f   %8.5f  %8.5f  %8.5f  %8.5f\n" % (x, s_s, s_p, s_x_s, s_x_p ) )

f.close()


