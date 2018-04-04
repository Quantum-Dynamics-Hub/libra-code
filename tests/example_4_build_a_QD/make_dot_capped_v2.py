#*********************************************************************************
#* Copyright (C) 2018 Alexey V. Akimov
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
import copy
#import capping
#import new_autoconnect

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from libra_py import *


Angst = 1.889725989

# Unit cell parameters
a = VECTOR(5.43, 0.0, 0.0) * Angst
b = VECTOR(0.0, 5.43, 0.0) * Angst
c = VECTOR(0.0, 0.0, 5.43) * Angst

# Create a Chemical system
Rdot = 3.00 * Angst
Nx, Ny, Nz = 3, 3, 3  # 

Rbond = 3.001 * Angst

L0, R0 = build.read_xyz("Si.xyz")
cap_dict = {"Si":"H"}
masses = [1.0]*5
MaxCoord = [4,1,1,1,1]

L_bulk, R_bulk = build.generate_replicas_xyz2(L0, R0, a, b, c, Nx, Ny, Nz)
MaxCoord_bulk = []
masses = []
PT_coord = { "Si": 4, "H":1 }
for i in L_bulk:
    MaxCoord_bulk.append(PT_coord[i])    # maximal coordination number refering to atoms in the bulk crystal
    masses.append(1.0)               # whatever mass
#res, line, pairs = autoconnect.autoconnect_pbc(R_bulk, MaxCoord_bulk, Rbond, Nx*a, Ny*b, Nz*c, "none", 0, 0)
res, line, pairs = new_autoconnect.autoconnect_pbc(R_bulk, L_bulk, MaxCoord_bulk, Rbond, Nx*a, Ny*b, Nz*c, "none", 0, 0, cap_dict)
print "\n Supercell is now connected \n"
#print line
#sys.exit(0)

L_sphere, R_sphere = build.crop_sphere_xyz2(L_bulk, R_bulk, Rdot)
MaxCoord_sphere = []
for i in L_sphere:
    MaxCoord_sphere.append(PT_coord[i]) # maximal coordination number refering to atoms in the sphere
#res, line1, pairs = autoconnect.autoconnect_pbc(R_sphere, MaxCoord_sphere, Rbond, Nx*a, Ny*b, Nz*c, "none", 0, 0)
res, line1, pairs = new_autoconnect.autoconnect_pbc(R_sphere, L_sphere, MaxCoord_sphere, Rbond, Nx*a, Ny*b, Nz*c, "none", 0, 0, cap_dict)
print "\n Sphere is cropped and connected \n"
print line1
#sys.exit(0)

######### Now, we will do the capping ########

L_capped, R_capped = capping.cap_system(L_sphere, R_sphere, res, MaxCoord_sphere, cap_dict)
MaxCoord_capped = []
for i in L_capped:
    MaxCoord_capped.append(PT_coord[i])          # maximal coordination number
res, line2, pairs = new_autoconnect.autoconnect_pbc(R_capped, L_capped, MaxCoord_capped, Rbond, Nx*a, Ny*b, Nz*c, "none", 0, 0, cap_dict)
print "\n Sphere is capped and connected \n"
print line2
#sys.exit(0)

syst = System()
build.add_atoms_to_system(syst, L_capped, R_capped, VECTOR(1.0, 1.0, 1.0), masses, "elements.dat")
syst.print_ent("Si-QD.ent") 
#res, line, pairs = new_autoconnect.autoconnect_pbc(R0, L0, MaxCoord_capped, Rbond, Nx*a, Ny*b, Nz*c, "none", 1, 0, cap_dict)
#print "Autodetermined connections: "
#print line

f = open("Si-QD.ent", "a")
f.write("END\n")
f.write(line1 + line2)
f.write("END\n")
f.close()

