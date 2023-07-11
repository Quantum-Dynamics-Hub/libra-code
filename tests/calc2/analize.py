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

######################################################################
#
# Analysis of azobenzene trajectory
#
######################################################################

import os
import sys
import math

# Fisrt, we add the location of the library to test to the PYTHON path
cwd = os.getcwd()
#cwd = "/projects/alexeyak/Software/libracode-code/tests/study2"

print "Current working directory", cwd

sys.path.insert(1,cwd+"/../../_build/src/mmath")
sys.path.insert(1,cwd+"/../../_build/src/dyn")
sys.path.insert(1,cwd+"/../../_build/src/chemobjects")
sys.path.insert(1,cwd+"/../../_build/src/hamiltonian")
sys.path.insert(1,cwd+"/../../_build/src/hamiltonian/Hamiltonian_Atomistic")
sys.path.insert(1,cwd+"/../../_build/src/hamiltonian/Hamiltonian_Atomistic/Hamiltonian_QM")
sys.path.insert(1,cwd+"/../../_build/src/hamiltonian/Hamiltonian_Atomistic/Hamiltonian_QM/Control_Parameters")
sys.path.insert(1,cwd+"/../../_build/src/hamiltonian/Hamiltonian_Atomistic/Hamiltonian_QM/Model_Parameters")
sys.path.insert(1,cwd+"/../../_build/src/hamiltonian/Hamiltonian_Atomistic/Hamiltonian_QM/Basis_Setups")
sys.path.insert(1,cwd+"/../../_build/src/dyn")
sys.path.insert(1,cwd+"/../../_build/src/qchem/qobjects")
sys.path.insert(1,cwd+"/../../_build/src/qchem/basis")
sys.path.insert(1,cwd+"/../../_build/src/converters")
sys.path.insert(1,cwd+"/../../_build/src/calculators")

print "\nTest 1: Importing the library and its content"
from cygmmath import *
from cygdyn import *
from cygchemobjects import *
from cyghamiltonian import *
from cyghamiltonian_qm import *
from cygcontrol_parameters import *
from cygmodel_parameters import *
from cygbasis_setups import *
from cygdyn import *
from cygqobjects import *
from cygbasis import *
from cygconverters import *
from cygcalculators import *



def angle(ri,rj,rk,rl):
    t = VECTOR()
    u = VECTOR()
    rp = VECTOR()

    direction = -1.0  # dihedral angle

     
    rij = ri - rj
    rkj = rk - rj
    rlk = rl - rk
    t.cross(rij, rkj)
    t = direction*t

    u.cross(rlk, rkj)
    modt = t.length()
    modu = u.length()

    phi = 0.0

    if((modt!=0.0) and (modu!=0.0)):
        cos_phi = (t*u)/(modt*modu)
        if(cos_phi>= 1.0):
            cos_phi=1.0
        elif(cos_phi<=-1.0):
            cos_phi=-1.0

        rp.cross(t,u)
        sgn = 0.0
        if rkj*rp <0.0:
            sgn = -1.0
        else:
            sgn = 1.0

        phi = sgn*math.acos(cos_phi)

    
    return phi



#C3 - N12 - N13 - C14
f = open("tsh_copy_0.xyz","r")
A = f.readlines()
f.close()

sz = len(A)
nsnaps = sz / (24 + 2)

print "number of snaps = ", nsnaps

r1 = VECTOR()
r2 = VECTOR()
r3 = VECTOR()
r4 = VECTOR()


for t in xrange(nsnaps):
    s_C1 = A[t*26+(3+2)].split()
    s_N2 = A[t*26+(12+2)].split()
    s_N3 = A[t*26+(13+2)].split()
    s_C4 = A[t*26+(14+2)].split()
 
    r1.x, r1.y, r1.z = float(s_C1[1]), float(s_C1[2]), float(s_C1[3])
    r2.x, r2.y, r2.z = float(s_N2[1]), float(s_N2[2]), float(s_N2[3])
    r3.x, r3.y, r3.z = float(s_N3[1]), float(s_N3[2]), float(s_N3[3])
    r4.x, r4.y, r4.z = float(s_C4[1]), float(s_C4[2]), float(s_C4[3])
    
    print t, angle(r1,r2,r3,r4)

    

