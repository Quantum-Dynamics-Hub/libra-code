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



def STO_1s(ksi, r):
    res = math.sqrt(ksi**3 / math.pi) * math.exp(-ksi*r)
    return res

def STO_2s(ksi, r):
    res = math.sqrt(ksi**5 / (3.0*math.pi) ) * r * math.exp(-ksi*r)
    return res

def STO_2p(ksi, r):
    res = math.sqrt(ksi**5 / (math.pi) ) * r * math.exp(-ksi*r)
    return res





print "\nTest 2: Construct AO objects: s-type primitives only"

# Pople values for C 2s orbitals
alp_2s = [0.994203, 0.231031,  0.0751386]
coeff_2s = [-0.09996723,  0.39951283,  0.70011547]


ksi = 1.3

aoa_2s_v1 = AO()
aoa_2s_v2 = AO()

for i in range(0,3):
    aoa_2s_v1.add_primitive(coeff_2s[i],PrimitiveG(0,0,0, ksi*ksi*alp_2s[i], VECTOR(0,0,0)))
    aoa_2s_v2.add_primitive(coeff_2s[i],PrimitiveG(0,0,0, (2.0*ksi)*(2.0*ksi)*alp_2s[i], VECTOR(0,0,0)))

print "norm of AO v1: <AO|AO> ", aoa_2s_v1.norm2()
print "norm of AO v2: <AO|AO> ", aoa_2s_v2.norm2()
aoa_2s_v1.normalize()
aoa_2s_v2.normalize()
print "norm of AO v1: <AO|AO> ", aoa_2s_v1.norm2()
print "norm of AO v2: <AO|AO> ", aoa_2s_v2.norm2()


print "\nTest 3: Test STO-3G"

f = open("sto-3g_2s.txt","w")
for i in range(0,50):
    x = 0.1*i
    res1 = aoa_2s_v1.compute(VECTOR(x,0.0,0.0))  # unnormalized sum of normalized Gaussians
    res2 = aoa_2s_v2.compute(VECTOR(x,0.0,0.0))  # as is multiply ksi by 2
    
    f.write("%8.5f   %8.5f  %8.5f  %8.5f  %8.5f\n" % (x, res1, STO_2s(ksi, x), res2, STO_2s(2.0*ksi, x)) )

f.close()



print "\nTest 4: Construct AO objects: p-type primitives only"

# Pople values for C 2p orbitals
alp_2p = [0.994203, 0.231031,  0.0751386]
coeff_2p = [0.15591627, 0.60768372, 0.39195739]


ksi = 1.3

aoa_2p_v1 = AO()
aoa_2p_v2 = AO()

for i in range(0,3):
    aoa_2p_v1.add_primitive(coeff_2p[i],PrimitiveG(1,0,0, ksi*ksi*alp_2p[i], VECTOR(0,0,0)))
    aoa_2p_v2.add_primitive(coeff_2p[i],PrimitiveG(1,0,0, (2.0*ksi)*(2.0*ksi)*alp_2p[i], VECTOR(0,0,0)))

print "norm of AO v1: <AO|AO> ", aoa_2p_v1.norm2()
print "norm of AO v2: <AO|AO> ", aoa_2p_v2.norm2()
aoa_2p_v1.normalize()
aoa_2p_v2.normalize()
print "norm of AO v1: <AO|AO> ", aoa_2p_v1.norm2()
print "norm of AO v2: <AO|AO> ", aoa_2p_v2.norm2()


print "\nTest 5: Test STO-3G"

f = open("sto-3g_2p.txt","w")
for i in range(0,50):
    x = 0.1*i
    res1 = aoa_2p_v1.compute(VECTOR(x,0.0,0.0))  # unnormalized sum of normalized Gaussians
    res2 = aoa_2p_v2.compute(VECTOR(x,0.0,0.0))  # as is multiply ksi by 2
    
    f.write("%8.5f   %8.5f  %8.5f  %8.5f  %8.5f\n" % (x, res1, STO_2p(ksi, x), res2, STO_2p(2.0*ksi, x)) )

f.close()




