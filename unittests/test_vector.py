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

cwd = os.getcwd()
print "Current working directory", cwd
sys.path.insert(1,cwd+"/../_build/src/math_linalg")

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    #from cyglibra_core import *
    from cyglinalg import *

elif sys.platform=="linux" or sys.platform=="linux2":
    #from liblibra_core import *
    from liblinalg import *


print "\nTest 2: Test default constructor"
print "v1 = VECTOR()"
v1 = VECTOR()  
print "v1 = ", v1, v1.x, v1.y, v1.z

print "\nTest 3: Test 3-argument constructor"
print "v2 = VECTOR(1.0, -2.0, 3.0)"
v2 = VECTOR(1.0, -2.0, 3.0)
print "v2 = ", v2, v2.x, v2.y, v2.z

print "\nTest 4: Test yet another constructor"
print "v3 = VECTOR(v2)"
v3 = VECTOR(v2)
print "v3 = ", v3, v3.x, v3.y, v3.z


print "\nTest 5: Test assignment"
print "v3 = VECTOR()"
print "v3 = v1"
v3 = VECTOR()
v3 = v1
print "v2 = ", v2, v2.x, v2.y, v2.z
print "So v2 stays the same, even though we \"modified\" v3"
print "v3 = ", v3, v3.x, v3.y, v3.z
print "v1 = ", v1, v1.x, v1.y, v1.z
print "change v1.."
print "v1.x, v1.y, v1.z = 0.0, -1.0, 2.0"
v1.x, v1.y, v1.z = 0.0, -1.0, 2.0
print "v3 = ", v3, v3.x, v3.y, v3.z
print "v1 = ", v1, v1.x, v1.y, v1.z
print "change v3.."
print "v3.x, v3.y, v3.z = 4.0, -5.0, 3.0"
v3.x, v3.y, v3.z = 4.0, -5.0, 3.0
print "v3 = ", v3, v3.x, v3.y, v3.z
print "v1 = ", v1, v1.x, v1.y, v1.z
print "But the vector assignment is by reference"


print "\nTest 6: Arithmetics:"
print "Reset vectors:"
print "v1 = VECTOR()"
print "v2 = VECTOR()"
print "v3 = VECTOR()"
v1 = VECTOR()  
v2 = VECTOR()  
v3 = VECTOR()  
print v1
print v2
print v3
print "v1.x, v1.y, v1.z = 1.0, 2.0, 1.0"
print "v2.x, v2.y, v2.z = -1.0, 3.0, 0.0"
print "v3.x, v3.y, v3.z = 0.0, 0.0, 0.0"
v1.x, v1.y, v1.z =  1.0, 2.0, 1.0
v2.x, v2.y, v2.z = -1.0, 3.0, 0.0
v3.x, v3.y, v3.z =  0.0, 0.0, 0.0

print "v3 = v1 + v2"
v3 = v1 + v2
print "v3 = ", v3, v3.x, v3.y, v3.z

print "v3 = v1 - v2"
v3 = v1 - v2
print "v3 = ", v3, v3.x, v3.y, v3.z

print "v3 = 2.0*v1 - v2*3.0"
v3 = 2.0*v1 - v2*3.0
print "v3 = ", v3, v3.x, v3.y, v3.z

print "v3 = v1/2.0"
v3 = v1/2.0
print "v3 = ", v3, v3.x, v3.y, v3.z

print "v3 = v1 + 1.0"
v3 = v1 + 1.0
print "v3 = ", v3, v3.x, v3.y, v3.z


print "\nTest 7: Member functions"
print "v1 = ",v1,  v1.x, v1.y, v1.z
print "v2 = ",v2,  v2.x, v2.y, v2.z
print "v3 = ",v3,  v3.x, v3.y, v3.z

print "v1.length() = ",v1.length()
print "math.sqrt(6.0) = ",math.sqrt(6.0)
print "v2.length2() = ",v2.length2()

print "t = v1.unit()"
t = v1.unit()
print "t = ", t, t.x,t.y,t.z
print "v1 = ", v1, v1.x,v1.y,v1.z
print "v1.normalize()"
v1.normalize()
print "t = ", t, t.x,t.y,t.z
print "v1 = ", v1, v1.x,v1.y,v1.z


print "v1.x, v1.y, v1.z = 1.0, 0.0, 0.0"
print "v2.x, v2.y, v2.z = 0.0, 1.0, 0.0"
print "v3.x, v3.y, v3.z = 0.0, 0.0, 0.0"
v1.x, v1.y, v1.z =  1.0, 0.0, 0.0
v2.x, v2.y, v2.z =  0.0, 1.0, 0.0
v3.x, v3.y, v3.z =  0.0, 0.0, 0.0
print "v3.cross(v1,v2)"
v3.cross(v1,v2)
print "v3 =", v3, v3.x, v3.y, v3.z

print "v3.cross(v2,v1)"
v3.cross(v2,v1)
print "v3 =", v3, v3.x, v3.y, v3.z

print "v3.cross(v1,v1)"
v3.cross(v1,v1)
print "v3 =", v3, v3.x, v3.y, v3.z



print "\nTest 8: Vector indexing"
vlst = VECTORList()
vlst = [v1, v2, v3]
print vlst
print vlst[0], vlst[0].x, vlst[0].y, vlst[0].z
print vlst[1]
print vlst[2]

