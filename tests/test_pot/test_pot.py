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



print "\nTest 2: SWITCH"
r1 = VECTOR(0.0, 0.0, 0.0)
r2 = VECTOR(1.0, 0.0, 0.0)

f = open("switch.txt","w")
for i in range(0,150):
    x = 0.1 * i
    r2.x = x   
    SW, dSW = 0.0, VECTOR(0.0, 0.0, 0.0)
    SW, dSW = SWITCH(r1, r2, 5.0, 10.0 )
    f.write("%8.5f  %8.5f  %8.5f\n" % (x, SW, dSW.x))
    #print x, SW, dSW.x
f.close()


print "\nTest 3: DOUBLE_SWITCH"
f = open("double_switch.txt","w")
for i in range(0,150):
    x = 0.1 * i
    SW, dSW = DOUBLE_SWITCH(x, 8.0, 1.0)
    f.write("%8.5f  %8.5f  %8.5f\n" % (x, SW, dSW))
    #print x, SW, dSW.x
f.close()



print "\nTest 4: Bond_Harmonic"
r1 = VECTOR()
r2 = VECTOR()

f = open("bond_harmonic.txt","w")
for i in range(-10,10):
    x = 0.1 * i
    r2.x = 1.5 + x   
    en, f1, f2 = 0.0, VECTOR(), VECTOR()
    en, f1, f2 = Bond_Harmonic(r1, r2, 100.0, 1.5 )
    f.write("%8.5f  %8.5f  %8.5f  %8.5f\n" % (x, en, f1.x, f2.x))
    #print x, SW, dSW.x
f.close()


print "\nTest 5: Bond_Quartic"
r1 = VECTOR()
r2 = VECTOR()

f = open("bond_quartic.txt","w")
for i in range(-10,20):
    x = 0.1 * i
    r2.x = 1.5 + x   
    en, f1, f2 = 0.0, VECTOR(), VECTOR()
    en, f1, f2 = Bond_Quartic(r1, r2, 100.0, 1.5 )
    f.write("%8.5f  %8.5f  %8.5f  %8.5f\n" % (x, en, f1.x, f2.x))
    #print x, SW, dSW.x
f.close()


print "\nTest 6: Bond_Morse"
r1 = VECTOR()
r2 = VECTOR()

f = open("bond_morse.txt","w")
for i in range(-10,100):
    x = 0.1 * i
    r2.x = 1.5 + x   
    en, f1, f2 = 0.0, VECTOR(), VECTOR()
    en, f1, f2 = Bond_Morse(r1, r2, 100.0, 1.5, 1.0)
    f.write("%8.5f  %8.5f  %8.5f  %8.5f\n" % (x, en, f1.x, f2.x))
    #print x, SW, dSW.x
f.close()



print "\nTest 7: Angle_Harmonic"
r1 = VECTOR(1.0, 0.0, 0.0)
r2 = VECTOR(0.0, 0.0, 0.0)
r3 = VECTOR(1.0, 0.0, 0.0)

f = open("angle_harmonic.txt","w")
f.write("theta(rad)   energy     f1.x    f1.y    f1.z   f2.x    f2.y    f2.z   f3.x    f3.y    f3.z\n")
for i in range(-2,180):
    theta = math.radians(i)
    r3.x = 1.0*math.cos(theta)
    r3.y = 1.0*math.sin(theta)

    en, f1, f2, f3 = 0.0, VECTOR(), VECTOR(), VECTOR()
    en, f1, f2, f3 = Angle_Harmonic(r1, r2, r3, 10.0, math.radians(90.0))
    f.write("%8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f\n" % (theta, en, f1.x, f1.y, f1.z, f2.x, f2.y, f2.z, f3.x, f3.y, f3.z))
    #print x, SW, dSW.x
f.close()



print "\nTest 8: Angle_Harmonic_Cos"
r1 = VECTOR(1.0, 0.0, 0.0)
r2 = VECTOR(0.0, 0.0, 0.0)
r3 = VECTOR(1.0, 0.0, 0.0)

f = open("angle_harmonic_cos.txt","w")
f.write("theta(rad)   energy     f1.x    f1.y    f1.z   f2.x    f2.y    f2.z   f3.x    f3.y    f3.z\n")
for i in range(-2,180):
    theta = math.radians(i)
    r3.x = 1.0*math.cos(theta)
    r3.y = 1.0*math.sin(theta)

    en, f1, f2, f3 = 0.0, VECTOR(), VECTOR(), VECTOR()
    en, f1, f2, f3 = Angle_Harmonic_Cos(r1, r2, r3, 10.0, math.cos(math.radians(90.0)), 2)
    f.write("%8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f\n" % (theta, en, f1.x, f1.y, f1.z, f2.x, f2.y, f2.z, f3.x, f3.y, f3.z))
    #print x, SW, dSW.x
f.close()





