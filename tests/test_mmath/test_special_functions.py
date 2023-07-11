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


print "\nTest 2: Testing FAST_POW"
for x in [[1.0, 1], [1.0, 0], [0.0, 0], [2.0, 0], [2.0, -2], [0.5, 2], [0.1,10]]:
    a = x[0]
    n = x[1]
    print "%8.3f^%3d = %8.5e" % (a, n, FAST_POW(a, n))

print "\nTest 3: Testing sinh_(x) = sinh(x)/x"
for i in range(-10,10):
    print "sinh_(%8.5f) = %8.5e" % (i*0.25, sinh_(i*0.25))

print "\nTest 4: Testing sin_(x) = sin(x)/x"
for i in range(-10,10):
    print "sin_(%8.5f) = %8.5e" % (i*0.25, sin_(i*0.25))

print "\nTest 5: Testing ERF and ERFC"
for i in range(0,20):
    x = i*0.25
    ef = ERF(x)
    efc = ERFC(x)
    su = ef + efc
    print "ERF(%8.5f) = %8.5e  ERFC(%8.5f) = %8.5e  sum = %8.5e" % (x, ef, x, efc, su)

