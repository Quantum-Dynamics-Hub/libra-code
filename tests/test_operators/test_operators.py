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



print "\nTest 2: Rotate"
r = [1.0, 0.0]
dphi = math.radians(10.0)
for i in range(0,10):
    r = rotate(r[0],r[1], dphi)
    phi = math.degrees(math.atan2(r[1],r[0]))
    print r, phi

print "\nTest 3: Shift"
x = 0.0
for i in range(0, 10):
    print i, x
    x = shift(x, 0.1)


print "\nTest 4: Scale"
x = 1.0
for i in range(0, 10):
    print i, x
    x = scale(x, -0.99)

    
#void shift(double&, double);
#double expt_shift(double x,double phi);

#void scale(double&, double);
#double expt_scale(double x,double phi);

