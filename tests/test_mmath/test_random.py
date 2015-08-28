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


print "\nTest 1: Importing the library and its content"
print "from cygmmath import *"
from cygmmath import *

print "\nTest 2: Constructor & uniform"
r = Random()
for i in range(0,10):
    print i, r.uniform(0.0, 1.0)

print "\nTest 3: exponential"
for i in range(0,10):
    print i, r.exponential(0.01)

print "\nTest 4: normal"
for i in range(0,10):
    print i, r.normal()

print "\nTest 5: gamma"
for i in range(0,10):
    print i, r.gamma(1.0)

print "\nTest 6: beta"
for i in range(0,10):
    print i, r.beta(1.0, 0.5)

print "\nTest 7: poiss1"
for i in range(0,10):
    print i, r.poiss1(2.0)





