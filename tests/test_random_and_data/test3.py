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
import math

# Fisrt, we add the location of the library to test to the PYTHON path
import sys

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from libra_py import *


def mapping(r, x0, N):
    res = []
    x = x0
    for i in xrange(N):
        res.append(x)
        x = 4.0*r*x*(1.0 - x)

    return res



y1 = mapping(0.98, 0.2, 10000)
dy1 = DATA(y1)

f = open("mapping.txt", "w")
for i in xrange(len(y1)):
    f.write("%8.5f  %8.5f \n" % (i, y1[i]))
f.close()

x = []
for i in range(0,100):
    x.append(i*0.01)
dens, cum = dy1.Calculate_Distribution(x)


f = open("case1.txt","w")
i = 0
sz = len(x)
for i in range(0, sz):
    f.write("%8.5f  %8.5f  %8.5f  \n" % (x[i], dens[i], cum[i] ) )
f.close()



