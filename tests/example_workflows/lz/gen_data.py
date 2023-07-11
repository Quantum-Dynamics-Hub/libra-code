#*********************************************************************************
#* Copyright (C) 2018  Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************
#
# This script is to generate some data for QSH testing
#

import os
import sys
import math

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *
from libra_py.units import inv_cm2Ha


def set1():

    prefix = "res/Hvib_"

    w1 = 1000.0 * inv_cm2Ha  # in cm^-1
    w2 = 2500.0 * inv_cm2Ha  # in cm^-1
    w3 = 3000.0 * inv_cm2Ha  # in cm^-1
    dt = 41.0

    for i in xrange(5000):
        H = CMATRIX(2,2)
        t = dt*i

        E1 = 0.0 # 0.005*math.sin(w1*t)
        E2 = 0.0025 + 0.001*math.cos(w2*t)
        d = 0.01*math.sin(w3*t)

        H.set(0, 0, E1*(1.0+0.0j));  H.set(0, 1, d*(0.0-1.0j))
        H.set(1, 0, d*(0.0-1.0j));   H.set(1, 1, E2*(1.0+0.0j))

        H.real().show_matrix(prefix+"%i_re" % (i))
        H.imag().show_matrix(prefix+"%i_im" % (i))
    
set1()
