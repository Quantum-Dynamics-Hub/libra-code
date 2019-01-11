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

import sys
import cmath
import math
import os

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from libra_py import models_LVC, units


def run_NEFGRL_populations(V, omega_DA, omega_nm, gamma_nm, req_nm, shift_NE, method, beta, dyn_type, dtau, tmax, dt, filename):

    
    f = open(filename, "w")
    f.close()

    nsteps = int(tmax/dt)+1
    summ, P, k = 0.0, 1.0, 0.0  # probability of donor state


    for step in xrange(nsteps):

        t = step*dt

        f = open(filename, "a")
        f.write("%8.5f  %8.5f  %8.5f \n" % (t, k, P))
        f.close()
 
        # k = k(t')
        k = NEFGRL_rate(t, V, omega_DA, omega_nm, gamma_nm, req_nm, shift_NE, method, beta, dyn_type, dtau)
        
        summ += k * dt; 
        P = math.exp(-summ)  # exp(- int dt' k(t'))

