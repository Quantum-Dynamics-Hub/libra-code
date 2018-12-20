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

fs = units.fs2au

params = models_LVC.get_LVC_set1()
 
T = 300.0 # K
beta = 1.0 / (units.kB * T)
tmax = 50 * fs
dt = 0.25 * fs
dtau = dt/5.0

gamma = 1.0
s = -1.0


method = 5


ndof = len(params["omega"])
print "ndof = ", ndof
omega = Py2Cpp_double(params["omega"])
coeff = Py2Cpp_double(params["coup"])
d1 = Py2Cpp_double(params["d1"])
d2 = Py2Cpp_double(params["d2"])


U = MATRIX(ndof, ndof)
omega_nm = normal_modes(omega, coeff, U)
print "U = "; U.show_matrix()


dE = LVC2GOA_dE(params["Delta1"], params["Delta2"], omega_nm, d1, d2)
print "dE = ", dE
print "Omega_DA = ", params["omega_DA"]

#req_nm = compute_req(omega, coeff, y0, U)
req_nm = LVC2GOA_req(omega_nm, d1, d2)
print "req_nm = ", Cpp2Py(req_nm)


print "Omega = ", params["omega"][0]
print "Er = ", params["Er"]
y0 = eq_shift(params["Er"], params["omega"][0])
print "y0 = ", y0


gamma_nm = compute_TT_scaled(U, gamma)
print "gamma_nm = ", Cpp2Py(gamma_nm)

shift_NE = compute_TT_scaled(U, s)
print "shift_NE = ", Cpp2Py(shift_NE)

V = coupling_Condon(gamma, dE, params["Er"], y0)
print "V = ", V

res = NEFGRL_population(V, params["omega_DA"], dtau, omega_nm, gamma_nm,req_nm, shift_NE, method, tmax, dt, beta)
res.show_matrix("res.txt")



