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
import fgr_py

fs = units.fs2au

params = models_LVC.get_LVC_set1b()  # parameters are the same as in Xiang's code
 
T = 100.0 # K
beta = 1.0 / (units.kB * T)
tmax = 10 * fs
dt = 0.1 * fs
dtau = dt/5.0
dyn_type = 1 # 0 - Condon, 1 - non-Condon

gamma = 1.0   # DA_coupling
s = -1.0


method = 0


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
dE = params["omega_DA"]

#req_nm = compute_req(omega, coeff, y0, U)
req_nm = LVC2GOA_req(omega_nm, d1, d2)
print "req_nm = ", Cpp2Py(req_nm)


print "Omega = ", params["omega"][0]
print "Er(from params) = ", params["Er"]
print "Er(computed) = ", reorganization_energy(omega_nm, req_nm)
y0 = eq_shift(params["Er"], params["omega"][0])
print "y0 = ", y0


#gamma_nm = compute_TT_scaled(U, gamma)
gamma_nm = Py2Cpp_double(params["coup"])
print "gamma_nm = ", Cpp2Py(gamma_nm)

#shift_NE = compute_TT_scaled(U, s)
shift_NE = doubleList()
sz = len(req_nm)
for i in xrange(sz):
    shift_NE.append(s * req_nm[i])

print "shift_NE = ", Cpp2Py(shift_NE)


#V = coupling_Condon(gamma, dE, params["Er"], y0)
#print "V = ", V
V = gamma

#res = NEFGRL_population(V, dE, omega_nm, gamma_nm,req_nm, shift_NE, method, beta, dyn_type, dtau, tmax, dt)
#res.show_matrix("_res.txt")

fgr_py.run_NEFGRL_populations(dE, V, omega_nm, gamma_nm,req_nm, shift_NE, method, beta, dyn_type, dtau, tmax, dt, "_res.txt")

"""
for step in xrange(50):
    t = step*dt
    k = NEFGRL_rate(t, V, dE, omega_nm, gamma_nm, req_nm, shift_NE, method, beta, dyn_type, dtau)
    print "Time [a.u.] = ", t, " rate constant[a.u.^-1] = ", k

    nomega = len(omega_nm)
    for w in range(nomega):
        tau = t
        integ = Integrand_NE_exact(t, tau, dE, omega_nm[w], shift_NE[w], req_nm[w], beta)
        lin = Linear_NE_exact(t, tau, gamma_nm[w], omega_nm[w], shift_NE[w], req_nm[w], beta)
        print w, omega_nm[w], integ, lin

"""




