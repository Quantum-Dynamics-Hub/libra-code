#*********************************************************************************
#* Copyright (C) 2019 Xiang Sun, Alexey V. Akimov
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

from libra_py import units
import fgr_py
fs = units.fs2au

tmax = 20.0
dt = 0.20
dtau = dt/50.0

gamma = 0.1
nomega = 500
w_c = 1.0
dw = 15.0*w_c/float(nomega)


for dyn_type in [0, 1]:  # Condon vs. non-Condon
    for method in [0, 1, 2, 3, 4, 5]:  #  All methods
        
        i1 = -1
        for w_DA in [0.0]: #[0.0, 2.0]:   # Donor-Acceptor energy gap
            i1 += 1
            i2 = -1
            for s in [-1.0]: # -1.0, 1.0, 3.0]:   # Noneq. initial shift of primary mode
                i2 += 1
                i3 = -1
                for beta in [1.0]: #, 2.0, 5.0]:   # themal energy
                    i3 += 1
                    i4 = -1
                    for etha in [0.5]: #, 1.0, 2.0]:  # friction
                        i4 += 1

                        print "-----------------------"
                        print "w_DA = ", w_DA
                        print "s = ", s
                        print "beta = ", beta
                        print "etha = ", etha
        
                        #============ Setup the parameters ===============
        
                        params = {}
                        params["Er"] = 0.5 * w_c
                        params["omega_DA"] = w_DA * w_c
                        params["omega"] = [0.5*w_c]
                        params["coup"] = [0.0]
        
                        for a in xrange(1, nomega):
                            w_a = a*dw
                            params["omega"].append(w_a)
        
                            J_a = etha * w_a * math.exp(-w_a/w_c)      # Eq. 45
                            c_a = math.sqrt((2.0/math.pi)*J_a*w_a*dw)  # Eq. 62
                            params["coup"].append(c_a)
        
                        #============ Setup the parameters =============== 
                        ndof = len(params["omega"])
                        print "ndof = ", ndof
                        omega = Py2Cpp_double(params["omega"])
                        coeff = Py2Cpp_double(params["coup"])
        
                        U = MATRIX(ndof, ndof)
                        omega_nm = normal_modes(omega, coeff, U)
                        #print "U = "; U.show_matrix()
        
                        dE = params["omega_DA"]
        
                        print "Omega = ", params["omega"][0]
                        print "Er = ", params["Er"]
                        y0 = eq_shift(params["Er"], params["omega"][0])
                        print "y0 = ", y0
        
                        req_nm = compute_req(omega, coeff, y0, U)
                        #print "req_nm = ", Cpp2Py(req_nm)
        
                        gamma_nm = compute_TT_scaled(U, gamma)
                        #print "gamma_nm = ", Cpp2Py(gamma_nm)
        
                        shift_NE = compute_TT_scaled(U, s)
                        #print "shift_NE = ", Cpp2Py(shift_NE)
        
                        V = coupling_Condon(gamma, dE, params["Er"], y0)
                        print "V = ", V
        
                        fgr_py.run_NEFGRL_populations(params["omega_DA"], V, omega_nm, gamma_nm,req_nm, shift_NE, method, beta, dyn_type, dtau, tmax, dt, "_res-%i-%i-%i-%i-%i-%i.txt" % (dyn_type, method, i1,i2,i3,i4) )
        


