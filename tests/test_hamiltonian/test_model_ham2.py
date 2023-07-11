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


print "\nTest 2: Object of Hamiltonian_Model class"
ham = Hamiltonian_Model(1)

print "\nTest 3: set_params()"
ham.set_rep(1)  # adiabatic
#ham.set_params([0.1, 0.1, 1.6, 1.0])

print "\nTest 3.1: set_q()"
ham.set_q([1.0])  

print "\nTest 3.2: set_v()"
ham.set_v([0.0])  


print "\nTest 4: compute_diabatic()"
ham.compute_diabatic()

print "\nTest 5: compute_adiabatic()"
ham.compute_adiabatic()

print "\nTest 6: DAC model: diabatic, adiabatic, etc."
f = open("dac_adia.txt","w")
f1 = open("dac_adia1.txt","w")
f2 = open("dac_adia2.txt","w")
f3 = open("dac_adia3.txt","w")

fd = open("dac_dia.txt","w")
fd1 = open("dac_dia1.txt","w")
fd2 = open("dac_dia2.txt","w")
fd3 = open("dac_dia3.txt","w")

for i in range(-100,100):
    x = 0.1 * i    
    ham.set_rep(1)
    ham.set_q([x])
    ham.set_v([1.0])
    ham.compute()

    # Adiabatic PES and couplings
    f.write("%8.5f  %8.5f  %8.5f  %8.5f\n" % (x, ham.H(0,0).real, ham.H(1,1).real, ham.H(0,1).real ) )
    f1.write("%8.5f  %8.5f  %8.5f  %8.5f\n" % (x, ham.dHdq(0,0, 0).real, ham.dHdq(1,1, 0).real, ham.dHdq(0,1, 0).real ) )
    f2.write("%8.5f  %8.5f  %8.5f  %8.5f\n" % (x, ham.D(0,0, 0).real, ham.D(1,1, 0).real, ham.D(0,1, 0).real ) )
    f3.write("%8.5f  %8.5f  %8.5f  %8.5f\n" % (x, ham.Hvib(0,0).real, ham.Hvib(1,1).real, ham.Hvib(0,1).imag ) )

    # Now print all similar results for diabatic representations
    ham.set_rep(0) 
    fd.write("%8.5f  %8.5f  %8.5f  %8.5f\n" % (x, ham.H(0,0).real, ham.H(1,1).real, ham.H(0,1).real ) )
    fd1.write("%8.5f  %8.5f  %8.5f  %8.5f\n" % (x, ham.dHdq(0,0, 0).real, ham.dHdq(1,1, 0).real, ham.dHdq(0,1, 0).real ) )
    fd2.write("%8.5f  %8.5f  %8.5f  %8.5f\n" % (x, ham.D(0,0, 0).real, ham.D(1,1, 0).real, ham.D(0,1, 0).real ) )
    fd3.write("%8.5f  %8.5f  %8.5f  %8.5f\n" % (x, ham.Hvib(0,0).real, ham.Hvib(1,1).real, ham.Hvib(0,1).real ) )

#    H2 = res[3]
#    f2.write("%8.5f  %8.5f  %8.5f  %8.5f\n" % (res[0], H2.get(0,0), H2.get(1,1), H2.get(0,1) ) )

f.close()
f1.close()
f2.close()
f3.close()

fd.close()
fd1.close()
fd2.close()
fd3.close()







