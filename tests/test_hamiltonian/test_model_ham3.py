#*********************************************************************************
#* Copyright (C) 2017 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/

#####################################################################################
#
#   This example shows how to perform calculations of the eigenvalues and eigenvectors
#   of a model Hamiltonian: setting up the Hamiltonian object, its parameters, doing
#   the calculations, accessing results, verifying that the calculations are as
#   we expect
#
#####################################################################################

import os
import sys
import math

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *

kb = 3.166811429e-6


# Dephasing by the time-dependent perturbation
ham_indx = 7 
T = 300.0
E0 = -1.0 * kb * T;  E1 =  1.0 * kb * T;  V01=  1.0 * kb * T  
D  =  0.0; L = 5.0; a = 0.1


# First, create the Hamiltonian
ham = Hamiltonian_Model(ham_indx)  # periodic SIN potential
ham.set_rep(1)                     # adiabatic
ham.set_params([E0, E1, V01, D, L ])
ham.set_v([ 0.0 ])
ham.set_q([ 0.0 ])
ham.compute()

C = ham.get_basis_transform()
H_dia = ham.get_ham_dia()
H_adi = ham.get_ham_adi()

print "H_dia = "; H_dia.show_matrix()
print "H_adi = "; H_adi.show_matrix()
print "C = "; C.show_matrix()
print "C.T() * H_dia * C = "; (C.T() * H_dia * C).show_matrix()  # <psi|H|psi> - adiabatic
print "C * H_adi * C.T() = "; (C * H_adi * C.T()).show_matrix()  # "inverse"


print "========================"

# Next, lets perturb the Hamiltonian by modifying the coupling (in the diabatic basis)
# See how the gap between adiabatic energies changes.

ham2 = Hamiltonian_Model(ham_indx)  # periodic SIN potential
ham2.set_rep(1)  # adiabatic
ham2.set_params([E0, E1, V01+0.0001, D, L ])
ham2.set_v([ 0.0 ])
ham2.set_q([ 0.0 ])
ham2.compute()

C2 = ham2.get_basis_transform()
H_dia2 = ham2.get_ham_dia()
H_adi2 = ham2.get_ham_adi()

print "H_dia2 = "; H_dia2.show_matrix()
print "H_adi2 = "; H_adi2.show_matrix()
print "C2 = "; C2.show_matrix()
print "C2.T() * H_dia2 * C2 = "; (C2.T() * H_dia2 * C2).show_matrix()  # <psi_tile|H_til|psi_til> - adiabatic
print "C2 * H_adi2 * C2.T() = "; (C2 * H_adi2 * C2.T()).show_matrix()  # "inverse"


