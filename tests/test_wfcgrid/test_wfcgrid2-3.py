#*********************************************************************************
#* Copyright (C) 2019 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
"""

 This example demonstrates how to use exact wfc calculation for 2-state Holstein

 Hamiltonian - in diabatic and adiabatic representations.

"""

import os
import sys
import math
import time
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *

import libra_py.models.Holstein as Holstein


def compute_model(q, params):

    full_id = Py2Cpp_int([0,0])   
    res = Holstein.Holstein2(q, params, full_id)

    return res


    
x0 = Py2Cpp_double([0.0])
p0 = Py2Cpp_double([0.0])
alphas = Py2Cpp_double([1.0])
nu = Py2Cpp_int([0])
num_el_st = 2
el_st = 0
masses = Py2Cpp_double([2000.0])

omega = alphas[0]/masses[0]
k = masses[0] * omega**2
  


# Initialize the grid and do the mappings:
# Wfcgrid2(vector<double>& rmin_, vector<double>& rmax_, vector<double>& dr_, int nstates_)
dx = 0.1
dt = 0.1
wfc = Wfcgrid2(Py2Cpp_double([-15.0]), Py2Cpp_double([15.0]),  Py2Cpp_double([dx]), num_el_st)

# Add a wavefunction in adiabatic basis 
wfc.add_wfc_HO(x0, p0, alphas, el_st, nu, 1.0+0.0j, 1)

# Get Ham in diabatic rep for all points 
# For this Hamiltonian, the adiabatic and diabatic kinetic energies are equal
wfc.update_Hamiltonian(compute_model, {"E_n": [0.0, 0.001], "x_n":[0.0, 0.0], "k_n":[0.001, 0.001]}, 0)

# For this one, they are not - likely because of the notable transformation matrix
#wfc.update_Hamiltonian(compute_model, {"E_n": [0.0, -0.001], "x_n":[0.0, 1.0], "k_n":[0.001, 0.001], "V":0.0001}, 0)

# Compute the dia-to-adi transform matrix
wfc.update_propagator_H(dt)

# Compute the diabatic wfc
wfc.update_diabatic() 

# Compute the reciprocals
wfc.update_reciprocal(0)  # in diabatic  rep
wfc.update_reciprocal(1)  # in adiabatic rep


print( "Norm (dia) = ", wfc.norm(0) )
print( "Norm (adi) = ", wfc.norm(1) )
print( "Ekin (dia) = ", wfc.e_kin(masses, 0) )
print( "Ekin (adi) = ", wfc.e_kin(masses, 1) )
print( "Expected kinetic energy = ", 0.5*alphas[0]/(2.0*masses[0]) )
print( "Epot (dia) = ", wfc.e_pot(0) )
print( "Epot (adi) = ", wfc.e_pot(1) )



