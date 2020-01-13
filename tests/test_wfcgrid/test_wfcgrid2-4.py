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

 This example demonstrates how to use SOFT integrator in diabatic and adiabatic reps

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
dt = 10.0
wfc = Wfcgrid2(Py2Cpp_double([-15.0]), Py2Cpp_double([15.0]),  Py2Cpp_double([dx]), num_el_st)

# Add a wavefunction in adiabatic basis 
wfc.add_wfc_HO(x0, p0, alphas, el_st, nu, 1.0+0.0j, 1)

# Get Ham in diabatic rep for all points 
wfc.update_Hamiltonian(compute_model, {"E_n": [0.0, -0.001], "x_n":[0.0, 1.0], "k_n":[0.001, 0.001], "V":0.0001}, 0)


# Compute the dia-to-adi transform matrix
wfc.update_propagator_H(0.5*dt)

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


# Need the reciprocal space propagator
wfc.update_propagator_K(dt, masses)

nsteps = 1000

f = open("dyn.txt", "w")

for step in range(nsteps):
    wfc.SOFT_propagate()
    wfc.update_reciprocal(0)  # in adiabatic rep

    # Diabatic is the rep used for propagation, so we need to 
    # convert wfcs into adiabatic one
    wfc.update_adiabatic()
    wfc.update_reciprocal(1)  # in adiabatic rep

    Ddia = wfc.get_den_mat(0)  # diabatic density matrix
    Dadi = wfc.get_den_mat(1)  # adiabatic density matrix

    p0_dia = Ddia.get(0,0).real
    p0_adi = Dadi.get(0,0).real

    # Use the adiabatic wavefunction to compute the properties!
    q = wfc.get_pow_q(1, 1).get(0).real
    p = wfc.get_pow_p(1, 1).get(0).real


    f.write( F"step= {step}  Ekin(dia)= {wfc.e_kin(masses, 0)} Ekin(adi)= { wfc.e_kin(masses, 1)} \
              Epot(dia)= {wfc.e_pot(0)} Epot(adi)= {wfc.e_pot(1)} \
              Etot(dia)= {wfc.e_tot(masses, 0)}  Etot(adi)= {wfc.e_tot(masses, 1)} \
              q= {q} p= {p} p0_dia= {p0_dia} p0_adi= {p0_adi}\n" )
            
f.close()

