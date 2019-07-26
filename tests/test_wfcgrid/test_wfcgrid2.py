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

 This example demonstrates the functions of the class Wfcgrid2
       
   defined in:   dyn/wfcgrid2/Wfcgrid2.h


"""

import os
import sys
import math
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *


class tmp:
    pass

def harmonic(q, params):
    """
    1D Harmonic potential 
    """
  
    x = q.get(0)
    k = params["k"]
   
    obj = tmp()
    obj.ham_dia = CMATRIX(1,1)    
    obj.ham_dia.set(0,0, 0.5*k*x**2)

    return obj
    

"""
 Here, we test simple Harmonic oscillator eigenfunctions and will 
 compare the energies as computed by Libra to the analytic results
"""

x0 = Py2Cpp_double([0.0])
p0 = Py2Cpp_double([0.0])
alphas = Py2Cpp_double([1.0])
nu = Py2Cpp_int([0])
num_el_st = 1
el_st = 0
rep = 0  
masses = Py2Cpp_double([2000.0])

omega = alphas[0]/masses[0]
k = masses[0] * omega**2

# For dynamics
dt = 1.0
nsteps = 100
  

# Initialize the grid and do the mappings:
# Wfcgrid2(vector<double>& rmin_, vector<double>& rmax_, vector<double>& dr_, int nstates_)
wfc = Wfcgrid2(Py2Cpp_double([-15.0]), Py2Cpp_double([15.0]),  Py2Cpp_double([0.01]), num_el_st)

# Add a wavefunction:
# void add_wfc_HO(vector<double>& x0, vector<double>& px0, vector<double>& alpha, int init_state, vector<int>& nu, complex<double> weight, int rep);
weight = 1.0+0.0j
wfc.add_wfc_HO(x0, p0, alphas, el_st, nu, weight, rep)

wfc.update_reciprocal(rep)

wfc.update_Hamiltonian(harmonic, {"k": k}, rep)

print "Norm = ", wfc.norm(rep)
print "Ekin = ", wfc.e_kin(masses, rep)
print "Expected kinetic energy = ", 0.5*alphas[0]/(2.0*masses[0])
print "Epot = ", wfc.e_pot(rep)
print "Expected potential energy = ", (0.5*k/alphas[0])*(0.5 + nu[0])
print "Etot = ", wfc.e_tot(masses, rep)
print "Expected total energy = ", omega*(0.5 + nu[0])


"""
 Compare the total energy of various HO states to the analytic results
"""

for n in [0, 1, 2, 3, 10, 20]:

    wfc = Wfcgrid2(Py2Cpp_double([-15.0]), Py2Cpp_double([15.0]),  Py2Cpp_double([0.01]), num_el_st)
    
    nu = Py2Cpp_int([n]) 
    wfc.add_wfc_HO(x0, p0, alphas, el_st, nu, 1.0+0.0j, rep)

    wfc.update_reciprocal(rep)
    wfc.update_Hamiltonian(harmonic, {"k": k}, rep)

    print "========== State %i ==============" % (n)
    print "Etot = ", wfc.e_tot(masses, rep)
    print "Expected total energy = ", omega*(0.5 + nu[0])



"""
 Here, lets demonstrate some dynamics using SOFT propagator
"""

wfc = Wfcgrid2(Py2Cpp_double([-15.0]), Py2Cpp_double([15.0]),  Py2Cpp_double([0.01]), num_el_st)
    
wfc.add_wfc_HO(x0, p0, alphas, el_st, Py2Cpp_int([0]) , 1.0+0.0j, rep)
wfc.add_wfc_HO(x0, p0, alphas, el_st, Py2Cpp_int([1]) , 1.0+0.0j, rep)

wfc.update_reciprocal(rep)
wfc.update_Hamiltonian(harmonic, {"k": k}, rep)


wfc.update_propagator_H(0.5*dt)
wfc.update_propagator_K(dt, masses)

for step in xrange(nsteps):
    wfc.SOFT_propagate()
    q = wfc.get_pow_q(0, 1).get(0).real
    p = wfc.get_pow_p(0, 1).get(0).real

    # Diabatic is the rep used for propagation, so we need to 
    # convert wfcs into adiabatic one
    wfc.update_adiabatic()

    Ddia = wfc.get_den_mat(0)  # diabatic density matrix
    Dadi = wfc.get_den_mat(1)  # adiabatic density matrix

    p0_dia = Ddia.get(0,0).real
    p0_adi = Dadi.get(0,0).real
    print "step= ", step, " Ekin= ", wfc.e_kin(masses, rep), " Epot= ", wfc.e_pot(rep), " Etot= ", wfc.e_tot(masses, rep), " q= ", q, " p= ", p, " p0_dia= ", p0_dia, " p0_adi= ", p0_adi
            





