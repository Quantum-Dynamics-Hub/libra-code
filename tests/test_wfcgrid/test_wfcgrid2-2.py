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

 This example demonstrates the Colbert-Miller approach to computing

 <psi|T|psi> matrix elements as implemented in the class Wfcgrid2
       
   defined in:   dyn/wfcgrid2/Wfcgrid2.h


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



class tmp:
    pass

#========================== 1D case ====================================


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
rep = 1
masses = Py2Cpp_double([2000.0])

omega = alphas[0]/masses[0]
k = masses[0] * omega**2
  


# Initialize the grid and do the mappings:
# Wfcgrid2(vector<double>& rmin_, vector<double>& rmax_, vector<double>& dr_, int nstates_)
dx = 0.1
wfc = Wfcgrid2(Py2Cpp_double([-15.0]), Py2Cpp_double([15.0]),  Py2Cpp_double([dx]), num_el_st)

# Add a wavefunction:
weight = 1.0+0.0j
wfc.add_wfc_HO(x0, p0, alphas, el_st, nu, weight, rep)

start = time.time()
wfc.update_reciprocal(rep)
end = time.time()
print("Time to run reciprocal = ", end - start)


wfc.update_Hamiltonian(harmonic, {"k": k}, rep)

psi = wfc.PSI_adi

start = time.time()
Tpsi = wfc.T_PSI_adi(Py2Cpp_int([1]), masses, 1.0+0.0j)
end = time.time()
print("Time to run T_PSI_adi = ", end - start)


# Integral Ekin = <psi|T|psi>
start = time.time()
npts = len(psi)
Ekin = 0.0+0.0j
f = open("wfc.txt", "w")
for ipt in range(npts):    
    f.write(F"i = {ipt} psi = {psi[ipt].get(0)}   Tpsi = {Tpsi[ipt].get(0)}\n")    
    Ekin += (psi[ipt].H() * Tpsi[ipt]).get(0, 0)     
f.close()

end = time.time()
print("Time to run integrals and printout = ", end - start)


Ekin = Ekin * dx

print("Ekin(Colbert-Miller) = ", Ekin)
print( "Norm = ", wfc.norm(rep) )
print( "Ekin = ", wfc.e_kin(masses, rep) )
print( "Expected kinetic energy = ", 0.5*alphas[0]/(2.0*masses[0]) )
print( "Epot = ", wfc.e_pot(rep) )


#========================== 2D case ====================================
#
#  WARNING!!! This test does work, but it is pretty slow because of the scaling
#
#  so use a setting with very small number of points (e.g. 128 per dimension)
#
#

def harmonic2D(q, params):
    """
    2D Harmonic potential 
    """
  
    x = q.get(0)
    y = q.get(1)
    kx = params["kx"]
    ky = params["ky"]
   
    obj = tmp()
    obj.ham_dia = CMATRIX(1,1)    
    obj.ham_dia.set(0, 0, (0.5*kx*x**2 + 0.5*ky*y**2)*(1.0+0.0j) )

    return obj
    

"""
 Here, we test simple Harmonic oscillator eigenfunctions and will 
 compare the energies as computed by Libra to the analytic results
"""

x0 = Py2Cpp_double([0.0, 0.0])
p0 = Py2Cpp_double([0.0, 0.0])
alphas = Py2Cpp_double([1.0, 2.0])
nu = Py2Cpp_int([0, 0])
num_el_st = 1
el_st = 0
rep = 1
masses = Py2Cpp_double([2000.0, 2000.0])

omega_x = alphas[0]/masses[0]
kx = masses[0] * omega_x**2

omega_y = alphas[1]/masses[1]
ky = masses[1] * omega_y**2

  


# Initialize the grid and do the mappings:
# Wfcgrid2(vector<double>& rmin_, vector<double>& rmax_, vector<double>& dr_, int nstates_)
dx = 0.1
dy = 0.1
wfc = Wfcgrid2(Py2Cpp_double([-5.0, -5.0]), Py2Cpp_double([5.0, 5.0]),  Py2Cpp_double([dx, dy]), num_el_st)


# Add a wavefunction:
weight = 1.0+0.0j
wfc.add_wfc_HO(x0, p0, alphas, el_st, nu, weight, rep)

start = time.time()
wfc.update_reciprocal(rep)
end = time.time()
print("Time to run reciprocal wfc calculations = ", end - start)
#sys.exit(0)

wfc.update_Hamiltonian(harmonic2D, {"kx": kx, "ky":ky}, rep)



psi = wfc.PSI_adi

start = time.time()
Tpsi = wfc.T_PSI_adi(Py2Cpp_int([1, 1]), masses, 1.0+0.0j)
end = time.time()
print("Time to run T_PSI_adi = ", end - start)




# Integral Ekin = <psi|T|psi>
start = time.time()
npts = len(psi)
Ekin = 0.0+0.0j
f = open("wfc2D.txt", "w")
for ipt in range(npts):    
    f.write(F"i = {ipt} psi = {psi[ipt].get(0)}   Tpsi = {Tpsi[ipt].get(0)}\n")    
    Ekin += (psi[ipt].H() * Tpsi[ipt]).get(0, 0)     
f.close()
Ekin = Ekin * dx * dy
end = time.time()
print("Time to run run integral & printout = ", end - start)


print("Ekin(Colbert-Miller) = ", Ekin)
print( "Norm = ", wfc.norm(rep) )
print( "Ekin = ", wfc.e_kin(masses, rep) )
print( "Expected kinetic energy = ", 0.5*alphas[0]/(2.0*masses[0]) + 0.5*alphas[1]/(2.0*masses[1]) )
print( "Epot = ", wfc.e_pot(rep) )


