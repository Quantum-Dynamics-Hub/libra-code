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


# First, create the Hamiltonian
ham = nHamiltonian(2,2,1)  

# Allocate memory and connect the allocated objects to the ham object
Hdia = CMATRIX(2,2);   Sdia = CMATRIX(2,2);    ham.set_ham_dia_by_ref(Hdia);    ham.set_ovlp_dia_by_ref(Sdia); 
invSdia = CMATRIX(2,2);
Hadi = CMATRIX(2,2);   ham.set_ham_adi_by_ref(Hadi);  
U = CMATRIX(2,2);      ham.set_basis_transform_by_ref(U); 
Cdia = CMATRIX(2,1);   ham.set_ampl_dia_by_ref(Cdia)
Cadi = CMATRIX(2,1);   ham.set_ampl_adi_by_ref(Cadi)

d1ham_dia = CMATRIXList(); d1ham_adi = CMATRIXList()
dc1_dia = CMATRIXList();   dc1_adi = CMATRIXList()
for i in [0]:
    d1ham_dia.append( CMATRIX(2,2) ); d1ham_adi.append( CMATRIX(2,2) )
    dc1_dia.append( CMATRIX(2,2) );   dc1_adi.append( CMATRIX(2,2) )

ham.set_d1ham_dia_by_ref(d1ham_dia);  ham.set_d1ham_adi_by_ref(d1ham_adi)
ham.set_dc1_dia_by_ref(dc1_dia);      ham.set_dc1_adi_by_ref(dc1_adi)
params = doubleList()

# Nuclear DOFs
mol = Nuclear(1)
mol.mass[0] = 2000.0
mol.q[0] = -1.1
mol.p[0] = 5.0

# Electronic State - For 1 state PES, set istate = 0
istate = 0

# Initialization
#model_SAC(Hdia, Sdia, d1ham_dia, dc1_dia, mol.q, params)
#model_DAC(Hdia, Sdia, d1ham_dia, dc1_dia, mol.q, params)
model_double_well(Hdia, Sdia, d1ham_dia, dc1_dia, mol.q, params)

print "\nBEGIN TEST 1"
print "\n #### Diabatic Representation #### "
print "Diabatic Hamiltonian"
Hdia.show_matrix()
print "istate = ", istate
epot = Hdia.get(istate, istate).real
mol.f[0] = -d1ham_dia[0].get(istate, istate).real
Sdia.show_matrix(), "\n"
print "Diabatic Potential = ", epot
print "Diabatic Force = ", mol.f[0], "\n"
#sys.exit(0)
print "END TEST 1\n"

print "\nBEGIN TEST 2"
print "\n#### Adiabatic Representation #### "
ham.compute_adiabatic(1)
print "Adiabatic Hamiltonian"
Hadi.show_matrix()
print "istate = ", istate
epot = Hadi.get(istate, istate).real
mol.f[0] = -d1ham_adi[0].get(istate, istate).real
print "Adiabatic Potential = ", epot
print "Adiabatic Force = ", mol.f[0], "\n"
print "Force Testing = ", (mol.q[0] - mol.q[0]**3)
#sys.exit(0)
print "END TEST 2 \n"

f = open("_nucl_nHam_2.1.txt","w")
dt = 1.0
for i in xrange(5000):

    # Start Nuclear dynamics
    mol.propagate_p(0.5*dt)
    mol.propagate_q(dt)

    # Select model
    #model_SAC(Hdia, Sdia, d1ham_dia, dc1_dia, mol.q, params)
    #model_DAC(Hdia, Sdia, d1ham_dia, dc1_dia, mol.q, params)
    model_double_well(Hdia, Sdia, d1ham_dia, dc1_dia, mol.q, params)

    # If adiabatic
    #ham.compute_adiabatic(1)
    #mol.f[0] = -d1ham_adi[0].get(istate, istate).real    
    #epot = Hadi.get(istate, istate).real

    # If diabatic
    mol.f[0] = -d1ham_dia[0].get(istate, istate).real  
    epot = Hdia.get(istate, istate).real

    mol.propagate_p(0.5*dt)
    # End Nuclear dynamics

    ekin = compute_kinetic_energy(mol)
 
    f.write("i= %3i q[0]= %8.5f p[0]= %8.5f  ekin= %8.5f  epot= %8.5f  etot= %8.5f\n" % (i, mol.q[0], mol.p[0], ekin, epot, ekin+epot))

      
f.close()

