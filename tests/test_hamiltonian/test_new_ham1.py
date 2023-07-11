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
import cmath
import math
import os
import sys
import unittest


cwd = os.getcwd()
print "Current working directory", cwd
sys.path.insert(1,cwd+"/../../_build/src/hamiltonian/nHamiltonian_Generic")
sys.path.insert(1,cwd+"/../../_build/src/converters")
sys.path.insert(1,cwd+"/../../_build/src/math_linalg")

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    #from cyglibra_core import *
    from cygconverters import *
    from cygnhamiltonian_generic import *
    from cyglinalg import *

elif sys.platform=="linux" or sys.platform=="linux2":
    #from liblibra_core import *
    from libconverters import *
    from libnhamiltonian_generic import *
    from liblinalg import *





Hdia = CMATRIX(2,2)
Hdia.set(0,0, 1.0+0.0j);  Hdia.set(0,1, 0.1+0.0j);
Hdia.set(1,0, 0.1+0.0j);  Hdia.set(1,1, 1.0+0.0j);

Sdia = CMATRIX(2,2)
Sdia.set(0,0, 1.0+0.0j);  Sdia.set(0,1, 0.0+0.0j);
Sdia.set(1,0, 0.0+0.0j);  Sdia.set(1,1, 1.0+0.0j);


Hadi = CMATRIX(2,2)
U = CMATRIX(2,2)


ham = nHamiltonian(2,2,1)  

opt = 1

if opt==0:
    print "Hdia\n"; ham.set_ham_dia_by_val(Hdia);  ham.get_ham_dia().show_matrix()
    print "Hadi\n"; ham.set_ham_adi_by_val(Hadi);  ham.get_ham_adi().show_matrix()
    print "Sdia\n"; ham.set_ovlp_dia_by_val(Sdia); ham.get_ovlp_dia().show_matrix()
    print "U\n"; ham.set_basis_transform_by_val(U); ham.get_basis_transform().show_matrix()

elif opt==1:
    print "Hdia\n"; ham.set_ham_dia_by_ref(Hdia);  ham.get_ham_dia().show_matrix()
    print "Hadi\n"; ham.set_ham_adi_by_ref(Hadi);  ham.get_ham_adi().show_matrix()
    print "Sdia\n"; ham.set_ovlp_dia_by_ref(Sdia); ham.get_ovlp_dia().show_matrix()
    print "U\n"; ham.set_basis_transform_by_ref(U); ham.get_basis_transform().show_matrix()


ham.compute_adiabatic(0); 

print "Hadi\n"; Hadi.show_matrix(); ham.get_ham_adi().show_matrix()
print Hadi
#Hadi.show_matrix_address()
#ham.get_ham_adi().show_matrix_address()
print "U\n"; U.show_matrix();  ham.get_basis_transform().show_matrix(); 


Hdia.set(0,1, 0.5+0.0j)
Hdia.set(1,0, 0.5+0.0j)
ham.compute_adiabatic(0); 

print "Hadi\n"; Hadi.show_matrix(); ham.get_ham_adi().show_matrix()
print "U\n"; U.show_matrix();  ham.get_basis_transform().show_matrix(); 
print "Hadi="; Hadi.show_matrix()
#Hadi.show_matrix_address()
#ham.get_ham_adi().show_matrix_address()


Cdia = CMATRIX(2,1)
Cdia.set(0, 1.0+0.0j); Cdia.set(1, 0.0+0.0j)
ham.set_ampl_dia_by_ref(Cdia)
print ham.Ehrenfest_energy_dia(); 



#U.show_matrix()






