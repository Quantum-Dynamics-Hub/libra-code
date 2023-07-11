#*********************************************************************************
#* Copyright (C) 2018 Brendan A. Smith, Alexey V. Akimov
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

"""
  This example illustrates how to run adiabatic dynamics 
  with a model potential using nHamiltonian class

  Note: in 1 (electronic) state case, the diabatic and adiabatic 
  representations are identical, so one can just populate the diabatic
  properties and don't have to do the "transformation" to the adiabatic
  representation. This is what the first example function demonstrates.

  Or, one can still do the conversion - the result should be the same.
  This is what the second example function demonstrates.

"""



# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *

def example_1():

    # First, create the Hamiltonian object
    ham = nHamiltonian(1,1,1)  

    # Allocate memory and connect the allocated objects to the ham object
    Hdia = CMATRIX(1,1);   Sdia = CMATRIX(1,1);    ham.set_ham_dia_by_ref(Hdia);    ham.set_ovlp_dia_by_ref(Sdia); 
    invSdia = CMATRIX(1,1);
    Hadi = CMATRIX(1,1);   ham.set_ham_adi_by_ref(Hadi);  
    U = CMATRIX(1,1);      ham.set_basis_transform_by_ref(U); 
    Cadi = CMATRIX(1,1);   ham.set_ampl_adi_by_ref(Cadi);
    Cdia = CMATRIX(1,1);   ham.set_ampl_dia_by_ref(Cdia);

    d1ham_dia = CMATRIXList(); d1ham_adi = CMATRIXList()
    dc1_dia = CMATRIXList();   dc1_adi = CMATRIXList()
    for i in [0]:
        d1ham_dia.append( CMATRIX(1,1) ); d1ham_adi.append( CMATRIX(1,1) )
        dc1_dia.append( CMATRIX(1,1) );   dc1_adi.append( CMATRIX(1,1) )

    ham.set_d1ham_dia_by_ref(d1ham_dia);  ham.set_d1ham_adi_by_ref(d1ham_adi)
    ham.set_dc1_dia_by_ref(dc1_dia);      ham.set_dc1_adi_by_ref(dc1_adi)


    # Nuclear DOFs
    mol = Nuclear(1)
    mol.mass[0] = 2000.0
    mol.q[0] = -1.1
    mol.p[0] = 5.0

    Cadi.set(0, 0, 1.0, 0.0)
    Cdia.set(0, 0, 1.0, 0.0)
    params = set_params_1S_1D_poly4("DW:")

    # Initialization
    model_1S_1D_poly4(Hdia, Sdia, d1ham_dia, dc1_dia, mol.q, params)

    f = open("_nucl_nHam_2.1.txt","w")
    f.close()

    dt = 1.0
    for i in xrange(1000):

        mol.propagate_p(0.5*dt)
        mol.propagate_q(dt)

        model_1S_1D_poly4(Hdia, Sdia, d1ham_dia, dc1_dia, mol.q, params)
        mol.f[0] = ham.forces_dia().get(0).real

        mol.propagate_p(0.5*dt)

        epot = Hdia.get(0,0).real
        ekin = compute_kinetic_energy(mol)
 
        f = open("_nucl_nHam_2.1.txt","a")
        f.write("i= %3i q[0]= %8.5f p[0]= %8.5f  ekin= %8.5f  epot= %8.5f  etot= %8.5f\n" % (i, mol.q[0], mol.p[0], ekin, epot, ekin+epot))      
        f.close()


def example_2():

    # First, create the Hamiltonian object
    ham = nHamiltonian(1,1,1)  

    # Allocate memory and connect the allocated objects to the ham object
    Hdia = CMATRIX(1,1);   Sdia = CMATRIX(1,1);    ham.set_ham_dia_by_ref(Hdia);    ham.set_ovlp_dia_by_ref(Sdia); 
    invSdia = CMATRIX(1,1);
    Hadi = CMATRIX(1,1);   ham.set_ham_adi_by_ref(Hadi);  
    U = CMATRIX(1,1);      ham.set_basis_transform_by_ref(U); 
    Cadi = CMATRIX(1,1);   ham.set_ampl_adi_by_ref(Cadi);
    Cdia = CMATRIX(1,1);   ham.set_ampl_dia_by_ref(Cdia);

    d1ham_dia = CMATRIXList(); d1ham_adi = CMATRIXList()
    dc1_dia = CMATRIXList();   dc1_adi = CMATRIXList()
    for i in [0]:
        d1ham_dia.append( CMATRIX(1,1) ); d1ham_adi.append( CMATRIX(1,1) )
        dc1_dia.append( CMATRIX(1,1) );   dc1_adi.append( CMATRIX(1,1) )

    ham.set_d1ham_dia_by_ref(d1ham_dia);  ham.set_d1ham_adi_by_ref(d1ham_adi)
    ham.set_dc1_dia_by_ref(dc1_dia);      ham.set_dc1_adi_by_ref(dc1_adi)


    # Nuclear DOFs
    mol = Nuclear(1)
    mol.mass[0] = 2000.0
    mol.q[0] = -1.1
    mol.p[0] = 5.0

    Cadi.set(0, 0, 1.0, 0.0)
    Cdia.set(0, 0, 1.0, 0.0)
    params = set_params_1S_1D_poly4("DW:")

    # Initialization
    model_1S_1D_poly4(Hdia, Sdia, d1ham_dia, dc1_dia, mol.q, params)

    f = open("_nucl_nHam_2.1.txt","w")
    f.close()

    dt = 1.0
    for i in xrange(1000):

        mol.propagate_p(0.5*dt)
        mol.propagate_q(dt)

        model_1S_1D_poly4(Hdia, Sdia, d1ham_dia, dc1_dia, mol.q, params)
        ham.compute_adiabatic(1)
        mol.f[0] = ham.forces_adi().get(0).real

        mol.propagate_p(0.5*dt)

        epot = Hadi.get(0,0).real
        ekin = compute_kinetic_energy(mol)
 
        f = open("_nucl_nHam_2.1.txt","a")
        f.write("i= %3i q[0]= %8.5f p[0]= %8.5f  ekin= %8.5f  epot= %8.5f  etot= %8.5f\n" % (i, mol.q[0], mol.p[0], ekin, epot, ekin+epot))      
        f.close()



#example_1()
example_2()

