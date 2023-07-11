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
sys.path.insert(1,cwd+"/../../_build/src/dyn")

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    #from cyglibra_core import *
    from cygconverters import *
    from cygnhamiltonian_generic import *
    from cyglinalg import *
    from cygdyn import *

elif sys.platform=="linux" or sys.platform=="linux2":
    #from liblibra_core import *
    from libconverters import *
    from libnhamiltonian_generic import *
    from liblinalg import *
    from libgdyn import *



def model1(Hdia, Sdia, d1ham_dia, dc1_dia, x,x0,k,D,V):
    """
              k*x^2         V
    Hdia =       V        k*(x-x0)^2 + D

    Sdia =  I

    Ddia  = 0.0

    """

    Sdia.set(0,0, 1.0+0.0j);  Sdia.set(0,1, 0.0+0.0j);
    Sdia.set(1,0, 0.0+0.0j);  Sdia.set(1,1, 1.0+0.0j);

    Hdia.set(0,0, k*x*x*(1.0+0.0j) );   Hdia.set(0,1, V*(1.0+0.0j));
    Hdia.set(1,0, V*(1.0+0.0j));        Hdia.set(1,1, (k*(x-x0)**2 + D)*(1.0+0.0j));


    for i in [0]:
        #  d Hdia / dR_0
        d1ham_dia[i].set(0,0, 2.0*k*x*(1.0+0.0j) );   d1ham_dia[i].set(0,1, 0.0+0.0j);
        d1ham_dia[i].set(1,0, 0.0+0.0j);   d1ham_dia[i].set(1,1,2.0*k*(x-x0)*(1.0+0.0j));

        #  <dia| d/dR_0| dia >
        dc1_dia[i].set(0,0, 0.0+0.0j);   dc1_dia[i].set(0,1, 0.0+0.0j);
        dc1_dia[i].set(1,0, 0.0+0.0j);   dc1_dia[i].set(1,1, 0.0+0.0j);

    

def model2(Hdia, Sdia, d1ham_dia, dc1_dia, x,x0,k,D,V):
    """
              k*x^2         V
    Hdia =       V        k*(x-x0)^2 + D

    Sdia =  I

    Ddia != 0.0, but Ddia + Ddia.H() = dSdia/dR, with dSdia = 0.0

    """

    Sdia.set(0,0, 1.0+0.0j);  Sdia.set(0,1, 0.0+0.0j);
    Sdia.set(1,0, 0.0+0.0j);  Sdia.set(1,1, 1.0+0.0j);

    Hdia.set(0,0, k*x*x*(1.0+0.0j) );   Hdia.set(0,1, V*(1.0+0.0j));
    Hdia.set(1,0, V*(1.0+0.0j));        Hdia.set(1,1, (k*(x-x0)**2 + D)*(1.0+0.0j));


    for i in [0]:
        #  d Hdia / dR_0
        d1ham_dia[i].set(0,0, 2.0*k*x*(1.0+0.0j) );   d1ham_dia[i].set(0,1, 0.0+0.0j);
        d1ham_dia[i].set(1,0, 0.0+0.0j);   d1ham_dia[i].set(1,1,2.0*k*(x-x0)*(1.0+0.0j));

        #  <dia| d/dR_0| dia >
        dc1_dia[i].set(0,0, 0.0+0.0j);   dc1_dia[i].set(0,1,-0.1+0.0j);
        dc1_dia[i].set(1,0, 0.1+0.0j);   dc1_dia[i].set(1,1, 0.0+0.0j);






def model3(Hdia, Sdia, d1ham_dia, dc1_dia, x,x0,k,D,V):
    """
              k*x^2         V
    Hdia =       V        k*(x-x0)^2 + D

    Sdia =    1             0.5*exp(-(x-x0)^2)
             0.5*exp(-(x-x0)^2)       1

    Ddia != 0.0, but Ddia + Ddia.H() = dSdia/dR, with dSdia !=0.0

    """

    Sdia.set(0,0, 1.0+0.0j);                          Sdia.set(0,1, math.exp(-(x-x0)**2)*(0.05+0.0j));
    Sdia.set(1,0, math.exp(-(x-x0)**2)*(0.05+0.0j));  Sdia.set(1,1, 1.0+0.0j);

    Hdia.set(0,0, k*x*x*(1.0+0.0j) );   Hdia.set(0,1, V*(1.0+0.0j));
    Hdia.set(1,0, V*(1.0+0.0j));        Hdia.set(1,1, (k*(x-x0)**2 + D)*(1.0+0.0j));


    for i in [0]:
        #  d Hdia / dR_0
        d1ham_dia[i].set(0,0, 2.0*k*x*(1.0+0.0j) );  d1ham_dia[i].set(0,1, 0.0+0.0j);
        d1ham_dia[i].set(1,0, 0.0+0.0j);             d1ham_dia[i].set(1,1,2.0*k*(x-x0)*(1.0+0.0j));

        #  <dia| d/dR_0| dia >
        d = -(x-x0)*math.exp(-(x-x0)**2)*(0.05+0.0j)
        dc1_dia[i].set(0,0, 0.0+0.0j);  dc1_dia[i].set(0,1, d);
        dc1_dia[i].set(1,0,-d);         dc1_dia[i].set(1,1, 0.0+0.0j);

    


def run_test_adi(model):

    ham = nHamiltonian(2,2,1)  


    # Allocate memory
    Hdia = CMATRIX(2,2);   Sdia = CMATRIX(2,2);    ham.set_ham_dia_by_ref(Hdia);    ham.set_ovlp_dia_by_ref(Sdia); 
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



    #  Set up the models and compute internal variables
    # Initialization
    Cadi.set(0, 1.0+0.0j); Cadi.set(1, 1.0+0.0j)        
    Cadi *= (1.0/math.sqrt(2.0))
  
    x0, k, D, V, dt  = 1.0, 1.0, -0.1, 0.1, 0.01
    q, p, m = 0.1, 0.0, 100.0

    if model==1:
        model1(Hdia, Sdia, d1ham_dia, dc1_dia, q, x0, k, D, V)
    elif model==2:
        model2(Hdia, Sdia, d1ham_dia, dc1_dia, q, x0, k, D, V)
    elif model==3:
        model3(Hdia, Sdia, d1ham_dia, dc1_dia, q, x0, k, D, V)


    ham.compute_adiabatic(1); 
    f = ham.Ehrenfest_forces_adi().get(0).real

    Etot = 0.5*p*p/m + ham.Ehrenfest_energy_adi().real

    Hvib = Hadi - 1.0j*(p/m)*dc1_adi[0]
#    tmp = U * Cadi;      Cdia.set(0, tmp.get(0)); Cdia.set(1, tmp.get(1));

    # Do the propagation
    for i in xrange(50000):

        propagate_electronic(0.5*dt, Cadi, Hvib)        
        p = p + 0.5*f*dt
        q = q + dt*p/m


        if model==1:
            model1(Hdia, Sdia, d1ham_dia, dc1_dia, q, x0, k, D, V)
        elif model==2:
            model2(Hdia, Sdia, d1ham_dia, dc1_dia, q, x0, k, D, V)
        elif model==3:
            model3(Hdia, Sdia, d1ham_dia, dc1_dia, q, x0, k, D, V)

        ham.compute_adiabatic(1); 
        f = ham.Ehrenfest_forces_adi().get(0).real

        Hvib = Hadi - 1.0j*(p/m)*dc1_adi[0]
 
        p = p + 0.5*f*dt
        propagate_electronic(0.5*dt, Cadi, Hvib)


        Etot = 0.5*p*p/m + ham.Ehrenfest_energy_adi().real

        dm = Cadi * Cadi.H()

        print i, Etot, 0.5*p*p/m, ham.Ehrenfest_energy_adi().real, dm.get(0,0).real, dm.get(1,1).real
    







def run_test_dia(model):

    ham = nHamiltonian(2,2,1)  


    # Allocate memory
    Hdia = CMATRIX(2,2);   Sdia = CMATRIX(2,2);    ham.set_ham_dia_by_ref(Hdia);    ham.set_ovlp_dia_by_ref(Sdia); 
    Hadi = CMATRIX(2,2);   ham.set_ham_adi_by_ref(Hadi);  
    U = CMATRIX(2,2);      ham.set_basis_transform_by_ref(U); 
    Cdia = CMATRIX(2,1);   ham.set_ampl_dia_by_ref(Cdia)
    Cadi = CMATRIX(2,1);   ham.set_ampl_adi_by_ref(Cadi)
    Dadi = CMATRIX(2,2);   ham.set_den_mat_adi_by_ref(Dadi)
    Ddia = CMATRIX(2,2);   ham.set_den_mat_dia_by_ref(Ddia)

    d1ham_dia = CMATRIXList(); d1ham_adi = CMATRIXList()
    dc1_dia = CMATRIXList();   dc1_adi = CMATRIXList()
    for i in [0]:
        d1ham_dia.append( CMATRIX(2,2) ); d1ham_adi.append( CMATRIX(2,2) )
        dc1_dia.append( CMATRIX(2,2) );   dc1_adi.append( CMATRIX(2,2) )
    
    ham.set_d1ham_dia_by_ref(d1ham_dia);  ham.set_d1ham_adi_by_ref(d1ham_adi)
    ham.set_dc1_dia_by_ref(dc1_dia);      ham.set_dc1_adi_by_ref(dc1_adi)



    #  Set up the models and compute internal variables
    # Initialization
    Cadi.set(0, 1.0+0.0j); Cadi.set(1, 1.0+0.0j)        
    Cadi *= (1.0/math.sqrt(2.0))
  
    x0, k, D, V, dt  = 1.0, 1.0, -0.1, 0.1, 0.01
    q, p, m = 0.1, 1.0, 100.0

    if model==1:
        model1(Hdia, Sdia, d1ham_dia, dc1_dia, q, x0, k, D, V)
    elif model==2:
        model2(Hdia, Sdia, d1ham_dia, dc1_dia, q, x0, k, D, V)
    elif model==3:
        model3(Hdia, Sdia, d1ham_dia, dc1_dia, q, x0, k, D, V)

    ham.compute_adiabatic(1); 
    ham.ampl_adi2dia()



    f = ham.Ehrenfest_forces_dia().get(0).real
#    f = ham.Ehrenfest_forces_adi().get(0).real
    Etot = 0.5*p*p/m + ham.Ehrenfest_energy_dia().real


    Hvib = Hdia - 1.0j*(p/m)*dc1_dia[0]   # non-Hermitian


    # Do the propagation
    for i in xrange(25000):

#        propagate_electronic_nonHermitian(0.5*dt, Cdia, Hvib)        
        propagate_electronic(0.5*dt, Cdia, Hvib)
        ham.ampl_dia2adi()


        p = p + 0.5*f*dt
        q = q + dt*p/m

        if model==1:
            model1(Hdia, Sdia, d1ham_dia, dc1_dia, q, x0, k, D, V)
        elif model==2:
            model2(Hdia, Sdia, d1ham_dia, dc1_dia, q, x0, k, D, V)
        elif model==3:
            model3(Hdia, Sdia, d1ham_dia, dc1_dia, q, x0, k, D, V)

        ham.compute_adiabatic(1); 
        ham.ampl_dia2adi()

        Hvib = Hdia - 1.0j*(p/m)*dc1_dia[0]


        f = ham.Ehrenfest_forces_dia().get(0).real
#        f = ham.Ehrenfest_forces_adi().get(0).real
 
        p = p + 0.5*f*dt
#        propagate_electronic_nonHermitian(0.5*dt, Cdia, Hvib)
        propagate_electronic(0.5*dt, Cdia, Hvib)
        ham.ampl_dia2adi()

        Dt = U.H() * dc1_dia[0].H() * U
        LHS = U.H() * d1ham_dia[0] * U - (Dt * Hadi + Hadi * Dt.H())        
        RHS = d1ham_adi[0] - (dc1_adi[0].H() * Hadi + Hadi * dc1_adi[0])

#        print i
#        LHS.show_matrix(); RHS.show_matrix()
#        print "Ehrenfest force adi = ", ( Cadi.H() * RHS * Cadi ).get(0).real
#        print "Ehrenfest force dia = ", ( Cdia.H() * LHS * Cdia ).get(0).real


        Etot = 0.5*p*p/m + ham.Ehrenfest_energy_dia().real

#        dm = Cdia * Cdia.H()
        dm = Cadi * Cadi.H()

        print i, Etot, 0.5*p*p/m, ham.Ehrenfest_energy_dia().real, dm.get(0,0).real, dm.get(1,1).real, ham.Ehrenfest_forces_dia().get(0).real, ham.Ehrenfest_forces_adi().get(0).real


model = 3

run_test_adi(model)
#run_test_dia(model)
        
