#*********************************************************************************
#* Copyright (C) 2017-2018 Alexey V. Akimov
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

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from libra_py import *




def model1(Hdia, Sdia, d1ham_dia, dc1_dia, x, params):
    """
              k*x^2         V
    Hdia =       V        k*(x-x0)^2 + D

    Sdia =  I

    Ddia  = 0.0

    """

    x0,k,D,V = params["x0"], params["k"], params["D"], params["V"]

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

    

def model2(Hdia, Sdia, d1ham_dia, dc1_dia, x, params):
    """
              k*x^2         V
    Hdia =       V        k*(x-x0)^2 + D

    Sdia =  I

    Ddia != 0.0, but Ddia + Ddia.H() = dSdia/dR, with dSdia = 0.0

    """

    x0,k,D,V = params["x0"], params["k"], params["D"], params["V"]

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



def model3(Hdia, Sdia, d1ham_dia, dc1_dia, x, params):
    """
              k*x^2         V
    Hdia =       V        k*(x-x0)^2 + D

    Sdia =    1                           0.05*exp(-(x-0.5*x0)^2)
             0.05*exp(-(x-0.5*x0)^2)               1

    Ddia != 0.0, but Ddia + Ddia.H() = dSdia/dR, with dSdia !=0.0

    """

    x0,k,D,V = params["x0"], params["k"], params["D"], params["V"]

    ex = math.exp(-(x-0.5*x0)**2)*(0.05+0.0j)

    Sdia.set(0,0, 1.0+0.0j);     Sdia.set(0,1, ex );
    Sdia.set(1,0, ex);           Sdia.set(1,1, 1.0+0.0j);

    Hdia.set(0,0, k*x*x*(1.0+0.0j) );   Hdia.set(0,1, V*(1.0+0.0j));
    Hdia.set(1,0, V*(1.0+0.0j));        Hdia.set(1,1, (k*(x-x0)**2 + D)*(1.0+0.0j));


    for i in [0]:
        #  d Hdia / dR_0
        d1ham_dia[i].set(0,0, 2.0*k*x*(1.0+0.0j) );  d1ham_dia[i].set(0,1, 0.0+0.0j);
        d1ham_dia[i].set(1,0, 0.0+0.0j);             d1ham_dia[i].set(1,1,2.0*k*(x-x0)*(1.0+0.0j));

        #  <dia| d/dR_0| dia >
        d = -(x-0.5*x0)*ex
        dc1_dia[i].set(0,0, 0.0+0.0j);  dc1_dia[i].set(0,1, d);
        dc1_dia[i].set(1,0, d);         dc1_dia[i].set(1,1, 0.0+0.0j);

    

def model4(Hdia, Sdia, d1ham_dia, dc1_dia, x, params):
    """
              k*cos(w*x)         V
    Hdia =       V        k*sin(w*x) + D

    Sdia =  I

    Ddia  = 0.0

    """

    k, w, V = params["k"], params["omega"], params["V"]

    Sdia.set(0,0, 1.0+0.0j);  Sdia.set(0,1, 0.0+0.0j);
    Sdia.set(1,0, 0.0+0.0j);  Sdia.set(1,1, 1.0+0.0j);

    Hdia.set(0,0, k*math.cos(x*w)*(1.0+0.0j) );   Hdia.set(0,1, V*(1.0+0.0j));
    Hdia.set(1,0, V*(1.0+0.0j));                  Hdia.set(1,1, k*math.sin(x*w)*(1.0+0.0j));


    for i in [0]:
        #  d Hdia / dR_0
        d1ham_dia[i].set(0,0,-w*k*math.sin(x*w)*(1.0+0.0j) );   d1ham_dia[i].set(0,1, 0.0+0.0j);
        d1ham_dia[i].set(1,0, 0.0+0.0j);                        d1ham_dia[i].set(1,1, w*k*math.cos(x*w)*(1.0+0.0j));

        #  <dia| d/dR_0| dia >
        dc1_dia[i].set(0,0, 0.0+0.0j);   dc1_dia[i].set(0,1, 0.0+0.0j);
        dc1_dia[i].set(1,0, 0.0+0.0j);   dc1_dia[i].set(1,1, 0.0+0.0j);




def compute_model(model, Hdia, Sdia, d1ham_dia, dc1_dia, q, params):

   if model==1:
       model1(Hdia, Sdia, d1ham_dia, dc1_dia, q, params)
   elif model==2:
       model2(Hdia, Sdia, d1ham_dia, dc1_dia, q, params)
   elif model==3:
       model3(Hdia, Sdia, d1ham_dia, dc1_dia, q, params)
   elif model==4:
       model4(Hdia, Sdia, d1ham_dia, dc1_dia, q, params)



def propagate_el(Cdia, Hvib, Sdia, Cadi, dt, rep):

    if rep==0:
        propagate_electronic(0.5*dt, Cdia, Hvib, Sdia)        
    elif rep==1:
        propagate_electronic(0.5*dt, Cadi, Hvib)        


def compute_etot(ham, p, m, rep):
    Etot = 0.0
    if rep==0:
        Etot = 0.5*p*p/m + ham.Ehrenfest_energy_dia().real 
    
    elif rep==1:
        Etot = 0.5*p*p/m + ham.Ehrenfest_energy_adi().real        

    return Etot


def compute_frc(ham, rep):
    f = None
    if rep==0:
        f = ham.Ehrenfest_forces_dia().get(0).real
    
    elif rep==1:
        f = ham.Ehrenfest_forces_adi().get(0).real

    return f


def compute_Hvib(Hdia, Hadi, dc1_dia, dc1_adi, p, m, rep):
    Hvib = None

    if rep==0:
        Hvib = (Hdia - 1.0j*(p/m)*dc1_dia[0])
        
    elif rep==1:
        Hvib = Hadi - 1.0j*(p/m)*dc1_adi[0]

    return Hvib


def run_test(model, rep, outname):
    """
    model - setup the Hamiltonian
    rep - 0 - diabatic, 1 - adiabatic
    outname - the name of the output file
    """

    ham = nHamiltonian(2,2,1)  

    # Allocate memory
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


    #  Set up the models and compute internal variables
    # Initialization
    # Model parameters 
    params = {}
    params["x0"], params["k"], params["D"], params["V"] = 1.0, 0.1, -0.1, 0.05
    params["omega"] = 0.25

    # Simulation parameters
    dt = 1.0

    # Dynamical variables and system-specific properties
    Cadi.set(0, 1.0+0.0j); Cadi.set(1, 1.0+0.0j);   Cadi *= (1.0/math.sqrt(2.0))  
    q, p, m = 0.1, 0.0, 100.0


    # Initial calculations
    compute_model(model, Hdia, Sdia, d1ham_dia, dc1_dia, q, params)        
    ham.compute_adiabatic(1); 
    ham.ampl_adi2dia()


    dm_dia, dm_adi = None, None


    Etot = compute_etot(ham, p, m, rep)
    f = compute_frc(ham, rep)
    Hvib = compute_Hvib(Hdia, Hadi, dc1_dia, dc1_adi, p, m, rep)


    out = open(outname, "w")
    out.close()


    # Do the propagation
    for i in xrange(500):

        propagate_el(Cdia, Hvib, Sdia, Cadi, dt, rep)
    
        p = p + 0.5*f*dt
        q = q + dt*p/m
    
    
        compute_model(model, Hdia, Sdia, d1ham_dia, dc1_dia, q, params)        
        ham.compute_adiabatic(1); 
        

        Etot = compute_etot(ham, p, m, rep)
        f = compute_frc(ham, rep)
        Hvib = compute_Hvib(Hdia, Hadi, dc1_dia, dc1_adi, p, m, rep)
            
        p = p + 0.5*f*dt

        if rep==0:
            Hvib = (Hdia - 1.0j*(p/m)*dc1_dia[0])        
        elif rep==1:
            Hvib = Hadi - 1.0j*(p/m)*dc1_adi[0]

        propagate_el(Cdia, Hvib, Sdia, Cadi, dt, rep)
     

        #=========== Properties ==========

        if rep==0:
            dm_dia = Sdia * Cdia * Cdia.H() * Sdia
            dm_adi = U.H() * dm_dia * U
           
        elif rep==1:
            dm_adi = Cadi * Cadi.H()
            su = Sdia * U
            dm_dia = su * dm_adi * su.H()

        Etot = compute_etot(ham, p, m, rep)


        out = open(outname, "a")
        ret = (i*dt, q, p, Etot, 0.5*p*p/m, dm_adi.get(0,0).real, dm_adi.get(1,1).real, dm_dia.get(0,0).real, dm_dia.get(1,1).real,
               Hdia.get(0,0).real, Hdia.get(1,1).real, Hadi.get(0,0).real, Hadi.get(1,1).real, dc1_adi[0].get(0,1).real
              )
        out.write("%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f\n" %  ret )
        out.close()



model = 4

run_test(model, 0, "0.txt")
#run_test(model, 1, "1.txt")

        
