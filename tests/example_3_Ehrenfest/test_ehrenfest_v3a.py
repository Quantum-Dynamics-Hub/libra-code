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

class tmp:
    pass    


def model1(q, params):
    """
              k*x^2         V
    Hdia =       V        k*(x-x0)^2 + D

    Sdia =  I

    Ddia  = 0.0

    """

    Hdia = CMATRIX(2,2)
    Sdia = CMATRIX(2,2)
    d1ham_dia = CMATRIXList();  d1ham_dia.append( CMATRIX(2,2))
    dc1_dia = CMATRIXList();  dc1_dia.append( CMATRIX(2,2))
  

    x = q.get(0)
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


    obj = tmp()
    obj.ham_dia = Hdia
    obj.ovlp_dia = Sdia
    obj.d1ham_dia = d1ham_dia
    obj.dc1_dia = dc1_dia

    return obj

    

def model2(q, params):
    """
              k*x^2         V
    Hdia =       V        k*(x-x0)^2 + D

    Sdia =  I

    Ddia != 0.0, but Ddia + Ddia.H() = dSdia/dR, with dSdia = 0.0

    """

    Hdia = CMATRIX(2,2)
    Sdia = CMATRIX(2,2)
    d1ham_dia = CMATRIXList();  d1ham_dia.append( CMATRIX(2,2))
    dc1_dia = CMATRIXList();  dc1_dia.append( CMATRIX(2,2))
  

    x = q.get(0)
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


    obj = tmp()
    obj.ham_dia = Hdia
    obj.ovlp_dia = Sdia
    obj.d1ham_dia = d1ham_dia
    obj.dc1_dia = dc1_dia

    return obj



def model3(q, params):
    """
              k*x^2         V
    Hdia =       V        k*(x-x0)^2 + D

    Sdia =    1                           0.05*exp(-(x-0.5*x0)^2)
             0.05*exp(-(x-0.5*x0)^2)               1

    Ddia != 0.0, but Ddia + Ddia.H() = dSdia/dR, with dSdia !=0.0

    """

    Hdia = CMATRIX(2,2)
    Sdia = CMATRIX(2,2)
    d1ham_dia = CMATRIXList();  d1ham_dia.append( CMATRIX(2,2))
    dc1_dia = CMATRIXList();  dc1_dia.append( CMATRIX(2,2))
  

    x = q.get(0)
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


    obj = tmp()
    obj.ham_dia = Hdia
    obj.ovlp_dia = Sdia
    obj.d1ham_dia = d1ham_dia
    obj.dc1_dia = dc1_dia

    return obj
    

def model4(q, params):
    """
              k*cos(w*x)         V
    Hdia =       V        k*sin(w*x) + D

    Sdia =  I

    Ddia  = 0.0

    """

    Hdia = CMATRIX(2,2)
    Sdia = CMATRIX(2,2)
    d1ham_dia = CMATRIXList();  d1ham_dia.append( CMATRIX(2,2))
    dc1_dia = CMATRIXList();  dc1_dia.append( CMATRIX(2,2))
  

    x = q.get(0)
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


    obj = tmp()
    obj.ham_dia = Hdia
    obj.ovlp_dia = Sdia
    obj.d1ham_dia = d1ham_dia
    obj.dc1_dia = dc1_dia

    return obj




def compute_model(q, params):

    model = params["model"]
    res = None

    if model==1:
        res = model1(q, params)
    elif model==2:
        res = model2(q, params)
    elif model==3:
        res = model3(q, params)
    elif model==4:
        res = model4(q, params)

    res.rep = params["rep"]
    return res
    


def compute_etot(ham, p, iM, rep):
    Etot = 0.0
    if rep==0:
        Etot = 0.5*(p.T()*iM*p).get(0) + ham.Ehrenfest_energy_dia().real 
    
    elif rep==1:
        Etot = 0.5*(p.T()*iM*p).get(0) + ham.Ehrenfest_energy_adi().real        

    return Etot





def run_test(model, rep, outname):
    """
    model - setup the Hamiltonian
    rep - 0 - diabatic, 1 - adiabatic
    outname - the name of the output file
    """

    ham = nHamiltonian(2,2,1)  

    # Allocate memory
    Hdia = CMATRIX(2,2);   Sdia = CMATRIX(2,2);    ham.set_ham_dia_by_ref(Hdia);    ham.set_ovlp_dia_by_ref(Sdia); 
    NACdia = CMATRIX(2,2); NACadi = CMATRIX(2,2);  ham.set_nac_dia_by_ref(NACdia);  ham.set_nac_adi_by_ref(NACadi);
    Hvibdia = CMATRIX(2,2); Hvibadi = CMATRIX(2,2);  ham.set_hvib_dia_by_ref(Hvibdia);  ham.set_hvib_adi_by_ref(Hvibadi);
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
    params["model"] = model
    params["rep"] = rep


    # Simulation parameters
    dt = 1.0

    # Dynamical variables and system-specific properties
    Cadi.set(0, 1.0+0.0j); Cadi.set(1, 1.0+0.0j);   Cadi *= (1.0/math.sqrt(2.0))  

    q = MATRIX(1,1); q.set(0, 0.1)
    p = MATRIX(1,1); p.set(0, 0.0)
    iM = MATRIX(1,1); iM.set(0, 1.0/100.0)


    # Initial calculations
    ham.compute_diabatic(compute_model, q, params)
    ham.compute_adiabatic(1); 
    ham.ampl_adi2dia()
    if rep==0:
        ham.compute_nac_dia(p, iM);  
        ham.compute_hvib_dia(); 
    elif rep==1:
        ham.compute_nac_adi(p, iM);
        ham.compute_hvib_adi(); 



    Etot = compute_etot(ham, p, iM, rep)
    dm_dia, dm_adi = None, None


    out = open(outname, "w")
    out.close()


    # Do the propagation
    for i in xrange(500):

        Ehrenfest(dt, q, p, iM, ham, compute_model, params, rep)

        """
        print q.get(0), p.get(0)
        print "Hdia = "; Hdia.show_matrix()
        print "Hadi = "; Hadi.show_matrix()
        print "Sdia = "; Sdia.show_matrix()
        print "NACdia = "; NACdia.show_matrix()
        print "NACadi = "; NACadi.show_matrix()
        print "Hvibdia = "; Hvibdia.show_matrix()
        print "Hvibadi = "; Hvibadi.show_matrix()
        """

        #sys.exit(0)
     
        #=========== Properties ==========
        if rep==0:
            dm_dia = Sdia * Cdia * Cdia.H() * Sdia
            dm_adi = U.H() * dm_dia * U
           
        elif rep==1:
            dm_adi = Cadi * Cadi.H()
            su = Sdia * U
            dm_dia = su * dm_adi * su.H()

        Etot = compute_etot(ham, p, iM, rep)


        out = open(outname, "a")
        ret = (i*dt, q.get(0), p.get(0), Etot, 0.5*(p.T()*iM*p).get(0), dm_adi.get(0,0).real, dm_adi.get(1,1).real, dm_dia.get(0,0).real, dm_dia.get(1,1).real,
               Hdia.get(0,0).real, Hdia.get(1,1).real, Hadi.get(0,0).real, Hadi.get(1,1).real, dc1_adi[0].get(0,1).real
              )
        out.write("%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f\n" %  ret )
        out.close()
                                                                    

model = 1

run_test(model, 0, "_0_new.txt")
run_test(model, 1, "_1_new.txt")

        
