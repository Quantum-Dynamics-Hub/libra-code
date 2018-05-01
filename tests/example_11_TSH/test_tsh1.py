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


"""
  This file demonstrates how to run an ensemble of TSH using many Hamiltonians 
"""

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *

class tmp:
    pass    


def Tully1(q, params):
    """
 
 
    """


    Hdia = CMATRIX(2,2)
    Sdia = CMATRIX(2,2)
    d1ham_dia = CMATRIXList();  d1ham_dia.append( CMATRIX(2,2) )
    dc1_dia = CMATRIXList();  dc1_dia.append( CMATRIX(2,2) )
  

    x = q.get(0)
    A = 0.01
    B = 1.6
    C = 0.005
    D = 1.0

    Sdia.set(0,0, 1.0+0.0j);  Sdia.set(0,1, 0.0+0.0j);
    Sdia.set(1,0, 0.0+0.0j);  Sdia.set(1,1, 1.0+0.0j);

    V11,dV11 = 0.0, 0.0
    if x>0:
        V11 = A*(1.0 - math.exp(-B*x))
        dV11 =  A*B*math.exp(-B*x)
    else:
        V11 = -A*(1.0 - math.exp(B*x))
        dV11 =  A*B*math.exp(B*x)
    
    V = C * math.exp(-D*x*x)
    dV = -2.0*x*C*D*math.exp(-D*x*x)

    Hdia.set(0,0, V11*(1.0+0.0j) );   Hdia.set(0,1, V*(1.0+0.0j));
    Hdia.set(1,0, V*(1.0+0.0j));      Hdia.set(1,1, V11*(-1.0+0.0j));


    for i in [0]:
        #  d Hdia / dR_0
        d1ham_dia[i].set(0,0, dV11*(1.0+0.0j) );  d1ham_dia[i].set(0,1, dV*(1.0+0.0j));
        d1ham_dia[i].set(1,0, dV*(1.0+0.0j));     d1ham_dia[i].set(1,1, dV11*(-1.0+0.0j));

        #  <dia| d/dR_0| dia >
        dc1_dia[i].set(0,0, 0.0+0.0j);   dc1_dia[i].set(0,1, 0.0+0.0j);
        dc1_dia[i].set(1,0, 0.0+0.0j);   dc1_dia[i].set(1,1, 0.0+0.0j);

    obj = tmp()
    obj.ham_dia = Hdia
    obj.ovlp_dia = Sdia
    obj.d1ham_dia = d1ham_dia
    obj.dc1_dia = dc1_dia

    return obj



def model1(q, params):
    """
              k*x^2         V
    Hdia =       V        k*(x-x0)^2 + D

    Sdia =  I

    Ddia  = 0.0

    """


    Hdia = CMATRIX(2,2)
    Sdia = CMATRIX(2,2)
    d1ham_dia = CMATRIXList();  d1ham_dia.append( CMATRIX(2,2) )
    dc1_dia = CMATRIXList();  dc1_dia.append( CMATRIX(2,2) )
  

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

#        sys.exit(0)

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




def compute_model(q, params, full_id):

    model = params["model"]
    res = None

    Id = Cpp2Py(full_id)
    indx = Id[-1]

    if model==1:
#        res = model1(q.col(indx), params)
        res = Tully1(q.col(indx), params)
    elif model==2:
        res = model2(q.col(indx), params)
    elif model==3:
        res = model3(q.col(indx), params)
    elif model==4:
        res = model4(q.col(indx), params)

    res.rep = params["rep"]    

    return res
    


def compute_etot(ham, p, iM, rep, nst, st):


    ntraj = p.num_of_cols
    ndof = p.num_of_rows

    epot, ekin = [], []    
    Epot, Ekin = 0.0, 0.0

    C = CMATRIX(nst, 1)
    C.set(st, 0, 1.0+0.0j)

    if rep==0:
        Epot = ham.Ehrenfest_energy_dia(C).real 
    elif rep==1:
        Epot = ham.Ehrenfest_energy_adi(C).real 

    Ekin = 0.0
    for dof in xrange(ndof):
        Ekin = Ekin + 0.5 * iM.get(dof, 0) * (p.get(dof, 0) ** 2)

    Etot = Ekin + Epot

    # Variances:
    dEkin, dEpot = 0.0, 0.0
    dEtot = dEkin + dEpot
   

    return Ekin, Epot, Etot, dEkin, dEpot, dEtot



def sample(x, mean_x, sigma_x, rnd):  
    nr, nc = x.num_of_rows, x.num_of_cols
    for i in range(nr):
        for j in range(nc):    
            x.set(i,j, mean_x.get(i,0) + sigma_x.get(i,0) * rnd.normal() )





def run_test(_q, _p, _Cadi, _iM, model, rep, outname, params1, rnd, st):
    """
    model - setup the Hamiltonian
    rep - 0 - diabatic, 1 - adiabatic
    outname - the name of the output file
    """

    # Create copies of the input dynamical variables, so we could run several run_test 
    # functions with the same input variables without worries that they will be altered
    # inside of run_test

    nel = _Cadi.num_of_rows;

    q = MATRIX(_q)
    p = MATRIX(_p)
    iM = MATRIX(_iM)

    Cadi = CMATRIX(_Cadi)
    Cdia = CMATRIX(nel, nel)


    # ======= Hierarchy of Hamiltonians =======
    ham = nHamiltonian(nel, nel, 1)
    ham.init_all(2)


    #  Set up the models and compute internal variables
    # Initialization
    # Model parameters 
    params = {"model":model, "rep":rep}

    # Simulation parameters
    dt = 1.0

    # Initial calculations
    ham.compute_diabatic(compute_model, q, params)
    ham.compute_adiabatic(1); 

    ham.ampl_adi2dia(Cdia, Cadi)


    if rep==0:
        ham.compute_nac_dia(p, iM);  
        ham.compute_hvib_dia(); 
    elif rep==1:
        ham.compute_nac_adi(p, iM);
        ham.compute_hvib_adi(); 

 


    Ekin, Epot, Etot, dEkin, dEpot, dEtot = compute_etot(ham, p, iM, rep, nel, st)
    dm_dia, dm_adi = None, None


    out = open(outname, "w")
    out.close()


    # Do the propagation
    for i in xrange(1000):

#        sys.exit(0)
        if rep==0:
            st = tsh0(dt, q, p, iM,  Cdia, st, ham, compute_model, params, params1, rnd)
        elif rep==1:
            st = tsh0(dt, q, p, iM,  Cadi, st, ham, compute_model, params, params1, rnd)


        #=========== Properties ==========
     
        if rep==0:
            ham.ampl_dia2adi(Cdia, Cadi)

            dm_dia = ham.get_ovlp_dia() * Cdia * Cdia.H() * ham.get_ovlp_dia()
            dm_adi = ham.get_basis_transform().H() * dm_dia * ham.get_basis_transform()
        
        elif rep==1:
            ham.ampl_adi2dia(Cdia, Cadi)
            dm_adi = Cadi * Cadi.H()
            su = ham.get_ovlp_dia() * ham.get_basis_transform()
            dm_dia = su * dm_adi * su.H()


        Ekin, Epot, Etot, dEkin, dEpot, dEtot = compute_etot(ham, p, iM, rep, nel, st)


        out = open(outname, "a")
        """
        ret = (i*dt, q.get(0), p.get(0), Etot, Ekin, dm_adi.get(0,0).real, dm_adi.get(1,1).real, dm_dia.get(0,0).real, dm_dia.get(1,1).real,
               Hdia.get(0,0).real, Hdia.get(1,1).real, Hadi.get(0,0).real, Hadi.get(1,1).real, dc1_adi[0].get(0,1).real
              )
        out.write("%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f\n" %  ret )
        """

        ret = (i*dt, q.get(0), p.get(0), 
               Ekin, Epot, Etot, dEkin, dEpot, dEtot,
               dm_adi.get(0,0).real, dm_adi.get(1,1).real, dm_dia.get(0,0).real, dm_dia.get(1,1).real, st
              )
        out.write("%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %5i\n" %  ret )


        out.close()


                                                                    


model = 1
rnd = Random()

q = MATRIX(1,1);   q.set(0, 0, -5.0) 
p = MATRIX(1,1);   p.set(0, 0, 20.0) 
iM = MATRIX(1,1);  iM.set(0,0, 1.0/2000.0)

Cdia, Cadi = CMATRIX(2,1), CMATRIX(2,1)

rep = 1  # representation for nuclear dynamics
istate = 0;  Cadi.set(istate, 0, 1.0+0.0j);




params1 = {"rep":1, "rep_sh":1, "tsh_method":0, "use_boltz_factor":0, 
           "Temperature":300.0, "do_reverse":1, "vel_rescale_opt":0 }

run_test(q, p, Cadi, iM, model, rep, "_0_new.txt", params1, rnd, istate)

        
