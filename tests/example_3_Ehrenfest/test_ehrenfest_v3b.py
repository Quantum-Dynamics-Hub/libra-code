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
  This file demonstrates how to run an ensemble of Ehrenfest trajectories
  using the hierarchy of Hamiltonians approach.
 
"""

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
        res = model1(q.col(indx), params)
    elif model==2:
        res = model2(q.col(indx), params)
    elif model==3:
        res = model3(q.col(indx), params)
    elif model==4:
        res = model4(q.col(indx), params)

    res.rep = params["rep"]    

    return res
    


def compute_etot(ham, p, Cdia, Cadi, iM, rep):


    ntraj = p.num_of_cols
    ndof = p.num_of_rows

    epot, ekin = [], []    
    Epot, Ekin = 0.0, 0.0
    for traj in xrange(ntraj):
        if rep==0:
            epot.append( ham.Ehrenfest_energy_dia(Cdia[traj], Py2Cpp_int([0,traj])).real )
            Epot = Epot + epot[traj]
        elif rep==1:
            epot.append( ham.Ehrenfest_energy_adi(Cdia[traj], Py2Cpp_int([0,traj])).real )
            Epot = Epot + epot[traj]

        tmp = 0.0
        for dof in xrange(ndof):
            tmp = tmp + 0.5 * iM.get(dof, 0) * (p.get(dof, traj) ** 2)
        ekin.append(tmp)
        Ekin = Ekin + ekin[traj]

    Ekin = Ekin / float(ntraj)
    Epot = Epot / float(ntraj)
    Etot = Ekin + Epot

    # Variances:
    dEkin, dEpot = 0.0, 0.0
    for traj in xrange(ntraj):
        dEkin = dEkin + (ekin[traj] - Ekin)**2
        dEpot = dEpot + (epot[traj] - Epot)**2

    dEtot = dEkin + dEpot

    dEkin = math.sqrt(dEkin/ float(ntraj))
    dEpot = math.sqrt(dEpot/ float(ntraj))
    dEtot = math.sqrt(dEtot/ float(ntraj))
    


    return Ekin, Epot, Etot, dEkin, dEpot, dEtot



def sample(x, mean_x, sigma_x, rnd):  
    nr, nc = x.num_of_rows, x.num_of_cols
    for i in range(nr):
        for j in range(nc):    
            x.set(i,j, mean_x.get(i,0) + sigma_x.get(i,0) * rnd.normal() )





def run_test(ndia, nadi, nnucl, ntraj, _q, _p, _Cdia, _Cadi, _iM, model, rep, outname):
    """
    model - setup the Hamiltonian
    rep - 0 - diabatic, 1 - adiabatic
    outname - the name of the output file
    """

    # Create copies of the input dynamical variables, so we could run several run_test 
    # functions with the same input variables without worries that they will be altered
    # inside of run_test

    q = MATRIX(_q)
    p = MATRIX(_p)
    iM = MATRIX(_iM)

    Cdia = CMATRIXList()
    for c in _Cdia:
        Cdia.append(CMATRIX(c))

    Cadi = CMATRIXList()
    for c in _Cadi:
        Cadi.append(CMATRIX(c))

    

    # ======= Hierarchy of Hamiltonians =======
    ham = nHamiltonian(ndia, nadi, nnucl)
    ham.init_all(2)
    print "id=", ham.id, " level=", ham.level

    ham1 = [] 
    for tr in xrange(ntraj):
        ham1.append( nHamiltonian(ndia, nadi, nnucl) )        
        print ham1[tr].id, ham1[tr].level
        ham1[tr].init_all(2)
        ham.add_child(ham1[tr])
        print Cpp2Py(ham1[tr].get_full_id())



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


    # Initial calculations
    ham.compute_diabatic(compute_model, q, params, 1)
    ham.compute_adiabatic(1, 1); 


    for traj in xrange(ntraj):
        indx = Py2Cpp_int([0,traj])
        ham.ampl_adi2dia(Cdia[traj], Cadi[traj],  indx )

        if rep==0:
            ham.compute_nac_dia(p, iM, indx);  
            ham.compute_hvib_dia(indx); 
        elif rep==1:
            ham.compute_nac_adi(p, iM, indx);
            ham.compute_hvib_adi(indx); 


    Ekin, Epot, Etot, dEkin, dEpot, dEtot = compute_etot(ham, p, Cdia, Cadi, iM, rep)
    dm_dia, dm_adi = None, None


    out = open(outname, "w")
    out.close()


    # Do the propagation
    for i in xrange(500):

        if rep==0:
            Ehrenfest1(dt, q, p, iM, Cdia, ham, compute_model, params, 0)
        elif rep==1:
            Ehrenfest1(dt, q, p, iM, Cadi, ham, compute_model, params, 1)



        #=========== Properties ==========
        for traj in xrange(ntraj):
            indx = Py2Cpp_int([0,traj])
     
            if rep==0:
                ham.ampl_dia2adi(Cdia[traj], Cadi[traj], indx)

                dm_dia = ham.get_ovlp_dia(indx) * Cdia[traj] * Cdia[traj].H() * ham.get_ovlp_dia(indx)
                dm_adi = ham.get_basis_transform(indx).H() * dm_dia * ham.get_basis_transform(indx)
           
            elif rep==1:
                ham.ampl_adi2dia(Cdia[traj], Cadi[traj], indx)
                dm_adi = Cadi[traj] * Cadi[traj].H()
                su = ham.get_ovlp_dia(indx) * ham.get_basis_transform(indx)
                dm_dia = su * dm_adi * su.H()


        Ekin, Epot, Etot, dEkin, dEpot, dEtot = compute_etot(ham, p, Cdia, Cadi, iM, rep)


        out = open(outname, "a")
        """
        ret = (i*dt, q.get(0), p.get(0), Etot, Ekin, dm_adi.get(0,0).real, dm_adi.get(1,1).real, dm_dia.get(0,0).real, dm_dia.get(1,1).real,
               Hdia.get(0,0).real, Hdia.get(1,1).real, Hadi.get(0,0).real, Hadi.get(1,1).real, dc1_adi[0].get(0,1).real
              )
        out.write("%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f\n" %  ret )
        """

        ret = (i*dt, q.get(0), p.get(0), 
               Ekin, Epot, Etot, dEkin, dEpot, dEtot,
               dm_adi.get(0,0).real, dm_adi.get(1,1).real, dm_dia.get(0,0).real, dm_dia.get(1,1).real
              )
        out.write("%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f\n" %  ret )


        out.close()


                                                                    


model = 1
ndia, nadi, nnucl, ntraj = 2, 2, 1, 1

rnd = Random()

# Dynamical variables and system-specific properties
mean_q = MATRIX(nnucl,1);   mean_q.set(0,0, 0.1)
sigma_q = MATRIX(nnucl,1);  sigma_q.set(0,0, 0.05)
mean_p = MATRIX(nnucl,1);   mean_p.set(0,0, 0.0)
sigma_p = MATRIX(nnucl,1);  sigma_p.set(0,0, 0.01)

q = MATRIX(nnucl,ntraj);  sample(q, mean_q, sigma_q, rnd)
p = MATRIX(nnucl,ntraj);  sample(p, mean_p, sigma_p, rnd)
iM = MATRIX(nnucl,1);     iM.set(0,0, 1.0/100.0)

Cdia, Cadi = CMATRIXList(), CMATRIXList()
for traj in xrange(ntraj):
    Cdia.append( CMATRIX(ndia,1) );    Cadi.append( CMATRIX(nadi,1) )                    
    Cadi[traj].set(0, 1.0+0.0j);    Cadi[traj].set(1, 1.0+0.0j);    Cadi[traj] *= (1.0/math.sqrt(2.0))  


run_test(ndia, nadi, nnucl, ntraj, q, p, Cdia, Cadi, iM, model, 0, "_0_new.txt")
run_test(ndia, nadi, nnucl, ntraj, q, p, Cdia, Cadi, iM, model, 1, "_1_new.txt")

        
