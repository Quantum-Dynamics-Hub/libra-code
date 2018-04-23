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


def model1(q, params, full_id):
    """
    Hdia = 0.5*k*x^2   
    Sdia =  1.0
    Ddia  = 0.0

    """

    Id = Cpp2Py(full_id)
    indx = Id[-1]

    Hdia = CMATRIX(1,1)
    Sdia = CMATRIX(1,1)
    d1ham_dia = CMATRIXList();  d1ham_dia.append( CMATRIX(1,1) )
    dc1_dia = CMATRIXList();  dc1_dia.append( CMATRIX(1,1) )
  
    x = q.col(indx).get(0)
    x0,k,D,V = params["x0"], params["k"], params["D"], params["V"]

    Hdia.set(0,0, k*x*x*(1.0+0.0j) )
    Sdia.set(0,0, 1.0+0.0j)

    for i in [0]:
        #  d Hdia / dR_0
        d1ham_dia[i].set(0,0, 2.0*k*x*(1.0+0.0j) )

        #  <dia| d/dR_0| dia >
        dc1_dia[i].set(0,0, 0.0+0.0j)


    obj = tmp()
    obj.ham_dia = Hdia
    obj.ovlp_dia = Sdia
    obj.d1ham_dia = d1ham_dia
    obj.dc1_dia = dc1_dia

    return obj
    



def compute_model(q, params, full_id):

    model = params["model"]
    res = None

    if model==1:
        res = model1(q, params, full_id)

#    res.rep = params["rep"]    
    return res
    


def compute_etot(ham, p, iM):

    ntraj = p.num_of_cols
    ndof = p.num_of_rows

    Epot, Ekin = 0.0, 0.0
    for traj in xrange(ntraj):
        Epot = Epot + ham.get_ham_adi(Py2Cpp_int([0,traj])).get(0,0).real

        for dof in xrange(ndof):
            Ekin = Ekin + 0.5 * iM.get(dof, 0) * (p.get(dof, traj) ** 2)

    Ekin = Ekin / float(ntraj) 
    Epot = Epot / float(ntraj) 
    Etot = Ekin + Epot
    
    return Ekin, Epot, Etot


def sample(x, mean_x, sigma_x, rnd):  
    nr, nc = x.num_of_rows, x.num_of_cols
    for i in range(nr):
        for j in range(nc):    
            x.set(i,j, mean_x.get(i,0) + sigma_x.get(i,0) * rnd.normal() )



def run_test(model, outname):
    """
    model - setup the Hamiltonian
    outname - the name of the output file
    """

    ndia, nadi, nnucl, ntraj = 1, 1, 1, 5

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
#    params["rep"] = rep


    # Simulation parameters
    dt = 1.0

    # Dynamical variables and system-specific properties
    mean_q = MATRIX(nnucl,1);   mean_q.set(0,0, 0.1)
    sigma_q = MATRIX(nnucl,1);  sigma_q.set(0,0, 0.05)
    mean_p = MATRIX(nnucl,1);   mean_p.set(0,0, 0.0)
    sigma_p = MATRIX(nnucl,1);  sigma_p.set(0,0, 0.01)

    rnd = Random()
    q = MATRIX(nnucl,ntraj);  sample(q, mean_q, sigma_q, rnd)
    p = MATRIX(nnucl,ntraj);  sample(p, mean_p, sigma_p, rnd)
    iM = MATRIX(nnucl,1);     iM.set(0,0, 1.0/100.0)


    # Initial calculations
    q.show_matrix()
    print "=== compute dia =="
#    ham.compute_diabatic(compute_model, q, params, 0)
    ham.compute_diabatic(compute_model, q, params, 1)

    print "=== compute adi =="
    ham.compute_adiabatic(1, 1); 

    for tr in xrange(ntraj):
        print Cpp2Py(ham1[tr].get_full_id())
        print "id = ", ham1[tr].id

    for i in xrange(ntraj):
        print i, "diabatic:", ham.get_ham_dia(Py2Cpp_int([0,i])).show_matrix()
        print i, "adiabatic:", ham.get_ham_adi(Py2Cpp_int([0,i])).show_matrix()



    Ekin, Epot, Etot = compute_etot(ham, p, iM)
    print Ekin, Epot, Etot

    out = open(outname, "w")
    out.close()


    # Do the propagation
    for i in xrange(250):

        Verlet1(dt, q, p, iM, ham, compute_model, params)
     
        #=========== Properties ==========

        Ekin, Epot, Etot = compute_etot(ham, p, iM)

        out = open(outname, "a")
        ret = (i*dt, q.get(0,1), p.get(0,1), Ekin, Epot, Etot )
        out.write("%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f \n" %  ret )
        out.close()
                                                                    

model = 1

run_test(model, "_0_new.txt")

        
