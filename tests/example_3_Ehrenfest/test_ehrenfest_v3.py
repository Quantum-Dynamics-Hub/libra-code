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

"""
  This file demonstrates how to run an ensemble of Ehrenfest trajectories
  using the hierarchy of Hamiltonians approach and the built in Ehrenfest1
  function
 
"""


import cmath
import math
import os
import sys


if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *




def compute_model(q, params, full_id):

    model = params["model"]
    res = None

    Id = Cpp2Py(full_id)
    indx = Id[-1]

    if model==1:
        res = models_Libra.model1(q.col(indx), params)
    elif model==2:
        res = models_Libra.model2(q.col(indx), params)
    elif model==3:
        res = models_Libra.model3(q.col(indx), params)
    elif model==4:
        res = models_Libra.model4(q.col(indx), params)

    res.rep = params["rep"]    

    return res
    




def run_Ehrenfest(ndia, nadi, nnucl, ntraj, _q, _p, _Cdia, _Cadi, _iM, model, rep, outname):
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
    Cdia = CMATRIX(_Cdia)
    Cadi = CMATRIX(_Cadi)

    

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
    ham.ampl_adi2dia(Cdia, Cadi, 0, 1)

    Cdia.show_matrix()
    Cadi.show_matrix()

    if rep==0:
        ham.compute_nac_dia(p, iM, 0, 1);  
        ham.compute_hvib_dia(1); 
    elif rep==1:
        ham.compute_nac_adi(p, iM, 0, 1);
        ham.compute_hvib_adi(1); 


    ham1[0].get_nac_adi().show_matrix()
    ham1[0].get_ham_adi().show_matrix()
    ham1[0].get_ham_dia().show_matrix()

    Ekin, Epot, Etot, dEkin, dEpot, dEtot = tsh.compute_etot(ham, p, Cdia, Cadi, iM, rep)


    out = open(outname, "w")
    out.close()

    # Do the propagation
    for i in xrange(500):

        if rep==0:
            Ehrenfest1(dt, q, p, iM, Cdia, ham, compute_model, params, rep)
        elif rep==1:
            Ehrenfest1(dt, q, p, iM, Cadi, ham, compute_model, params, rep)


        #=========== Properties ==========
        dm_dia, dm_adi = tsh.compute_dm(ham, Cdia, Cadi, rep, 1)

        Ekin, Epot, Etot, dEkin, dEpot, dEtot = tsh.compute_etot(ham, p, Cdia, Cadi, iM, rep)



        out = open(outname, "a")
        ret = (i*dt, q.get(0), p.get(0), 
               Ekin, Epot, Etot, dEkin, dEpot, dEtot,
               dm_adi.get(0,0).real, dm_adi.get(1,1).real, dm_dia.get(0,0).real, dm_dia.get(1,1).real
              )
        out.write("%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f\n" %  ret )
        out.close()



def run_test():                                                                    

    model = 4
    ndia, nadi, nnucl, ntraj = 2, 2, 1, 1

    rnd = Random()

    # Dynamical variables and system-specific properties
    mean_q = MATRIX(nnucl,1);   mean_q.set(0,0, 0.1)
    sigma_q = MATRIX(nnucl,1);  sigma_q.set(0,0, 0.05)
    mean_p = MATRIX(nnucl,1);   mean_p.set(0,0, 0.0)
    sigma_p = MATRIX(nnucl,1);  sigma_p.set(0,0, 0.01)

    q = MATRIX(nnucl,ntraj);  tsh.sample(q, mean_q, sigma_q, rnd)
    p = MATRIX(nnucl,ntraj);  tsh.sample(p, mean_p, sigma_p, rnd)
    iM = MATRIX(nnucl,1);     iM.set(0,0, 1.0/100.0)

    Cdia, Cadi = CMATRIX(ndia, ntraj), CMATRIX(nadi, ntraj)
    for traj in xrange(ntraj):
        Cadi.set(0, traj, 1.0+0.0j);
        Cadi.set(1, traj, 1.0+0.0j);
        Cadi *= (1.0/math.sqrt(2.0))  

    run_Ehrenfest(ndia, nadi, nnucl, ntraj, q, p, Cdia, Cadi, iM, model, 0, "_0_new.txt")
    run_Ehrenfest(ndia, nadi, nnucl, ntraj, q, p, Cdia, Cadi, iM, model, 1, "_1_new.txt")

        
run_test()
