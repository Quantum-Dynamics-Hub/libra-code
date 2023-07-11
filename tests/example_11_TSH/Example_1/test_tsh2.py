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
  focus is on the Tully's models
"""

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *

#import energy

class tmp:
    pass    



def compute_model(q, params, full_id):

    model = params["model"]
    res = None

    Id = Cpp2Py(full_id)
    indx = Id[-1]

    if model==1:
        res = models_Tully.Tully1(q.col(indx), params)
    elif model==2:
        pass
    elif model==3:
        pass

    res.rep = params["rep"]    

    return res
    




def get_probabilities(ham, states):
    """
    states = [list of ntraj integers] matrix with states
    """

    nst = ham.nadi
    ntraj = len(states)

    res = MATRIX(nst, 1)

    for traj in xrange(ntraj):        
        res.add(states[traj], 0, 1.0)

    res = res * (1.0/float(ntraj))

    return res



def run_test(ndia, nadi, nnucl, ntraj, _q, _p, _Cdia, _Cadi, _iM, model, rep, outname, params1, rnd, states):
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
    params = {"model":model, "rep":rep}

    # Simulation parameters
    dt = 0.5


    # Initial calculations
    ham.compute_diabatic(compute_model, q, params, 1)
    ham.compute_adiabatic(1, 1); 
    ham.ampl_adi2dia(Cdia, Cadi, 0, 1)


    if rep==0:
        ham.compute_nac_dia(p, iM, 0, 1);  
        ham.compute_hvib_dia(1); 
    elif rep==1:
        ham.compute_nac_adi(p, iM, 0, 1);
        ham.compute_hvib_adi(1); 


#    sys.exit(0)
    Ekin, Epot, Etot, dEkin, dEpot, dEtot = tsh.compute_etot_tsh(ham, p, Cdia, Cadi, states, iM, rep) 
    print Ekin, Epot, Etot, dEkin, dEpot, dEtot



    out = open(outname, "w")
    out.close()


    # Do the propagation
    for i in xrange(400):


        if rep==0:
            tsh1(dt, q, p, iM,  Cdia, states, ham, compute_model, params, params1, rnd)
        elif rep==1:
            tsh1(dt, q, p, iM,  Cadi, states, ham, compute_model, params, params1, rnd, 1, 1)


        #=========== Properties ==========
        if rep==0:
            ham.ampl_dia2adi(Cdia, Cadi, 0, 1)
        elif rep==1:
            ham.ampl_adi2dia(Cdia, Cadi, 0, 1)


        dm_dia, dm_adi = tsh.compute_dm(ham, Cdia, Cadi, rep, 1)
        Ekin, Epot, Etot, dEkin, dEpot, dEtot = tsh.compute_etot_tsh(ham, p, Cdia, Cadi, states, iM, rep)
        pops = get_probabilities(ham, states)

        out = open(outname, "a")


        ind = 0.0
        for tr in xrange(ntraj):
            ind = ind + ham1[tr].get_ordering_adi()[0]
        ind = ind/float(ntraj)

        ret = (i*dt, q.get(0), p.get(0), 
               Ekin, Epot, Etot, dEkin, dEpot, dEtot,
               dm_adi.get(0,0).real, dm_adi.get(1,1).real, dm_dia.get(0,0).real, dm_dia.get(1,1).real, pops.get(0,0), pops.get(1,0), ind
              )
        out.write("%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f\n" %  ret )

        out.close()


                                                                    


model = 1
ndia, nadi, nnucl, ntraj = 2, 2, 1, 25

rnd = Random()

# Dynamical variables and system-specific properties
mean_q = MATRIX(nnucl,1);   mean_q.set(0,0, -2)
sigma_q = MATRIX(nnucl,1);  sigma_q.set(0,0, 0.0)
mean_p = MATRIX(nnucl,1);   mean_p.set(0,0, 20.0)
sigma_p = MATRIX(nnucl,1);  sigma_p.set(0,0, 0.0)

q = MATRIX(nnucl,ntraj);  tsh.sample(q, mean_q, sigma_q, rnd)
p = MATRIX(nnucl,ntraj);  tsh.sample(p, mean_p, sigma_p, rnd)
iM = MATRIX(nnucl,1);     iM.set(0,0, 1.0/2000.0)

istate = 0
Cdia, Cadi = CMATRIX(ndia, ntraj), CMATRIX(nadi, ntraj)
states = intList() #CMATRIX(ndia, ntraj)

for traj in xrange(ntraj):
    Cadi.set(istate, traj, 1.0+0.0j);  
    states.append(istate) #set(istate, traj, 1.0+0.0j)


rep = 1
params1 = {"rep":rep, "rep_sh":1, "tsh_method":0, "use_boltz_factor":0,
           "Temperature":300.0, "do_reverse":1, "vel_rescale_opt":0 }
run_test(ndia, nadi, nnucl, ntraj, q, p, Cdia, Cadi, iM, model, rep, "_0_new.txt", params1, rnd, states)



#rep = 0
#params1 = {"rep":rep, "rep_sh":0, "tsh_method":0, "use_boltz_factor":0,
#           "Temperature":300.0, "do_reverse":1, "vel_rescale_opt":0 }
#run_test(ndia, nadi, nnucl, ntraj, q, p, Cdia, Cadi, iM, model, rep, "_0_new.txt", params1, rnd, states)


        
