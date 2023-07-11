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
    



def get_probabilities(q, states, params, ndia, nadi, rep, ham, Cdia, Cadi):
    """
    Computes the scattering probabilities
    q = [ndof x ntraj] coordinates
    states = [list of ntraj integers]
    params = defines conditions
    """

    nst = ham.nadi
    ntraj = len(states)

    pop_transm = MATRIX(nst, 1)  # transmitted
    pop_refl = MATRIX(nst, 1)  # reflected

    act_dof = params["act_dof"]
    left_boundary = params["left_boundary"]
    right_boundary = params["right_boundary"]


    ntransm, nrefl = 0.0, 0.0
    for traj in xrange(ntraj):

        if q.get(act_dof, traj) < left_boundary:
            pop_refl.add(states[traj], 0, 1.0)
            nrefl += 1.0

        if q.get(act_dof, traj) > right_boundary:
            pop_transm.add(states[traj], 0, 1.0)
            ntransm += 1.0         
 
    if ntransm > 0.0:
        pop_transm = pop_transm / ntransm

    if nrefl > 0.0:
        pop_refl = pop_refl / nrefl



    dm_dia_tr, dm_adi_tr = CMATRIX(ndia, ndia), CMATRIX(nadi, nadi)
    dm_dia_re, dm_adi_re = CMATRIX(ndia, ndia), CMATRIX(nadi, nadi)

    for traj in xrange(ntraj):
        indx = Py2Cpp_int([0,traj])
    
        if rep==0:
            dm_tmp = ham.get_ovlp_dia(indx) * Cdia.col(traj) * Cdia.col(traj).H() * ham.get_ovlp_dia(indx)

            if q.get(act_dof, traj) < left_boundary:
                dm_dia_re = dm_dia_re + dm_tmp
                dm_adi_re = dm_adi_re + ham.get_basis_transform(indx).H() * dm_tmp * ham.get_basis_transform(indx)

            if q.get(act_dof, traj) > right_boundary:
                dm_dia_tr = dm_dia_tr + dm_tmp
                dm_adi_tr = dm_adi_tr + ham.get_basis_transform(indx).H() * dm_tmp * ham.get_basis_transform(indx)

       
        elif rep==1:
            c = Cadi.col(traj)
            M = ham.get_ordering_adi(Py2Cpp_int([0, traj]))
            iM = inverse_permutation(M)
            c.permute_rows(iM)
            dm_tmp = c * c.H()
            S = ham.get_ovlp_dia(indx)
            U = ham.get_basis_transform(indx)     
#            correct_phase(U)
            su = S * U


            su = ham.get_ovlp_dia(indx) * ham.get_basis_transform(indx)

            if q.get(act_dof, traj) < left_boundary:
                dm_adi_re = dm_adi_re + dm_tmp
                dm_dia_re = dm_dia_re + su * dm_tmp * su.H()

            if q.get(act_dof, traj) > right_boundary:
                dm_adi_tr = dm_adi_tr + dm_tmp
                dm_dia_tr = dm_dia_tr + su * dm_tmp * su.H()

        
    dm_adi_tr = dm_adi_tr / float(ntransm)
    dm_dia_tr = dm_dia_tr / float(ntransm)
    dm_adi_re = dm_adi_re / float(nrefl)
    dm_dia_re = dm_dia_re / float(nrefl)


    return pop_refl, pop_transm, dm_adi_re.real(), dm_dia_re.real(), dm_adi_tr.real(), dm_dia_tr.real()



def run_test(dt, md_run, ndia, nadi, nnucl, ntraj, _q, _p, _Cdia, _Cadi, _iM, model, rep, params1, rnd, _states):
    """
    model - setup the Hamiltonian
    rep - 0 - diabatic, 1 - adiabatic
    """

    # Create copies of the input dynamical variables, so we could run several run_test 
    # functions with the same input variables without worries that they will be altered
    # inside of run_test

    q = MATRIX(_q)
    p = MATRIX(_p)
    iM = MATRIX(_iM)
    Cdia = CMATRIX(_Cdia)
    Cadi = CMATRIX(_Cadi)
    states = intList()
    for tr in xrange(ntraj):
        states.append(_states[tr])
    

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


    Ekin, Epot, Etot, dEkin, dEpot, dEtot = tsh.compute_etot_tsh(ham, p, Cdia, Cadi, states, iM, rep) 
    print Ekin, Epot, Etot, dEkin, dEpot, dEtot


    # Do the propagation
    for i in xrange(md_run):

        if rep==0:
            tsh1(dt, q, p, iM,  Cdia, states, ham, compute_model, params, params1, rnd)
        elif rep==1:
            tsh1(dt, q, p, iM,  Cadi, states, ham, compute_model, params, params1, rnd, 1, 1)
 

        #=========== Properties ==========
        if rep==0:
            ham.ampl_dia2adi(Cdia, Cadi, 0, 1)
        elif rep==1:
            ham.ampl_adi2dia(Cdia, Cadi, 0, 1)


    Ekin, Epot, Etot, dEkin, dEpot, dEtot = tsh.compute_etot_tsh(ham, p, Cdia, Cadi, states, iM, rep) 
    params_observ = {"act_dof":0, "left_boundary": -10.0, "right_boundary": 5.0 }

    return get_probabilities(q, states, params_observ, ndia, nadi, rep, ham, Cdia, Cadi)
                                                                    


def run_scattering():

    model = 1
    ndia, nadi, nnucl, ntraj = 2, 2, 1, 100

    rnd = Random()

    # Dynamical variables and system-specific properties
    mean_q = MATRIX(nnucl,1);   mean_q.set(0,0, -5)
    sigma_q = MATRIX(nnucl,1);  sigma_q.set(0,0, 0.0)
    mean_p = MATRIX(nnucl,1);   mean_p.set(0,0, 20.0)
    sigma_p = MATRIX(nnucl,1);  sigma_p.set(0,0, 0.0)

    q = MATRIX(nnucl,ntraj);  tsh.sample(q, mean_q, sigma_q, rnd)
    iM = MATRIX(nnucl,1);     iM.set(0,0, 1.0/2000.0)

    rep = 1
    istate = 0
    Cdia, Cadi = CMATRIX(ndia, ntraj), CMATRIX(nadi, ntraj)
    states = intList()
    


    for traj in xrange(ntraj):
        Cadi.set(istate, traj, 1.0+0.0j);  
        states.append(istate)

    params1 = {"rep":rep, "rep_sh":1, "tsh_method":0, "use_boltz_factor":0,
               "Temperature":300.0, "do_reverse":1, "vel_rescale_opt":0 }


    out = open("_scattering.txt", "w")
    out.close()

    p0 = 1.0
    while p0<50.0:

        mean_p.set(0, 0, p0)
        p = MATRIX(nnucl,ntraj);  tsh.sample(p, mean_p, sigma_p, rnd)

        dt = 10.0 / p0
        md_run = 5000

        pops = run_test(dt, md_run, ndia, nadi, nnucl, ntraj, q, p, Cdia, Cadi, iM, model, rep, params1, rnd, states)
 
        out = open("_scattering.txt", "a")
        rec = (p0, pops[0].get(0), pops[0].get(1),  pops[1].get(0), pops[1].get(1),
              pops[2].get(0,0), pops[2].get(1,1), pops[3].get(0,0), pops[3].get(1,1),
              pops[4].get(0,0), pops[4].get(1,1), pops[5].get(0,0), pops[5].get(1,1))
        out.write("%8.5f Refl(0)= %8.5f Refl(1)= %8.5f Trans(0)= %8.5f Trans(1)= %8.5f \
                         Refl_adi(0)= %8.5f Refl_adi(1)= %8.5f  Refl_dia(0)= %8.5f Refl_dia(1)= %8.5f \
                         Trans_adi(0)= %8.5f Trans_adi(1)= %8.5f  Trans_dia(0)= %8.5f Trans_dia(1)= %8.5f\n" %  rec)
        out.close()

        p0 = p0 + 0.5
        

run_scattering()
