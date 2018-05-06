#*********************************************************************************
#* Copyright (C) 2017-2018  Brendan A. Smith, Alexey V. Akimov
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
import aux_functs

class tmp:
    pass    

def compute_model(q, params, full_id):

    model = params["model"]
    res = None

    Id = Cpp2Py(full_id)
    indx = Id[-1]

    if model==1:
        res = models_Libra.model5(q.col(indx), params)
    if model==2:
        res = models_Martens.Martens1(q.col(indx), params)
    if model==3:
        res = models_Martens.Martens2(q.col(indx), params)

#    res.rep = params["rep"]    
    return res  

def run_test0():
    """
    This example runs ETHD dyanimcs in the 2D symmetrical double well potential 
    It also has the added functionality of being able to print ensemble distirbutions
    """
    
    ndia, nadi, nnucl, ntraj = 1, 1, 2, 500

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

    # Set up the models and compute internal variables
    # Initialization
    # Model parameters 
    params = { "model":1 }

    # Simulation parameters
    dt = 1.0

    # Dynamical variables and system-specific properties
    mean_q = MATRIX(nnucl,1);   
    sigma_q = MATRIX(nnucl,1);  
    mean_p = MATRIX(nnucl,1);   
    sigma_p = MATRIX(nnucl,1);  
    iM = MATRIX(nnucl,1);

    for i in xrange(nnucl):
        mean_q.set(i,0, -1.0)  
        sigma_q.set(i,0, 0.05) 
        mean_p.set(i,0, 0.0)  
        sigma_p.set(i,0, 0.0)
        iM.set(i,0, 1.0/2000.0)

    rnd = Random()
    q = MATRIX(nnucl,ntraj);  aux_functs.sample(q, mean_q, sigma_q, rnd)
    p = MATRIX(nnucl,ntraj);  aux_functs.sample(p, mean_p, sigma_p, rnd)   

    # Initial calculations
    q.show_matrix()

    # Compute Initial trajectory probability distributions for all dof
    #bin(q, -2.0, 2.0, 0.01)

    ham.compute_diabatic(compute_model, q, params, 1)
    ham.compute_adiabatic(1, 1);
    ham.add_ethd_adi(q, iM, 1)

    os.system("mkdir _2D_dist")
    out1 = open("_output.txt", "w"); out1.close()   

    # Do the propagation
    for i in xrange(100):

        aux_functs.bin2(q, -2.0, 2.0, 0.1, -2.0, 2.0, 0.1, "_2D_dist/_2D_distrib_"+str(i)+"_.txt")

        Verlet1(dt, q, p, iM, ham, compute_model, params, 1)

        #=========== Properties ==========

        Ekin, Epot, Etot = aux_functs.compute_etot(ham, p, iM)

        # Print the ensemble average - kinetic, potential, and total energies
        # Print the tunneling information. Here, we count each trajectory across the barrier.
        out1 = open("_output.txt", "a")
        out1.write( " %8.5f  %8.5f  %8.5f  %8.5f\n" % ( i*dt, Ekin, Epot, Etot ) )
        out1.close()


def run_test1(Nsnaps,Nsteps):

    """
    This example runs ETHD dyanimcs in Model 2 of Martens' paper: 
    L. Wang, C.C. Martens,and Y. Zheng, "Entangled trajectory molecular dynamics in 
    multidimensional systems: Two-dimensional quantum tunneling through 
    the Eckart barrier" J. Chem. Phys. 137, 34113 (2012)

    Like run_test0(), it too has the added functionality of being able to print ensemble 
    distirbutions
    """
    
    ndia, nadi, nnucl, ntraj = 1, 1, 2, 500

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

    # Set up the models and compute internal variables
    # Initialization
    # Model parameters 
    params = { "model":2, "ndof":nnucl, "ntraj":ntraj }

    # Simulation parameters
    dt = 1.0

    # The following function call should be uncommented
    # only if you previously ran an initial calculation calling
    # aux_functs.extract_q_p_info(q,p). 
    # aux_functs.get_q_p_info(params) retrieves the previous q and
    # p coordinates form a prevous smulation. Before uncommenting it,
    # be sure that the directory _q_p_info is present in your folder 

    #q, p = aux_functs.get_q_p_info(params)


    # Dynamical variables and system-specific properties
    mean_q = MATRIX(nnucl,1);   
    sigma_q = MATRIX(nnucl,1);
    mean_p = MATRIX(nnucl,1);   
    sigma_p = MATRIX(nnucl,1);  

    # In this example, the mean position of the position 
    # dof are not equal. The uncertainty is asscoaited with 
    # a minimal uncertainty gaussian with a mass of 2,000 (a.u)
    # and omega = 0.004, which is associated with the ground 
    # state of the harmonic oscillator
    # Ex) sigma_q = sqrt( hbar/(2*m*w) ) 
    mean_q.set(0,0, -1.0)  
    mean_q.set(1,0,  0.0)
    sigma_q.set(0,0, 0.25)
    sigma_q.set(1,0, 0.25) 

    # Calculate sigma_p from sigma_q choice. sigma_p will
    # be the minimum unertainty value based on the uncertainty
    # principle.
    mean_p.set(0,0,3.0)
    mean_p.set(1,0,0.0)
    sigma_p.set(0,0,0.5/sigma_q.get(0,0))
    sigma_p.set(1,0,0.5/sigma_q.get(1,0))
  
    rnd = Random()
    q = MATRIX(nnucl,ntraj);  tsh.sample(q, mean_q, sigma_q, rnd)
    p = MATRIX(nnucl,ntraj);  tsh.sample(p, mean_p, sigma_p, rnd)   

    # Make lists containing the positions and momenta for each trajectory,
    # to be used in later simulations, when extra trajecotry ensemble
    # replicas are needed. If running the first of many simulations in
    # which you would like  to reuse the identical ensemble, uncomment.
    aux_functs.extract_q_p_info(q,p)

    # Initial calculations
    print "\n", "Showing q_matrix", "\n"
    q.show_matrix()

    # Set mass matricies  
    iM = MATRIX(nnucl,1);
    iM.set(0,0, 1.0/2000.0)
    iM.set(1,0, 1.0/2000.0)

    # Check potential energy surface
    aux_functs.check_potential(q, params, -7.0, 7.0, 0.1, -5.0, 5.0, 0.1, "_pes_scan.txt")

    # Compute Initial trajectory probability distributions for all position dof
    # Momenta distrbutions coming soon.
    aux_functs.bin(q, -8.0, 8.0, 0.05,"_q")
    aux_functs.bin(p, -8.0, 8.0, 0.05,"_p")

    ham.compute_diabatic(compute_model, q, params, 1)
    ham.compute_adiabatic(1, 1);
    # If running an ETHD simulation, uncomment
    ham.add_ethd_adi(q, iM, 1)

    os.system("mkdir _2D_dist")
    out1 = open("_output.txt", "w"); out1.close()   
    # Do the propagation
    for i in xrange(Nsnaps):

        # Count the number of trajectories that cross the barrier 
        react_prob = aux_functs.traj_counter(q, 0.0, 0)

        aux_functs.bin2(q, -7.0, 7.0, 0.1, -5.0, 5.0, 0.1, "_2D_dist/_2D_distrib_"+str(i)+"_.txt")

        #=========== Properties ==========

        Ekin, Epot, Etot = aux_functs.compute_etot(ham, p, iM)

        # Print the ensemble average - kinetic, potential, and total energies
        # Print the tunneling information. Here, we count each trajectory across the barrier.
        out1 = open("_output.txt", "a")
        out1.write( " %8.5f %8.5f %8.5f %8.5f %8.5f\n" % ( i*dt*Nsteps, Ekin, Epot, Etot, react_prob ) )
        out1.close()
        
        for j in xrange(Nsteps):
            Verlet1(dt, q, p, iM, ham, compute_model, params, 1)

run_test1(100, 1)



