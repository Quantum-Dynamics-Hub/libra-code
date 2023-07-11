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
from pyscf import gto, scf
from pyscf import ci, grad


class tmp:
    pass    


def model_cisd(q, params, full_id):
    """

    """

    Id = Cpp2Py(full_id)
    indx = Id[-1]

    r = q.col(indx) 
    x1 = r.get(0); y1 = r.get(1); z1 = r.get(2)
    x2 = r.get(3); y2 = r.get(4); z2 = r.get(5)
    r = (1.0/1.889725989) * r  # convert to Angst

    mol = gto.M(atom=[['H', (x1,y1,z1)], ['H', (x2,y2,z2)]],  basis="631g")
    mf = scf.RHF(mol)
    mf.verbose=-1
    res_mf = mf.kernel()
    print "HF energy = ", res_mf
#   force = mf.nuc_grad_method().kernel()
#    force = mf.grad()

    force = grad.RHF(mf).kernel()
    print force

    """
    myci = ci.cisd.CISD(mf)
    myci.verbose = -1
    myci.nstates = 1
    res_ci = myci.kernel()
    print "CISD correlation energy = ", res_ci[0]
    print "CISD vectors = ", res_ci[1]

    force = myci.nuc_grad_method().kernel(state=0)
    print force
    """


    obj = tmp()
    obj.ham_dia = CMATRIX(1,1)
    obj.ovlp_dia = CMATRIX(1,1)
    d1ham_dia = CMATRIXList(); 
    dc1_dia = CMATRIXList(); 

    for i in xrange(6):
        d1ham_dia.append( CMATRIX(1,1) )
        dc1_dia.append( CMATRIX(1,1) )
  

    obj.ham_dia.set(0,0, (res_mf)*(1.0+0.0j) )
#    obj.ham_dia.set(0,0, (res_mf + res_ci[0])*(1.0+0.0j) )
    obj.ovlp_dia.set(0,0, 1.0+0.0j)

    #  d Hdia / dR_0
    d1ham_dia[0].set(0,0, (force[0][0])*(1.0+0.0j) )
    d1ham_dia[1].set(0,0, (force[0][1])*(1.0+0.0j) )
    d1ham_dia[2].set(0,0, (force[0][2])*(1.0+0.0j) )
    d1ham_dia[3].set(0,0, (force[1][0])*(1.0+0.0j) )
    d1ham_dia[4].set(0,0, (force[1][1])*(1.0+0.0j) )
    d1ham_dia[5].set(0,0, (force[1][2])*(1.0+0.0j) )


    for i in xrange(6):
        #  <dia| d/dR_0| dia >
        dc1_dia[i].set(0,0, 0.0+0.0j)

    obj.d1ham_dia = d1ham_dia
    obj.dc1_dia = dc1_dia
    
    return obj
    



def compute_model(q, params, full_id):

    res = None

    model = 1 
    if model==1:
        res = model_cisd(q, params, full_id)

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



def run_test(model, outname, output_1, output_2, output_3):
    """
    model - setup the Hamiltonian
    outname - the name of the output file
    """

    ndia, nadi, nnucl, ntraj = 1, 1, 6, 3

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

    # Simulation parameters
    dt = 10.0

    # Dynamical variables and system-specific properties
    mean_q = MATRIX(nnucl,1);   
    mean_q.set(0,0, 1.2)
    mean_q.set(1,0, 0.0)
    mean_q.set(2,0, 0.0)
    mean_q.set(3,0, 0.0)
    mean_q.set(4,0, 0.0)
    mean_q.set(5,0, 0.0)

    iM = MATRIX(nnucl,1);
    for i in xrange(6):
        iM.set(i,0, 1.0/1836.0)

    for i in xrange(6):
        sigma_q = MATRIX(nnucl,1);  sigma_q.set(i,0, 0.05)
        mean_p = MATRIX(nnucl,1);   mean_p.set(i,0, 0.0)
        sigma_p = MATRIX(nnucl,1);  sigma_p.set(i,0, 0.01)

    mean_q = mean_q * 1.889725989
    sigma_q = sigma_q * 1.889725989

    rnd = Random()
    q = MATRIX(nnucl,ntraj);  sample(q, mean_q, sigma_q, rnd)
    p = MATRIX(nnucl,ntraj);  sample(p, mean_p, sigma_p, rnd)



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

    out1 = open("_output.txt", "w"); out1.close()   
    out2 = open("_phase_space.txt", "w"); out2.close()  
    out3 = open("_pos_space.txt", "w"); out3.close()  

    # Do the propagation
    for i in xrange(250):

        Verlet1(dt, q, p, iM, ham, compute_model, params, 0)
     
        #=========== Properties ==========

        Ekin, Epot, Etot = compute_etot(ham, p, iM)

        # Print the ensemble average - kinetic, potential, and total energies
        # Print the tunneling information. Here, we count each trajectory across the barrier.
        out1 = open(output_1, "a")
        out1.write( " %8.5f  %8.5f  %8.5f  %8.5f\n" % ( i*dt, Ekin, Epot, Etot ) )
        out1.close()

        # Print the phase space information
        out2 = open(output_2, "a") 
        for j in range(ntraj):
            out2.write(" %8.5f  %8.5f" % ( q.get(0,j), p.get(0,j) ) ),
        out2.write( "\n" )
        out2.close()  

        # Print the position versus time infromation
        out3 = open(output_3, "a")
        out3.write( " %8.5f" % (i*dt) ),
        for j in range(ntraj):
            out3.write( " %8.5f" % ( q.get(0,j) ) )
        out3.write( "\n" )
        out3.close()  


model = 1
run_test(model, "_0_new.txt", "_output.txt", "_phase_space.txt", "_pos_space.txt")

        
