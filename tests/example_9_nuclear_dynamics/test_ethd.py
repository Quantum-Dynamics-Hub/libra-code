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
    
def model2(q, params, full_id):
    """
    Symmetric Double Well Potential
    Hdia = 0.25*q^4 - 0.5*q^2   
    Sdia = 1.0
    Ddia = 0.0
    """

    Id = Cpp2Py(full_id)
    indx = Id[-1]

    Hdia = CMATRIX(1,1)
    Sdia = CMATRIX(1,1)
    d1ham_dia = CMATRIXList();  d1ham_dia.append( CMATRIX(1,1) )
    dc1_dia = CMATRIXList();  dc1_dia.append( CMATRIX(1,1) )

    x = q.col(indx).get(0)

    Hdia.set(0,0, (0.25*x*x*x*x - 0.5*x*x)*(1.0+0.0j))
    Sdia.set(0,0, 1.0+0.0j)

    for i in [0]:
        #  d Hdia / dR_0
        d1ham_dia[i].set(0,0, (x*x*x - x)*(1.0+0.0j))

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
    if model==2:
        res = model2(q, params, full_id)

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

    Epot += ham.get_ham_adi().get(0,0).real
    Epot /= float(ntraj)  
    Ekin /= float(ntraj)

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

    ndia, nadi, nnucl, ntraj = 1, 1, 1, 100

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
    sigma_q = MATRIX(nnucl,1);  sigma_q.set(0,0, 1.0)
    mean_p = MATRIX(nnucl,1);   mean_p.set(0,0, 0.0)
    sigma_p = MATRIX(nnucl,1);  sigma_p.set(0,0, 0.01)

    rnd = Random()
    q = MATRIX(nnucl,ntraj);  sample(q, mean_q, sigma_q, rnd)
    p = MATRIX(nnucl,ntraj);  sample(p, mean_p, sigma_p, rnd)
    iM = MATRIX(nnucl,1);     iM.set(0,0, 1.0/2000.0)


    # Initial calculations
    q.show_matrix()
    print "=== compute dia =="
#    ham.compute_diabatic(compute_model, q, params, 0)
    ham.compute_diabatic(compute_model, q, params, 1)

#
#    ham.add_ethd_dia(q, iM, 1)
#    sys.exit(0)


    print "=== compute adi =="
    ham.compute_adiabatic(1, 1); 
    ham.add_ethd_adi(q, iM, 1)

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
    for i in xrange(100):

        Verlet1(dt, q, p, iM, ham, compute_model, params, 1)
     
        #=========== Properties ==========

        Ekin, Epot, Etot = compute_etot(ham, p, iM)

        out = open(outname, "a")
        ret = (i*dt, q.get(0,1), p.get(0,1), Ekin, Epot, Etot )
        out.write("%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f \n" %  ret )
        out.close()


def run_test1(model, outname):
    """
    model - setup the Hamiltonian
    outname - the name of the output file
    """

    ndia, nadi, nnucl, ntraj = 1, 1, 1, 100

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
    params = {}
    params["x0"], params["k"], params["D"], params["V"] = 1.0, 0.05, -0.1, 0.05
    params["omega"] = 0.25
    params["model"] = model
    #params["rep"] = rep

    # Simulation parameters
    dt = 1.0
    if model == 2:
        barrier = 0.0

    # Dynamical variables and system-specific properties
    mean_q = MATRIX(nnucl,1);   mean_q.set(0,0, -1.1)
    sigma_q = MATRIX(nnucl,1);  sigma_q.set(0,0, 0.04)
    mean_p = MATRIX(nnucl,1);   mean_p.set(0,0, 0.0)
    sigma_p = MATRIX(nnucl,1);  sigma_p.set(0,0, 0.0)

    rnd = Random()
    q = MATRIX(nnucl,ntraj);  #sample(q, mean_q, sigma_q, rnd)
    p = MATRIX(nnucl,ntraj);  sample(p, mean_p, sigma_p, rnd)

    qq = [-1.09854, -1.12192, -1.09621, -1.05124, -1.16032, -1.06631, -1.13264, -1.08911, -1.13068, -1.12903, -1.16492, -1.13773, -1.06396, -1.11805, -1.09131, -1.16439, -1.19177, -1.08141, -1.17604, -1.09649, -1.12952, -1.11991, -1.11633, -1.14776, -1.12007, -1.0591, -1.13377, -1.10087, -1.08388, -1.13641, -1.15817, -1.085, -1.07592, -1.08342, -1.04719, -1.0848, -1.12494, -1.15889, -1.08761, -1.13408, -1.06406, -1.08297, -1.1232, -1.09591, -1.08338, -1.05569, -1.03208, -1.13144, -1.13133, -1.10966, -1.07806, -1.1481, -1.04755, -1.13187, -1.07395, -1.06856, -1.18011, -1.15939, -1.13139, -1.05896, -1.17767, -1.08951, -1.03073, -1.05692, -1.122, -1.17254, -1.14556, -1.08375, -1.11523, -1.08684, -1.13929, -1.07307, -1.11072, -1.06146, -1.09339, -1.10263, -1.04003, -1.15927, -1.07781, -1.11834, -1.10949, -0.99579, -1.07997, -1.0763, -1.16783, -1.11118, -1.18549, -1.09557, -1.09824, -1.09558, -1.11735, -1.113, -1.17004, -1.09844, -1.14848, -1.10039, -1.11218, -1.01426, -1.10208, -1.03295]

    for i in xrange(len(qq)):
        q.set(0,i,qq[i])
    iM = MATRIX(nnucl,1);     iM.set(0,0, 1.0/2000.0)

    # Initial calculations
    q.show_matrix()
    print "=== compute dia =="
    #ham.compute_diabatic(compute_model, q, params, 0)
    ham.compute_diabatic(compute_model, q, params, 1)

    print "=== compute adi =="
    ham.compute_adiabatic(1, 1);
    ham.add_ethd_adi(q, iM, 1)

    for tr in xrange(ntraj):
        print Cpp2Py(ham1[tr].get_full_id())
        print "id = ", ham1[tr].id

    for i in xrange(ntraj):
        print i, "diabatic:", ham.get_ham_dia(Py2Cpp_int([0,i])).show_matrix()
        print i, "adiabatic:", ham.get_ham_adi(Py2Cpp_int([0,i])).show_matrix()

    Ekin, Epot, Etot = compute_etot(ham, p, iM)

    os.system("mkdir energy")
    os.system("mkdir phase_space")
    os.system("mkdir pos_space")
    os.system("mkdir tunnel")

    e = open("energy/energy.txt", "w")
    r = open("phase_space/phase_space.txt", "w")
    g = open("pos_space/pos_space.txt", "w")
    hh = open("tunnel/tunnel.txt", "w")

    # Do the propagation
    for i in xrange(4300):

        Verlet1(dt, q, p, iM, ham, compute_model, params, 1)

        #=========== Properties ==========

        Ekin, Epot, Etot = compute_etot(ham, p, iM)

        e.write( " %8.5f  %8.5f  %8.5f  %8.5f\n" % ( i*dt, Ekin, Epot, Etot ) )

        for j in range(ntraj):
            r.write(" %8.5f  %8.5f" % ( q.get(0,j), p.get(0,j) ) ),
        r.write( "\n" )

        g.write( " %8.5f" % (i*dt) ),
        for j in range(ntraj):
            g.write( " %8.5f" % ( q.get(0,j) ) )
        g.write( "\n" )

        count = 0.0
        for j in range(ntraj):
            if q.get(0,j) > barrier:
                count = count + 1.0
        count = count/ntraj
        hh.write( " %8.5f  %8.5f " % ( i*dt, count ) )
        hh.write( "\n" )

model = 2

#run_test(model, "_0_new.txt")
run_test1(model, "_0_new.txt")
