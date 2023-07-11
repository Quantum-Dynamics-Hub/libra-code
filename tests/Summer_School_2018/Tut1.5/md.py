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

from common import compute_etot, print_xyz
from Hamiltonian import compute_model

kb = (1.9872065E-3 / 627.5094709)  # in Ha now


def run_nve(_nnucl, _ntraj, _q, _p, iM, model, params):
    """
    model - setup the Hamiltonian
    outname - the name of the output file
    """

    ndia, nadi, nnucl, ntraj = 1, 1, _nnucl, _ntraj
    q = MATRIX(_q)
    p = MATRIX(_p)

    # ======= Hierarchy of Hamiltonians =======
    ham = nHamiltonian(ndia, nadi, nnucl)
    ham.init_all(2)

    ham1 = [] 
    for tr in xrange(ntraj):
        ham1.append( nHamiltonian(ndia, nadi, nnucl) )        
        ham1[tr].init_all(2)
        ham.add_child(ham1[tr])


    #  Set up the models and compute internal variables
    # Initial calculations
    ham.compute_diabatic(compute_model, q, params, 1)
    ham.compute_adiabatic(1, 1); 

    Ekin, Epot, Etot = compute_etot(ham, p, iM)

    out1 = open(params["out_energy"], "w"); out1.close()   
    out2 = open(params["out_phase_space"], "w"); out2.close()  
    out3 = open(params["out_positions"], "w"); out3.close()  

    dt = params["dt"]

    Q, P = [], []
    Q.append( MATRIX(q) )
    P.append( MATRIX(p) )


    f = open(params["trajectory_filename"], "w")    
    f.close()


    for i in xrange(params["nsteps"]):

        print_xyz(params["label"], q, params["trajectory_index"], params["trajectory_filename"], i)
        Verlet1(dt, q, p, iM, ham, compute_model, params, 0)

#        sys.exit(0)

        """
        if params["is_periodic"] == 1:
            for dof in xrange(nnucl):
                print dof
                for tr in xrange(ntraj):     
                    if(q.get(dof, tr) > TV):
                        pass #q.add(dof, tr, -1.0*TV)
                    elif(q.get(dof, tr) < 0.0):
                        pass #q.add(dof, tr, TV)
        sys.exit(0)
        """            

        Q.append( MATRIX(q) )
        P.append( MATRIX(p) )

     
        #=========== Properties ==========

        Ekin, Epot, Etot = compute_etot(ham, p, iM)

        Tcurr = 2.0 * Ekin / (nnucl * kb)

        # Print the ensemble average - kinetic, potential, and total energies
        # Print the tunneling information. Here, we count each trajectory across the barrier.
        out1 = open(params["out_energy"], "a")
        out1.write( " %8.5f  %8.5f  %8.5f  %8.5f  %8.5f\n" % ( i*dt, Ekin, Epot, Etot, Tcurr ) )
        out1.close()

        # Print the phase space information
        out2 = open(params["out_phase_space"], "a") 
        for j in range(ntraj):
            out2.write(" %8.5f  %8.5f" % ( q.get(0,j), p.get(0,j) ) ),
        out2.write( "\n" )
        out2.close()  

        # Print the position versus time infromation
        out3 = open(params["out_positions"], "a")
        out3.write( " %8.5f" % (i*dt) ),
        for j in range(ntraj):
            out3.write( " %8.5f" % ( q.get(0,j) ) )
        out3.write( "\n" )
        out3.close()  


    return Q, P




def run_nvt(_nnucl, _ntraj, _q, _p, iM, model, params):
    """
    model - setup the Hamiltonian
    outname - the name of the output file
    """

    Q, P = [], []

    ndia, nadi, nnucl, ntraj = 1, 1, _nnucl, _ntraj
    q = MATRIX(_q)
    p = MATRIX(_p)

    # ======= Hierarchy of Hamiltonians =======
    ham = nHamiltonian(ndia, nadi, nnucl)
    ham.init_all(2)

    ham1 = [] 
    for tr in xrange(ntraj):
        ham1.append( nHamiltonian(ndia, nadi, nnucl) )        
        ham1[tr].init_all(2)
        ham.add_child(ham1[tr])


    #  Set up the models and compute internal variables


    # Initial calculations
    ham.compute_diabatic(compute_model, q, params, 1)
    ham.compute_adiabatic(1, 1); 

    therms = Thermostat(params)
    therms.set_Nf_t(nnucl*ntraj)
    therms.set_Nf_r(0)
    therms.set_Nf_b(0)
    therms.init_nhc()
    Ebath = therms.energy()

    therms.show_info()

#    sys.exit()


    Ekin, Epot, Etot = compute_etot(ham, p, iM)


    out1 = open(params["out_energy"], "w"); out1.close()   
    out2 = open(params["out_phase_space"], "w"); out2.close()  
    out3 = open(params["out_positions"], "w"); out3.close()  

    dt = params["dt"]

    Q.append( MATRIX(q) )
    P.append( MATRIX(p) )



    f = open(params["trajectory_filename"], "w")    
    f.close()

    # Do the propagation
    for i in xrange(params["nsteps"]):

        print_xyz(params["label"], q, params["trajectory_index"], params["trajectory_filename"], i)
        
        Verlet1_nvt(dt, q, p, iM, ham, compute_model, params, 0, therms)

        Q.append( MATRIX(q) )
        P.append( MATRIX(p) )
     
        #=========== Properties ==========

        Ekin, Epot, Etot = compute_etot(ham, p, iM)
        Ebath = therms.energy() / float(ntraj)

        Tcurr = 2.0 * Ekin / (nnucl * kb)

        # Print the ensemble average - kinetic, potential, and total energies
        # Print the tunneling information. Here, we count each trajectory across the barrier.
        out1 = open(params["out_energy"], "a")
        out1.write( " %8.5f  %8.5f  %8.5f  %8.5f   %8.5f  %8.5f  %8.5f\n" % ( i*dt, Ekin, Epot, Etot, Ebath, Etot+Ebath, Tcurr ) )
        out1.close()

        # Print the phase space information
        out2 = open(params["out_phase_space"], "a") 
        for j in range(ntraj):
            out2.write(" %8.5f  %8.5f" % ( q.get(0,j), p.get(0,j) ) ),
        out2.write( "\n" )
        out2.close()  

        # Print the position versus time infromation
        out3 = open(params["out_positions"], "a")
        out3.write( " %8.5f" % (i*dt) ),
        for j in range(ntraj):
            out3.write( " %8.5f" % ( q.get(0,j) ) )
        out3.write( "\n" )
        out3.close()  

    return Q, P


