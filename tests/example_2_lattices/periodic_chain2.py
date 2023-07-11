#*********************************************************************************
#* Copyright (C) 2017 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/

import os
import sys
import math
import copy

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from libra_py import *

"""
  This example demonstrates the calculation of the frequnecy analysis
  for 1D N-particle chains. 

  We run an MD simulation to find and collect a "dipole moment" ACF. It's fourier transform
  gives the frequencies. Note that here we quote the dipole moment because the charges assigned
  to the particles are somewhat random - but depending in their sizes, one may or may not be
  able to "see" particular frequencies in the spectrum. Play with the charges to find out.
  Generally, less simmetries there are, more frequencies are observed

  We then run the Hessian calculations to find out the true Harmonic normal modes. 
  For the 2-particle system, we can compare it the the analytic result.

  You can turn the periodic interaction on/off (presently they are off).

  All masses are assumed to be the same.

  In the MD, momenta are initialized randomly, but only in 1D (along x axis), so we don't get other
  modes in the spectrum. You can initialize the velocities along y and z directions to see what happens.





"""

inv_cm2ev = (1.0/8065.54468111324)
ev2Ha = (1.0/27.211)    # 27.2 ev is 1 Ha 
inv_cm2Ha = inv_cm2ev * ev2Ha


def init_system(N, dist_init, mass, scl):

    rnd = Random()

    # Associate the coordinates with a molecular system
    U = Universe() 
    LoadPT.Load_PT(U, "elements.dat", 0)
    syst = System()  # chemical system

    dist = 0.0
    for i in xrange(N):
        R = VECTOR(dist, 0.0, 0.0)
        dist = dist + dist_init
        syst.CREATE_ATOM( Atom(U,  {"Atom_element":"H","Atom_cm_x":R.x,"Atom_cm_y":R.y,"Atom_cm_z":R.z})  )

        # Masses
        syst.Atoms[i].Atom_RB.set_mass(mass)

        # Momenta
        p = VECTOR(scl.x*rnd.normal(), scl.y*rnd.normal(), scl.z*rnd.normal())
        syst.Atoms[i].Atom_RB.set_momentum(p);


    return syst



def init_interactions(interactions, syst, Hess):
    """
    interactions - (list of Interaction_N_Body or derived objects). Represents interactions
    syst - (System) - Chemical system object
    Hess - (MATRIX) - A Hessian matrix

    This function sets energy, atomic stress, all forces and Hessian matrix to zero
    """

    for inter in interactions:
        inter.energy = 0.0
        inter.stress_at *= 0.0

    syst.zero_atom_forces()

    Hess *= 0.0



def set_interactions_1D(syst, Hess, t0, t, K, dist, is_periodic):
    """
    syst - (System) - Chemical system object
    Hess - (MATRIX) - Hessian 
    t0 - (VECTOR) - this is the zero vector, has to be defined outside, since the interactions will
    be using a pointer to that object, so the object has to exist for the entire duration of the
    simulation. This is also true for other veriables used here: R, F, Hess, t
    t0 - (VECTOR) - the 1D unit cell periodic translation vector
    K - (double) - force constant for the Harmonic potential of the type U = k*(r-r_0)^2
    dist - (double) - r_0 in the above equation
    is_periodic (int) - a flag to tell whether the first and the last particles in the chain interact
    via a spring

    Define interactions in the system. 
    At this point, we construct them as a chain of 1D Harmonic oscillators, which may or may not
    by periodically interacting.    
    Once the variables are set, the all the calculations are also run.

    """

    N = syst.Number_of_atoms

    inter = []
        
    for i in xrange(N-1):
        inter.append( Bond_Interaction() )
        inter[i].set_functional("Harmonic")
        inter[i].set_params({"K":K, "r0":dist })
        inter[i].set_coords(syst.Atoms[i].Atom_RB.rb_cm, syst.Atoms[i+1].Atom_RB.rb_cm)
        inter[i].set_transl(t0,t0)
        inter[i].set_forces(syst.Atoms[i].Atom_RB.rb_force, syst.Atoms[i+1].Atom_RB.rb_force)
        inter[i].set_hessian(Hess, Py2Cpp_int([3*i, 3*i+1, 3*i+2, 3*(i+1), 3*(i+1)+1, 3*(i+1)+2]))

        inter[i].compute()
        #print i, inter[i].energy


    if is_periodic:
        # Periodic interaction
        inter.append( Bond_Interaction() )
        inter[N-1].set_functional("Harmonic")
        inter[N-1].set_params({"K":K, "r0":dist })
        inter[N-1].set_coords(syst.Atoms[0].Atom_RB.rb_cm, syst.Atoms[N-1].Atom_RB.rb_cm)
        inter[N-1].set_transl(t, t0)
        inter[N-1].set_forces(syst.Atoms[0].Atom_RB.rb_force, syst.Atoms[N-1].Atom_RB.rb_force)
        inter[N-1].set_hessian(Hess, Py2Cpp_int([0, 1, 2, 3*(N-1), 3*(N-1)+1, 3*(N-1)+2]))

        inter[N-1].compute()
        #print N-1, inter[N-1].energy


    return inter



def potential(interactions, syst, Hess):
    """
    interactions - (list of Interaction_N_Body or derived objects). Represents interactions
    syst - (System) - Chemical system object
    Hess - (MATRIX) - A Hessian matrix

    Compute potential energy and all related variables. This function also
    takes care of zeroing all the forces, stress, and Hessian before all the
    interactions are computed.
    """

    init_interactions(interactions, syst, Hess)

    energy = 0.0

    N = len(interactions)
    for i in xrange(N):
        interactions[i].compute()
        #print i, interactions[i].energy

        energy = energy + interactions[i].energy
    return energy    



def kinetic(syst):
    """
    syst - (System) - Chemical system object

    Compute the kinetic energy    
    """

    res = syst.ekin_tr_atom()

    return res





def md( params ):
    """
    params - a dictionary of the control parameters

    Run an MD simulation of the chain of N particles
    """



    N = params["N"]
    dist_init = params["dist_init"]
    scl = params["scl"]
    mass = params["mass"]

    K = params["K"]
    dist_eq = params["dist_eq"]
    is_periodic = params["is_periodic"]

    dt = params["dt"]
    Nsteps = params["Nsteps"]

    efile = params["energy_filename"]
    tfile = params["trajectory_filename"]
    acffile = params["acf_filename"]
    specfile = params["spectrum_filename"]
   
    # "Charges"
    q = []
    if N==2:
        q = [1.0, -1.0]    
    elif N==4:
        q = [1.0,  -1.0,  3.0, -3.0]
    elif N==6:
        q = [1.0,  -1.0,  3.0, -3.0, 6.0, -6.0]
    else:
        q = [0.0]*N


    # Hessian
    Hess = MATRIX(3*N,3*N)

    # Translation vectors
    t0, t = VECTOR(0.0, 0.0, 0.0),  VECTOR(N*dist_eq, 0.0, 0.0)


    # Create the molecular system
    syst = init_system(N, dist_init, mass, scl)

           
    # Initialize the interactions
    inter = set_interactions_1D(syst, Hess, t0, t, K, dist_eq, is_periodic)

    # Compute all the energies
    epot = potential(inter, syst, Hess)
    ekin = kinetic(syst)
    etot = ekin + epot
    tim = 0.0


    f = open(efile, "w")
    f.write("%8.5f %8.5f %8.5f %8.5f \n" % (tim, ekin, epot, etot))  


    ##============ MD ==============================
    mu = []

    for j in xrange(Nsteps):

        #====== MD: Verlet propagation ==========
        for i in xrange(N):
            syst.Atoms[i].Atom_RB.apply_force(0.5*dt)

            iM = syst.Atoms[i].Atom_RB.rb_iM
            Pi = syst.Atoms[i].Atom_RB.rb_p
            syst.Atoms[i].Atom_RB.shift_position(dt*iM*Pi)


        epot = potential(inter, syst, Hess)

        for i in xrange(N):
            syst.Atoms[i].Atom_RB.apply_force(0.5*dt)

        ekin = kinetic(syst)
        etot = ekin + epot
        tim = tim + dt

        f.write("%8.5f %8.5f %8.5f %8.5f \n" % (tim, ekin, epot, etot))

        #====== Compute dipole moment and collect it =======
        
        ave_mu = VECTOR(0.0, 0.0, 0.0)

        for i in xrange(N):
            ave_mu = ave_mu + q[i]*syst.Atoms[i].Atom_RB.rb_cm
        ave_mu = ave_mu / float(N)
        mu.append(ave_mu)


        #====== Print out the trajectories ===============
        syst.print_xyz(tfile, j)


    f.close()


    # Compute the ACF and spectra    
    acf_vector.recipe1(mu, dt*0.02419, 5000.0, 1.0, acffile, specfile, 0)



def normal_modes( params ):
    """
    params - a dictionary of the control parameters

    Compute the normal modes of the system
    """


    N = params["N"]
    dist_init = params["dist_init"]
    scl = params["scl"]
    mass = params["mass"]

    K = params["K"]
    dist_eq = params["dist_eq"]
    is_periodic = params["is_periodic"]

    Hess_file = params["Hess_filename"]
    nm_file = params["normal_modes_filename"]

    # Coordinates, momenta, forces, masses
    syst = init_system(N, dist_init, mass, scl)

    # Mass matrix
    G = MATRIX(3*N, 3*N)
    for i in xrange(3*N):
        G.set(i,i, (1.0/math.sqrt(mass)) )


    # Hessian
    Hess = MATRIX(3*N,3*N)

    # Translation vectors
    t0, t = VECTOR(0.0, 0.0, 0.0),  VECTOR(N*dist_eq, 0.0, 0.0)



    # Initialize the interactions
    inter = set_interactions_1D(syst, Hess, t0, t, K, dist_eq, is_periodic)

    # Compute all the energies
    epot = potential(inter, syst, Hess)

    Hess.show_matrix(Hess_file)

    dynMat = CMATRIX(Hess)    

    C = CMATRIX(3*N, 3*N)
    omega = CMATRIX(3*N, 3*N)
    I = CMATRIX(3*N, 3*N); I.identity();


    solve_eigen(dynMat, I, omega, C, 0)

    
    line = ""
    for i in xrange(3*N):
        line = line + "Mode %5i " % (i)
        if omega.get(i,i).real > 0.0:
            line = line + "Real %8.5f cm^-1 : " % (math.sqrt(omega.get(i,i).real/mass) /inv_cm2Ha) 
        else:
            line = line + "Imag %8.5f cm^-1 : " % (math.sqrt(abs (omega.get(i,i).real)/mass) /inv_cm2Ha) 

        for j in xrange(3*N):
            line = line + " %8.5f " % (C.get(j,i).real)

        line = line + "\n"

    f = open(nm_file, "w")        
    f.write(line)
    f.close()



###===================== Test calculations ============================   

## 2-particle system
params = {"N":2, "dist_init":2.1, "scl":VECTOR(0.0, 0.0, 0.0), "mass":2000.0,
         "K":0.1, "dist_eq":2.0, "is_periodic":0,         
         "dt":20.0, "Nsteps":500, 
         "trajectory_filename":"test1/chain.xyz", "energy_filename":"test1/energy.txt", "acf_filename":"test1/acf.txt", "spectrum_filename":"test1/spectrum.txt"
         }

md(params)

params["dist_init"] = 2.0
params["Hess_filename"] = "test1/Hessian.txt"
params["normal_modes_filename"] = "test1/Modes.txt"

normal_modes(params)


# We use 2K because the potential is of the form K*(r-r0)^2, and not 1/2*K*(r-r0)^2
print "Frequency = ", 2.0*math.sqrt(0.1/2000.0), " rad/a.u. (a.u.)"
print "Frequency = ", 2.0*math.sqrt(0.1/2000.0) / inv_cm2Ha, " cm^-1"




