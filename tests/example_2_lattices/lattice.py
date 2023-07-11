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
  This file contains the functions and procedures for frequnecy analysis and md of 1, 2, and 3D lattices
  or clusters

  We run an MD simulation to find and collect a "dipole moment" ACF. It's fourier transform
  gives the frequencies. Note that here we quote the dipole moment because the charges assigned
  to the particles are somewhat random - but depending in their sizes, one may or may not be
  able to "see" particular frequencies in the spectrum. Play with the charges to find out.
  Generally, less simmetries there are, more frequencies are observed

  We then run the Hessian calculations to find out the true Harmonic normal modes. 

  You can turn the periodic interaction on/off (presently they are off).

  In the MD, momenta are initialized randomly


"""

inv_cm2ev = (1.0/8065.54468111324)
ev2Ha = (1.0/27.211)    # 27.2 ev is 1 Ha 
inv_cm2Ha = inv_cm2ev * ev2Ha



def compute_S_H_el(bas, H_dia, H_adi, orb_dia, orb_adi, Nocc):
    ##
    # Compute S and H matrices for the chosen space of active fragment MOs (diabatic MOs)
    #
    # \param[in] bas - (list of AO objects of length Nao) AO basis 
    # \param[in] H_dia - diabatic Hamiltonian in the AO basis (H0) - MATRIX(Nao,Nao) 
    # \param[in] H_adi - adiabatic Hamiltonian in the AO basis (H1 = H0 + H') - MATRIX(Nao,Nao) 
    # \param[in] orb_dia - Indices of the MOs from the ham_dia set included in the active space
    # \param[in] orb_adi - Indices of the MOs from the ham_adi set included in the active space

    # Dimensions
    Nao = len(bas)           # the number of AO basis functions
    N_dia = len(orb_dia)     # the total number of diabatic MOs
    N_adi = len(orb_adi)     # the total number of adiabatic MOs


    # Compute the AO overlap
    um = MATRIX(Nao, Nao);  um.Init_Unit_Matrix(1.0)
    Sao = MATRIX(Nao, Nao); 

    dist2 = 10000.0
    MO_overlap(Sao, bas, bas, um, um, Py2Cpp_int(range(0,Nao)), Py2Cpp_int(range(0,Nao)), dist2) 

    # Compute the MOs
    # Complex-valued version
    res_dia = Fock_to_P(CMATRIX(H_dia), CMATRIX(Sao), Nocc, 1.0, 0.001, 1e-8, 0)    
    res_adi = Fock_to_P(CMATRIX(H_adi), CMATRIX(Sao), Nocc, 1.0, 0.001, 1e-8, 0)    

    E_dia = res_dia[0] # MO energies
    E_adi = res_adi[0] # MO energies

    C_dia = res_dia[1] # MO-LCAO
    C_adi = res_adi[1] # MO-LCAO


    # Result matrices
    # PSI = sum_i { C_i^dia * |dia_i> } = sum_j { C_j^adi * |adi_j> }
    # <adi_k|PSI> = sum_i { C_i^dia * <adi_k|dia_i> } = C_k^adi, so C^adi = L * C_dia
    # <dia_k|PIS> = C_k^dia = sum_j { C_j^adi * <dia_k|adi_j> } = L.H() * C^adi = C^dia
    L = CMATRIX(N_adi, N_dia)  # L = <adi|dia>  - projector

    # P = <adi|dia>
    MO_overlap(L, ham_old.basis_ao, bas, C_adi, C_dia, Py2Cpp_int(orb_adi), Py2Cpp_int(orb_dia), dist2) 

    
    # I = <dia|adi><adi|dia>
    I = L.H() * L

    # H = L^+ * E_adi * L  - this is the approximation of the H_adi in the diabatic basis
    e_adi = CMATRIX(N_adi, N_adi)
    pop_submatrix(E_adi, e_adi, orb_adi, orb_adi)
    H_el = L.H() * e_adi * L

    # E_dia
    e_dia = CMATRIX(N_dia, N_dia)
    pop_submatrix(E_dia, e_dia, orb_dia, orb_dia)


    # Perturbation in the dia basis
    V = C_dia.H() * ( H_adi - H_dia ) * C_dia   

    V_dia = CMATRIX(N_dia, N_dia)
    pop_submatrix(V, V_dia, orb_dia, orb_dia)


    t.stop();  print "Overlap and electronic Hamiltonian = ", t.show()

    return L, H_el, I, V_dia, e_dia, e_adi, C_dia, C_adi



def compute_Hvib(ham_old, ham_cur, orb_dia, orb_adi, dt, orb_prime):
    ##         
    # Computes the vibronic Hamiltonian
    # \param[in] ham_old (Hamiltonian object) - unperturbed at t-dt
    # \param[in] ham_cur (Hamiltonian object) - unperturbed at t
    # \param[in] orb_dia (list of ints) - indices of the orbitals included in the diabatic active space
    # \param[in] orb_adi (list of ints) - indices of the orbitals included in the adiabatic active space
    # \param[in] dt (float) Time step (in a.u.)
    # \param[in] orb_prime (list of ints) - indiced of the orbitals coupled by the perturbation
    # \param[in] H_prime (float) - the magnitude of the perturbation

    N_dia = len(orb_dia)
    N_adi = len(orb_adi)

    es_old = ham_old.get_electronic_structure()
    es_cur = ham_cur.get_electronic_structure()

    Nocc = es_old.Nocc_alp
    H_dia = CMATRIX(es_old.get_Fao_alp() )
    H_adi = perturb(H_dia, orb_prime, H_prime)
    L_old, H_el_old, I_old, V_old, e_dia_old, e_adi_old, C_dia_old, C_adi_old  = compute_S_H_el(ham_old.basis_ao, H_dia, H_adi, orb_dia, orb_adi, Nocc)


    Nocc = es_cur.Nocc_alp    
    H_dia = CMATRIX( es_cur.get_Fao_alp() )
    H_adi = perturb(H_dia, orb_prime, H_prime)
    L_cur, H_el_cur, I_cur, V_cur, e_dia_cur, e_adi_cur, C_dia_cur, C_adi_cur = compute_S_H_el(ham_cur.basis_ao, H_dia, H_adi, orb_dia, orb_adi, Nocc)


    # <MO(t)|MO(t-dt)>

    D_dia = CMATRIX(N_dia,N_dia)
    MO_overlap(L, ham_old.basis_ao, ham_cur.basis_ao, C_old, C_cur,
               Py2Cpp_int(orb_dia), Py2Cpp_int(orb_dia), 1e+10
              )  # 
    D_dia = (0.5/dt)*(D_dia - D_dia.H())


    D_adi = CMATRIX(N_adi,N_adi)
    MO_overlap(D_adi, 
               ham_old.basis_ao, ham_cur.basis_ao, 
               C_adi_old, C_adi_cur,
               Py2Cpp_int(orb_adi), Py2Cpp_int(orb_adi), dist2
              )  # 
    D_adi = (0.5/dt)*(D_adi - D_adi.H())


    Hel = 0.5*(H_el_old + H_el_cur) 
    V = 0.5*(V_old + V_cur)
    Edia = 0.5*(e_dia_old + e_dia_cur)
    Eadi = 0.5*(e_adi_old + e_adi_cur)
        
    Hvib1 = Hel - 1.0j*D_dia        # spectral
    Hvib2 = Edia + V - 1.0j*D_dia   # perturbative
    Hvib3 = Eadi - 1.0j*D_adi       # direct adiabatic


    datautils.show_matrix_splot(Hvib1.imag(), "_Hvib1_im.dat")
    datautils.show_matrix_splot(Hvib2.imag(), "_Hvib2_im.dat")
    datautils.show_matrix_splot(Hvib3.imag(), "_Hvib3_im.dat")

    datautils.show_matrix_splot(Hvib1.real(), "_Hvib1_re.dat")
    datautils.show_matrix_splot(Hvib2.real(), "_Hvib2_re.dat")
    datautils.show_matrix_splot(Hvib3.real(), "_Hvib3_re.dat")

    

    return Hvib1, Hvib2, Hvib3, L_cur




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




def set_interactions_2B(syst, Hess, pairs, K, dist):
    """
    This is a general setter of the 2-body interactions

    syst - (System) - Chemical system object
    Hess - (MATRIX) - Hessian 
    pairs - ([int, int, VECTOR, VECTOR]) - defines the interacting pair, the two integers are the
         indices of the atoms, the two vectors are their periodic translations
    K - (double) - force constant for the Harmonic potential of the type U = k*(r-r_0)^2
    dist - (double) - r_0 in the above equation

    Define interactions in the system. 
    Once the variables are set, the all the calculations are also run.

    """

    N_int = len(pairs)

    inter = []
        
    for i in xrange(N_int):
        I, J = pairs[i][0], pairs[i][1]
        inter.append( Bond_Interaction() )
        inter[i].set_functional("Harmonic")
        inter[i].set_params({"K":K, "r0":dist })
        inter[i].set_coords(syst.Atoms[I].Atom_RB.rb_cm, syst.Atoms[J].Atom_RB.rb_cm)
        inter[i].set_transl(pairs[i][2], pairs[i][3])
        inter[i].set_forces(syst.Atoms[I].Atom_RB.rb_force, syst.Atoms[J].Atom_RB.rb_force)
        inter[i].set_hessian(Hess, Py2Cpp_int([3*I, 3*I+1, 3*I+2, 3*J, 3*J+1, 3*J+2]))

#    for i in xrange(N_int):
#        inter[i].compute()
#        print i, inter[i].energy


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





def md(syst, R, MaxCoord, params):
    """
    params - a dictionary of the control parameters

    Run an MD simulation of the chain of N particles
    """


    Nx, Ny, Nz = params["Nx"], params["Ny"], params["Nz"]
    a, b, c = params["a"], params["b"], params["c"]
    Rcut = params["Rcut"]
    pbc_opt = params["pbc_opt"]

    K = params["K"]
    dist_eq = params["dist_eq"]
    is_periodic = params["is_periodic"]

    dt = params["dt"]
    Nsteps = params["Nsteps"]
    is_elstr = params["is_elstr"]  # whether to do electronic structure calculations or not
    is_mu_acf = params["is_mu_acf"]

    efile = params["energy_filename"]
    tfile = params["trajectory_filename"]
    acffile = params["acf_filename"]
    specfile = params["spectrum_filename"]


    # Initialize the QM Hamiltonian at two timesteps
    prms = Control_Parameters(); 
    control_filename = "control_parameters_eht.dat"

    
    res, line, pairs = autoconnect.autoconnect_pbc(R, MaxCoord, Rcut, Nx*a, Ny*b, Nz*c, pbc_opt, 0, 0)
    Nat = len(R)  # number of atoms
   
    # "Charges"
    q = [0.0]*Nat


    # Hessian
    Hess = MATRIX(3*Nat,3*Nat)

    # Initialize the interactions
    inter = set_interactions_2B(syst, Hess, pairs, params["K"], dist_eq)

    # Compute all the energies
    epot = potential(inter, syst, Hess)
    ekin = kinetic(syst)
    etot = ekin + epot
    tim = 0.0


    #======= Electronic structure =======================
    ham_cur, ham_old = None, None

    if is_elstr==1:
        ham_cur = listHamiltonian_QM(control_filename, syst )  # this includes the overlap calculation
        ham_cur.compute_scf(syst)                              # this includes the core Hamiltonian calculations

        es = ham_cur.get_electronic_structure()
        homo = es.Nocc_alp - 1  # index of the HOMO orbital in the adiabatic sub-system
        print "Adiabatic sub-system: number of electrons (alpha) = ", es.Nocc_alp, " homo = ", homo       



    f = open(efile, "w")
    f.write("%8.5f %8.5f %8.5f %8.5f \n" % (tim, ekin, epot, etot))  

    ##============ MD ==============================
    mu = []

    for j in xrange(Nsteps):

        #====== MD: Verlet propagation ==========
        for i in xrange(Nat):
            syst.Atoms[i].Atom_RB.apply_force(0.5*dt)

            iM = syst.Atoms[i].Atom_RB.rb_iM
            Pi = syst.Atoms[i].Atom_RB.rb_p
            syst.Atoms[i].Atom_RB.shift_position(dt*iM*Pi)

        epot = potential(inter, syst, Hess)

        for i in xrange(Nat):
            syst.Atoms[i].Atom_RB.apply_force(0.5*dt)


        ekin = kinetic(syst)
        etot = ekin + epot
        tim = tim + dt

        f.write("%8.5f %8.5f %8.5f %8.5f \n" % (tim, ekin, epot, etot))


        #=============== Now handle the electronic structure ==========================
        if is_elstr==1:
            ham_old = listHamiltonian_QM(ham_cur)        
        
            # Update positions of the basis AOs!!!
            for o in xrange(ham_cur.Norb):
                at_indx = ham_cur.ao_to_atom_map[o]  # indices of the atoms in the system bound to this Hamiltonin object, that
                                                     # is syst_adi !!!  so we will need another mapping of the indices of those atoms to
                                                     # the indices to the overall system atoms
                ham_cur.basis_ao[o].set_position(syst.Atoms[at_indx].Atom_RB.rb_cm)

            # Update overlaps and core Hamiltonian
            ham_cur.compute_overlap(syst)
            ham_cur.compute_scf(syst)       # this includes the core Hamiltonian calculations
            es = ham_cur.get_electronic_structure() 

            H_dia = CMATRIX( es.get_Fao_alp() )
            Hvib, L = compute_Hvib(ham_old, ham_cur, orb_dia, orb_adi, dt)


        #====== Compute dipole moment and collect it =======
        if is_mu_acf==1:        
            ave_mu = VECTOR(0.0, 0.0, 0.0)

            for i in xrange(Nat):
                ave_mu = ave_mu + q[i]*syst.Atoms[i].Atom_RB.rb_cm
            ave_mu = ave_mu / float(Nat)
            mu.append(ave_mu)


        #====== Print out the trajectories ===============
        syst.print_xyz(tfile, j)


    f.close()


    # Compute the ACF and spectra    
    if is_mu_acf==1:
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

rnd = Random()

# Unit cell parameters
a = VECTOR(2.0, 0.0, 0.0)
b = VECTOR(0.0, 2.0, 0.0)
c = VECTOR(0.0, 0.0, 2.0)

# Create a Chemical system
Nx, Ny, Nz = 5, 5, 5  # 

syst, R, MaxCoord = System(), [], []
build.add_atom_to_system(syst, R, MaxCoord, Nx,Ny,Nz, a, b, c, VECTOR(0.0, 0.0, 0.0), "H", 2000.0, VECTOR(1.0, 1.0, 1.0), 6, rnd)


# Setup parameters of MD

params = {"Nx":Nx, "Ny":Ny, "Nz":Nz,
          "a":a, "b":b, "c":c, "Rcut":2.001, "pbc_opt":"abc",
          "K":0.1, "dist_eq":2.0, "is_periodic":0,         
          "dt":20.0, "Nsteps":500, "is_elstr":0, "is_mu_acf":0,
          "trajectory_filename":"test1/chain.xyz", "energy_filename":"test1/energy.txt", "acf_filename":"test1/acf.txt", "spectrum_filename":"test1/spectrum.txt"
         }

md(syst, R, MaxCoord, params)



sys.exit()

params["dist_init"] = 2.0
params["Hess_filename"] = "test1/Hessian.txt"
params["normal_modes_filename"] = "test1/Modes.txt"

#normal_modes(params)


# We use 2K because the potential is of the form K*(r-r0)^2, and not 1/2*K*(r-r0)^2
#print "Frequency = ", 2.0*math.sqrt(0.1/2000.0), " rad/a.u. (a.u.)"
#print "Frequency = ", 2.0*math.sqrt(0.1/2000.0) / inv_cm2Ha, " cm^-1"

