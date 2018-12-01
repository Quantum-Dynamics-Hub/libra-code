#*********************************************************************************
#* Copyright (C) 2017 Brendan A. Smith, Wei Li, Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
#
#  
#
import sys
import cmath
import math
import os

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
#from libra_py import *

import mapping
import libra_py.workflows.common_utils as comn
import libra_py.tsh as tsh
import libra_py.units as units
import libra_py.hungarian as hungarian



def show_matrix_splot(X, filename):
    ncol, nrow = X.num_of_cols, X.num_of_rows

    line = ""
    for i in xrange(nrow):
        for j in xrange(ncol):
            val = X.get(i,j)
            if i==j:
                val = 0.0
            line = line + "%4i %4i %8.5f \n" % (i, j, val)
        line = line + "\n"

    f = open(filename, "w")
    f.write(line)
    f.close()
 
    

def add_printout(i, pop, filename):
    # pop - CMATRIX(nstates, 1)

    f = open(filename,"a")
    line = "step= %4i " % i    

    tot_pop = 0.0
    for st in xrange(pop.num_of_cols):
        pop_o = pop.get(st,st).real
        tot_pop = tot_pop + pop_o
        line = line + " P(%4i)= %8.5f " % (st, pop_o)
    line = line + " Total= %8.5f \n" % (tot_pop)
    f.write(line)
    f.close()



def get_data(params):
    """
    Read the "elementary" overlaps and energies in the basis of KS orbitals

    Required parameter keys:

    params["norbitals"]        [int] - how many lines/columns in the file - the total number of spin-orbitals
    params["active_space"]     [list of ints] - which orbitals we care about (indexing starts with 0)
    params["S_re_prefix"]      [string] - prefixes of the files with real part of the MO overlaps at time t
    params["S_re_suffix"]      [string] - suffixes of the files with real part of the MO overlaps at time t
    params["S_im_prefix"]      [string] - prefixes of the files with imaginary part of the MO overlaps at time t
    params["S_im_suffix"]      [string] - suffixes of the files with imaginary part of the MO overlaps at time t
    params["St_re_prefix"]     [string] - prefixes of the files with real part of the MO overlaps at times t and t' = t+dt
    params["St_re_suffix"]     [string] - suffixes of the files with real part of the MO overlaps at times t and t' = t+dt
    params["St_im_prefix"]     [string] - prefixes of the files with imaginary part of the MO overlaps at times t and t' = t+dt
    params["St_im_suffix"]     [string] - suffixes of the files with imaginary part of the MO overlaps at times t and t' = t+dt
    params["E_re_prefix"]      [string] - prefixes of the files with real part of the MO energies at time t
    params["E_re_suffix"]      [string] - suffixes of the files with real part of the MO energies  at time t
    params["E_im_prefix"]      [string] - prefixes of the files with imaginary part of the MO energies at time t
    params["E_im_suffix"]      [string] - suffixes of the files with imaginary part of the MO energies  at time t
    params["nsteps"]           [int] - how many files to read, starting from index 0
    params["is_pyxaid_format"] [Boolean] - whether the input in the old PYXAID format (just orbital space matrices)

    """

    critical_params = ["norbitals", "active_space", "S_re_prefix", "S_im_prefix",
    "St_re_prefix", "St_im_prefix", "E_re_prefix", "E_im_prefix", "nsteps" ]

    default_params = { "is_pyxaid_format":False, "S_re_suffix":"_re", "S_suffix":"_im",
    "St_re_suffix":"_re", "St_im_suffix":"_im", "E_re_suffix":"_re", "E_im_suffix":"_im"}

    comn.check_input(params, default_params, critical_params)

    norbitals = params["norbitals"]  # the number of orbitals in the input files
    active_space = params["active_space"]
    nstates = len(active_space)
    nsteps = params["nsteps"]

    S, St, E = [], [], [] 

    for i in range(0,nsteps):

        filename_re = params["S_re_prefix"]+str(i)+params["S_re_suffix"]
        filename_im = params["S_im_prefix"]+str(i)+params["S_im_suffix"]
        s = comn.get_matrix(norbitals, norbitals, filename_re, filename_im, active_space )
        if params["is_pyxaid_format"]==True:
            s = comn.orbs2spinorbs(s)
        S.append(s)

        filename_re = params["St_re_prefix"]+str(i)+params["St_re_suffix"]
        filename_im = params["St_im_prefix"]+str(i)+params["St_im_suffix"]
        st = comn.get_matrix(norbitals, norbitals, filename_re, filename_im, active_space )
        if params["is_pyxaid_format"]==True:
            st = comn.orbs2spinorbs(st)
        St.append(st)

        filename_re = params["E_re_prefix"]+str(i)+params["E_re_suffix"]
        filename_im = params["E_im_prefix"]+str(i)+params["E_im_suffix"]
        e = comn.get_matrix(norbitals, norbitals, filename_re, filename_im, active_space )
        if params["is_pyxaid_format"]==True:
            e = comn.orbs2spinorbs(e)
        E.append(e)



    return S, St, E


def apply_state_reordering(St, E, params):
    """
    St [list of CMATRIX]
    E [list of CMATRIX] 
    """

    critical_params = [ ]
    default_params = { "do_state_reordering":2, "state_reordering_alpha":0.0 }
    comn.check_input(params, default_params, critical_params)


    nsteps = len(St)
    nstates = St[0].num_of_cols


    # Initialize the cumulative permutation as the identity permutation
    perm_cum = intList() # cumulative permutation
    for a in xrange(nstates):
        perm_cum.append(a)

    # Current permutation
    perm_t = intList() 
    for a in xrange(nstates):
        perm_t.append(a)


    for i in range(0, nsteps):
    
        ### Perform state reordering (must be done before the phase correction) ###
    
        if params["do_state_reordering"]==1:
            """
            A simple approach based on permuations - but this is not robust
            may have loops
            """
            perm_t = get_reordering(St[i])

            # apply the cumulative permutation  
            update_permutation(perm_t, perm_cum)

            # apply the permutation
            # Because St = <psi(t)|psi(t+dt)> - we permute only columns
            St[i].permute_cols(perm_cum)

            E[i].permute_cols(perm_cum)
            E[i].permute_rows(perm_cum)
    
    
        elif params["do_state_reordering"]==2:
            """
            The Hungarian approach
            """

            # S' = Pn^+ * S
            St[i].permute_rows(perm_t)  # here we use the old value, perm_t = P_n 

            # construct the cost matrix
            alp = 0.0
            if "state_reordering_alpha" in params.keys():
                alp = params["state_reordering_alpha"]
    
            cost_mat = MATRIX(nstates, nstates)
            for a in xrange(nstates):
                for b in xrange(nstates):
                    s = St[i].get(a,b)
                    s2 =  (s*s.conjugate()).real
                    dE = (E[i].get(a,a)-E[i].get(b,b)).real
                    val = s2 * math.exp(-alp*dE**2)
                    cost_mat.set(a,b, val)
    
            # run the assignment calculations
            res = hungarian.maximize(cost_mat)
    
            # convert the list of lists into the permutation object
            for r in res:
                perm_t[r[0]] = r[1]   # now, this becomes a new value: perm_t = P_{n+1}
    
            # apply the permutation
            # Because St = <psi(t)|psi(t+dt)> - we permute only columns
            St[i].permute_cols(perm_t)

            #E[i].permute_cols(perm_t)
            #E[i].permute_rows(perm_t)






def apply_phase_correction(St, params):
    """
    Perform the phase correction according to:
        
    Akimov, A. V. J. Phys. Chem. Lett, 2018, 9, 6096

    """

    critical_params = [ ]
    default_params = { "do_phase_correction":1 }
    comn.check_input(params, default_params, critical_params)


    nsteps = len(St)
    nstates = St[0].num_of_cols

    ### Initiate the cumulative phase correction factors ###    
    cum_phase = CMATRIX(nstates,1)  # F(n-1)  cumulative phase
    for a in xrange(nstates):
        cum_phase.set(a, 0, 1.0+0.0j)


    for i in range(0, nsteps):

        if params["do_phase_correction"]==1:
            ### Compute the instantaneous phase correction factors ###
            phase_i = compute_phase_corrections(St[i])   # f(i)
        
            ### Correct the overlap matrix ###
            for a in xrange(nstates):   
                for b in xrange(nstates): 
                    fab = cum_phase.get(b) * cum_phase.get(b).conjugate() * phase_i.get(b).conjugate()
                    St[i].scale(a,b, fab)
        
            ### Update the cumulative phase correction factors ###
            for a in xrange(nstates):
                cum_phase.scale(a, 0, phase_i.get(a))
          


def sac_matrices(coeff):

    n_chi = len(coeff)
    n_phi = len(coeff[0])

    P2C = CMATRIX(n_phi, n_chi)
    for j in xrange(n_chi):
        for i in xrange(n_phi):
            P2C.set(i,j,coeff[j][i]*(1.0+0.0j) )


    # Normalize the Chi wavefunctions #
    norm = P2C.H() * P2C
    for i in xrange(n_chi):
        if norm.get(i,i).real > 0.0:
            P2C.scale(-1, i, 1.0/math.sqrt(norm.get(i,i).real) )
        else:
            print "Error in CHI normalizaiton: some combination gives zero norm\n"
            sys.exit(0)


    return P2C



def compute_Hvib(basis, St_ks, E_ks, dE, dt):
    """
    Basis - list of list of lists of integers
    St_ks - time overlap (CMATRIX) of the KS spin-orbitals
    E_ks - energies of KS spin-orbitals at the mid-point
    dE - energy corrections to the SD orbitals ("scissor operator")
    dt - the timestep for MD integrations

    Returns: The Vibronic Hamiltonian
    """

    St    = mapping.ovlp_mat_arb(basis, basis, St_ks) 
    H_el  = mapping.energy_mat_arb(basis, E_ks, dE)
    H_vib = H_el - (0.5j/dt)*(St-St.H())

    return H_vib




def run_namd(params):
    """
    The main procedure to run NA-MD calculations according to the PYXAID2 workflow

    === Required parameter keys: ===

    params["Phi_basis"]        [list of lists of ints] - define the Slater Determinants basis
    params["P2C"]              [list of lists of complex] - define the superpositions to SDs to get spin-adapted functions
    params["Phi_dE"]           [list of doubles] - define corrections of the SAC state energies
    params["T"]                [double] - temperature of nuclear/electronic dynamics [in K, default: 300.0]
    params["num_sh_traj"]      [int] - the number of stochastic surface hopping trajectories [default: 1000]
    params["sh_method"]        [int] - how to compute surface hopping probabilities [default: 1]
                               Options:
                               0 - MSSH
                               1 - FSSH
    params["do_collapse"]      [int] - decoherence method [default: 1]
                               Options:
                               0 - no decoherence
                               1 - ID-A
    params["dt"]               [double] - nuclear dynamics integration time step [in a.u. of time, default: 41.0]

    params["Boltz_opt"]        [int] - option to select a probability of hopping acceptance [default: 3]
                               Options:
                               0 - all proposed hops are accepted - no rejection based on energies
                               1 - proposed hops are accepted with exp(-E/kT) probability - the old (hence the default approach)
                               2 - proposed hops are accepted with the probability derived from Maxwell-Boltzmann distribution - more rigorous
                               3 - generalization of "1", but actually it should be changed in case there are many degenerate levels
    params["init_Chi"]         [int] - index of the initial spin-adapted state [default: 0]
    params["outfile"]          [string] - the name of the file where to print populations and energies of states [default: "_out.txt"]    


    === Required by the get_data() ===

    params["norbitals"]        [int] - how many lines/columns in the file - the total number of spin-orbitals
    params["active_space"]     [list of ints] - which orbitals we care about (indexing starts with 0)
    params["S_re_prefix"]      [string] - prefixes of the files with real part of the MO overlaps at time t
    params["S_re_suffix"]      [string] - suffixes of the files with real part of the MO overlaps at time t
    params["S_im_prefix"]      [string] - prefixes of the files with imaginary part of the MO overlaps at time t
    params["S_im_suffix"]      [string] - suffixes of the files with imaginary part of the MO overlaps at time t
    params["St_re_prefix"]     [string] - prefixes of the files with real part of the MO overlaps at times t and t' = t+dt
    params["St_re_suffix"]     [string] - suffixes of the files with real part of the MO overlaps at times t and t' = t+dt
    params["St_im_prefix"]     [string] - prefixes of the files with imaginary part of the MO overlaps at times t and t' = t+dt
    params["St_im_suffix"]     [string] - suffixes of the files with imaginary part of the MO overlaps at times t and t' = t+dt
    params["E_re_prefix"]      [string] - prefixes of the files with real part of the MO energies at time t
    params["E_re_suffix"]      [string] - suffixes of the files with real part of the MO energies  at time t
    params["E_im_prefix"]      [string] - prefixes of the files with imaginary part of the MO energies at time t
    params["E_im_suffix"]      [string] - suffixes of the files with imaginary part of the MO energies  at time t
    params["nsteps"]           [int] - how many files to read, starting from index 0
    params["is_pyxaid_format"] [Boolean] - whether the input in the old PYXAID format (just orbital space matrices)

     
    """

    critical_params = [ "P2C", "Phi_basis", "Phi_dE", ]
    default_params = { "T":300.0, "num_sh_traj":1000, 
                       "sh_method":1, "do_collapse":1, "dt":41.0, "Boltz_opt":3,
                       "init_Chi":0, "outfile":"_out.txt" }
    comn.check_input(params, default_params, critical_params)


    rnd = Random()
    T = params["T"]
    kT = units.kB * T
    num_sh_traj = params["num_sh_traj"]
    sh_method = params["sh_method"]
    do_collapse = params["do_collapse"]
    dt = params["dt"]
    bolt_opt = params["Boltz_opt"]

    """
    1. Read in the "elementary" overlaps and energies in the basis of KS orbitals
    2. Apply state reordering to KS
    3. Apply phase correction to KS
    4. Construct the Hvib in the basis of Slater determinants
    5. Convert the Hvib to the basis of symmery-adapted configurations (SAC)
    """

    P2C = sac_matrices(params["P2C"])
    
    S_dia_ks, St_dia_ks, E_dia_ks = get_data(params)  
    apply_state_reordering(St_dia_ks, E_dia_ks, params)    
    apply_phase_correction(St_dia_ks, params)

    nsteps = len(St_dia_ks)
    nstates = len(params["P2C"])
    H_vib = []
    for i in xrange(nsteps):
        # Hvib in the basis of SDs
        hvib = compute_Hvib(params["Phi_basis"], St_dia_ks[i], E_dia_ks[i], params["Phi_dE"], dt) 

        # SAC
        H_vib.append( P2C.H() * hvib * P2C )



    #========== Compute decoherence times  ===============
    tau, tau_inv = comn.decoherence_times(H_vib, 1)



    #========== Initialize the wavefunction amplitudes ===============
    init_Chi = params["init_Chi"]    

    # TD-SE coefficients
    Coeff_Chi = [];  
    Coeff_Phi = [];

    for tr in xrange(num_sh_traj):

        Coeff_Chi.append(CMATRIX(nstates, 1)); 
        Coeff_Chi[tr].set(init_Chi, 1.0, 0.0)
        Coeff_Phi.append(P2C*Coeff_Chi[tr]) 

        if tr==0:
            print "Coeff_Chi = "; Coeff_Chi[tr].show_matrix()
            print "Coeff_Phi = "; Coeff_Phi[tr].show_matrix()


    #========== Initialize the SE density matrices ===============
    # Density matrices
    denmat_Chi_se = tsh.amplitudes2denmat(Coeff_Chi)
    denmat_Phi_se = tsh.amplitudes2denmat(Coeff_Phi)

    print "denmat_Chi_se = "; denmat_Chi_se[0].show_matrix()
    print "denmat_Phi_se = "; denmat_Phi_se[0].show_matrix()

    #========== Initialize the state for adiabatic sh ===============
    istate  = []  # in the Chi basis (referring to spin-adaapted states)
    for tr in xrange(num_sh_traj):
        prob = tsh.denmat2prob(denmat_Chi_se[tr])
        st = tsh.set_random_state(prob, rnd.uniform(0.0, 1.0)) 
        istate.append(st)

    #========== Initialize output files ===============
    f = open("_pop_Chi_se.txt","w");   f.close()
    f = open("_pop_Chi_sh.txt","w");   f.close()
    f = open("_pop_Phi_se.txt","w");   f.close()
    f = open("_pop_Phi_sh.txt","w");   f.close()


    #=============== Entering Dynamics !! =========================  
    #=============== Now handle the electronic structure ==========
    for i in xrange(nsteps):

        pops = tsh.compute_sh_statistics(nstates, istate)
        comn.printout(i*params["dt"], pops, H_vib[i], params["outfile"])

        #============== TD-SE and surface hopping =================
        for tr in xrange(num_sh_traj):

            # Coherent evolution for Chi
            propagate_electronic(dt, Coeff_Chi[tr], H_vib[i])  # propagate in the diabatic basis (Chi)

            ksi  = rnd.uniform(0.0, 1.0);
            ksi2 = rnd.uniform(0.0, 1.0)

            # Surface hopping in Chi basis
            istate[tr] = tsh.hopping(Coeff_Chi[tr], H_vib[i], istate[tr], sh_method, do_collapse, ksi, ksi2, dt, T, bolt_opt)
            Coeff_Phi[tr] = P2C*Coeff_Chi[tr]


        #============== Analysis of Dynamics  =====================
        # Update SE-derived density matrices
        denmat_Chi_se = tsh.amplitudes2denmat(Coeff_Chi)
        denmat_Phi_se = tsh.amplitudes2denmat(Coeff_Phi)

        # Use SE-derived density matrices to make SH-derived density matrices
        denmat_Chi_sh = []
        denmat_Phi_sh = []
        for tr in xrange(num_sh_traj):

            denmat_Chi_sh.append(CMATRIX(nstates, nstates))
            denmat_Chi_sh[tr].set(istate[tr],istate[tr], 1.0, 0.0)
            denmat_Phi_sh.append( P2C*denmat_Chi_sh[tr]*P2C.H() )

        ave_pop_Chi_sh, ave_pop_Chi_se = tsh.ave_pop(denmat_Chi_sh, denmat_Chi_se)
        ave_pop_Phi_sh, ave_pop_Phi_se = tsh.ave_pop(denmat_Phi_sh, denmat_Phi_se)

        #===================== Print out ==============================
        add_printout(i, ave_pop_Chi_sh, "_pop_Chi_sh.txt")
        add_printout(i, ave_pop_Chi_se, "_pop_Chi_se.txt")

        add_printout(i, ave_pop_Phi_sh, "_pop_Phi_sh.txt")
        add_printout(i, ave_pop_Phi_se, "_pop_Phi_se.txt")


