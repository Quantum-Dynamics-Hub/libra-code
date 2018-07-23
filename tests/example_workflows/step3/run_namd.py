#*********************************************************************************
#* Copyright (C) 2017 Alexey V. Akimov, Wei Li
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

from libra_py import *
sys.path.insert(1,"/projects/academic/alexeyak/brendan/libra/libra-code/workflows/pyxaid2")
import mapping
#from PYXAID2 import *  #for test purpose, I run this file in a different folder, so this line is needed


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




def Chi2Phi(Chi_basis):
    """
    Chi - is the input of the spin-adapted SDs in the format:
    
    [ [  [c0, SD0  ],  [c1, SD1  ],  [c2, SD2  ],     ...     ],     
      [  [c0, SD0' ],  [c1, SD1' ],  [c2, SD2' ],     ...     ],
      ...
    ]


    |Chi> = |Phi> * T

    """

    nst_dia_sac = len(Chi_basis)  # length of the spin-adapted basis
    
    # First, lets find the unique determinants and make them ordered
    Phi_basis = []  # spin-diabatic SDs

    for i in xrange(nst_dia_sac):
        sz = len(Chi_basis[i])

        for j in xrange(sz):
            sd = sorted(Chi_basis[i][j][1])
 
            if sd not in Phi_basis:
                Phi_basis.append(sd)

    nst_dia = len(Phi_basis)    

    # Now determine the input-based superposition coefficients
    T = CMATRIX(nst_dia, nst_dia_sac)

    for i in xrange(nst_dia_sac):
        sz = len(Chi_basis[i])

        for j in xrange(sz):
            sd = sorted(Chi_basis[i][j][1])

            k = Phi_basis.index(sd)

            T.set(k,i, Chi_basis[i][j][0]) 


    return Phi_basis, T

    
def get_matrix(nrows, ncols, i, prefix_re, suffix_re, prefix_im, suffix_im):

    filename_re = prefix_re + str(i) + suffix_re
    filename_im = prefix_im + str(i) + suffix_im

    print "\nfilename_re\n"
    print filename_re
    print "\nfilename_im\n"
    print filename_im

    x_re = MATRIX(nrows, ncols); x_re.Load_Matrix_From_File(filename_re)
    x_im = MATRIX(nrows, ncols); x_im.Load_Matrix_From_File(filename_im)


    return CMATRIX(x_re, x_im)


def compute_Hvib_dia(Phi_basis, St_dia_ks, E_dia_ks, dE, dt ):
    """
    Psi_basis - spin-adiabatic SDs (list of lists of integers)

    St_adi_ks - time overlap (CMATRIX) of the spin-adiabatic KS orbitals
    E_adi_ks - energies of spin-adiabatic KS orbitals at the mid-point
    dE - energy corrections to the spin-adiabatic SD orbitals ("scissor operator")
    dt - the timestep for MD integrations

    """

    nst_dia = len(Phi_basis)

    Stdia = mapping.ovlp_mat_dd(Phi_basis, Phi_basis, St_dia_ks) # nst_adi x nst_adi, one of the terms in Eq. 15
    H_el = mapping.energy_mat_dd(Phi_basis, E_dia_ks, dE)
    H_vib = H_el - (0.5j/dt)*(Stdia-Stdia.H())
    return H_vib


def compute_Hvib_adi(Psi_basis, St_adi_ks, E_adi_ks, dE, dt ):
    """
    Psi_basis - spin-adiabatic SDs (list of lists of integers)

    St_adi_ks - time overlap (CMATRIX) of the spin-adiabatic KS orbitals
    E_adi_ks - energies of spin-adiabatic KS orbitals at the mid-point
    dE - energy corrections to the spin-adiabatic SD orbitals ("scissor operator")
    dt - the timestep for MD integrations

    """

    nst_adi = len(Psi_basis)

    Stadi = mapping.ovlp_mat_aa(Psi_basis, Psi_basis, St_adi_ks) # nst_adi x nst_adi, one of the terms in Eq. 15
    H_el = mapping.energy_mat_aa(Psi_basis, E_adi_ks, dE)
    H_vib = H_el - (0.5j/dt)*(Stdia-Stdia.H())
    return H_vib



def compute_L(Phi_basis, Psi_basis, S_adi_ks, S_dia_adi_ks, S_dia_ks):
    # 
    # Phi_basis - spin-diabatic SDs (list of lists of integers)
    # Psi_basis - spin-adiabatic SDs (list of lists of integers)
    # S_adi_ks - overlaps of spin-adiabatic KS orbitals
    # S_dia_adi_ks - overlaps of spin-diabatic and spin-adiabatic KS orbitals
    # S_dia_ks - overlaps of spin-diabatic KS orbitals


    L = mapping.ovlp_mat_da(Phi_basis, Psi_basis, S_dia_adi_ks)       # nst_dia x nst_adi, Eq. 25
    
    # These may be needed later: for orthogonalizaiton!
    #S_adi = mapping.ovlp_mat_aa(Psi_basis, Psi_basis, S_adi_ks)   # nst_adi x nst_adi
    #S_dia = mapping.ovlp_mat_dd(Phi_basis, Phi_basis, S_dia_ks)   # nst_dia x nst_dia

    return L  #, S_adi, S_dia




def run_namd(params):

    rnd = Random()
    T = params["T"]
    init_time = params["init_time"]  # integer
    dt = params["dt"]
    len_traj = params["len_traj"]
    num_sh_traj = params["num_sh_traj"]
    sh_method = params["sh_method"]
    do_collapse = params["do_collapse"]

    psi_dia_ks = params["psi_dia_ks"];   nst_dia_ks = 2*len(psi_dia_ks)

    Chi_basis  = params["Chi_basis"];    nst_dia_sac = len(Chi_basis)
    Phi_basis, C2P = Chi2Phi(Chi_basis); nst_dia = len(Phi_basis) # C2P is the matrix T from Eq. 3, (nst_dia x nst_dia_sac)

    print "\nPrinting Chi_basis"
    print Chi_basis
    print "\nPrinting Phi_basis"
    print Phi_basis
    print "\nPrinting C2P"
    C2P.show_matrix()

    # ------------------read and store the projection and energies------------------
    H_vib = []

    for i in range(0,len_traj):

        ##############################################################################
        # Read in the "elementary" overlaps and energies - in the basis of KS orbitals
        ##############################################################################       

        re_pr = params["E_dia_ks_re_prefix"] 
        re_sf = params["E_dia_ks_re_suffix"] 
        im_pr = params["E_dia_ks_im_prefix"] 
        im_sf = params["E_dia_ks_im_suffix"] 
        E_dia_ks = get_matrix(nst_dia_ks, nst_dia_ks, i, re_pr, re_sf, im_pr, im_sf )

        re_pr = params["S_dia_ks_re_prefix"] 
        re_sf = params["S_dia_ks_re_suffix"] 
        im_pr = params["S_dia_ks_im_prefix"] 
        im_sf = params["S_dia_ks_im_suffix"] 
        S_dia_ks = get_matrix(nst_dia_ks, nst_dia_ks, i, re_pr, re_sf, im_pr, im_sf )

        re_pr = params["St_dia_ks_re_prefix"] 
        re_sf = params["St_dia_ks_re_suffix"] 
        im_pr = params["St_dia_ks_im_prefix"] 
        im_sf = params["St_dia_ks_im_suffix"] 
        St_dia_ks = get_matrix(nst_dia_ks, nst_dia_ks, i, re_pr, re_sf, im_pr, im_sf )

        # Printing what we just extracted
        if i == 0:

            print "\nE_dia_ks = "
            E_dia_ks.real().show_matrix()

            print "\nS_dia_ks = "
            S_dia_ks.real().show_matrix()  

            print "\nSt_dia_ks = "
            St_dia_ks.real().show_matrix()

        # Now, let's compute the vibronic Hamiltonian in the SD basis
        Hvib_dia = compute_Hvib_dia(Phi_basis, St_dia_ks, E_dia_ks, params["Phi_dE"], params["dt"] )

        # Print the vibronic Hamiltonian in the SD basis 
        print "\nHvib_dia.imag()\n"
        Hvib_dia.imag().show_matrix()

        print "\nHvib_dia.real()\n"
        Hvib_dia.real().show_matrix()
        #"""
   

        H_vib.append ( compute_Hvib_dia(Phi_basis, St_dia_ks, E_dia_ks, params["Phi_dE"], params["dt"] ) )
       

    init_Chi = params["init_Chi"]    
    #========== Initialize the wavefunction amplitudes ===============
    # TD-SE coefficients
    Coeff_Chi = [];  
    Coeff_Phi = [];

    for tr in xrange(num_sh_traj):

        Coeff_Chi.append(CMATRIX(nst_dia_sac, 1)); 
        Coeff_Chi[tr].set(init_Chi, 1.0, 0.0)
        Coeff_Phi.append(C2P * Coeff_Chi[tr]) 

        if tr==0:
            print "Coeff_Chi = "; Coeff_Chi[tr].show_matrix()
            print "Coeff_Phi = "; Coeff_Phi[tr].show_matrix()

    #========== Initialize the SE density matrices ===============
    # Density matrices

    denmat_Chi_se = tsh.amplitudes2denmat(Coeff_Chi)
    denmat_Phi_se = tsh.amplitudes2denmat(Coeff_Phi)

    print "denmat_Chi_se = "; denmat_Chi_se[0].show_matrix()
    print "denmat_Phi_se = "; denmat_Phi_se[0].show_matrix()

    #========== Initialize the SH density matrices ===============
    denmat_Chi_sh = []
    denmat_Phi_sh = []

    # initialize the state for adiabatic sh
    istate  = []  # in the Psi basis (adiabatic SDs)

    for tr in xrange(num_sh_traj):

        prob = tsh.denmat2prob(denmat_Chi_se[tr])
        st = tsh.set_random_state(prob, rnd.uniform(0.0, 1.0)) 
        istate.append(st)

        denmat_Phi_sh.append(CMATRIX(nst_dia, nst_dia))
        denmat_Phi_sh[tr].set(st,st, 1.0, 0.0)

        #denmat_Phi_sh.append( L[init_time] * denmat_Psi_sh[tr] * L[init_time].H())
        denmat_Chi_sh.append( C2P.H() * denmat_Phi_sh[tr] * C2P )


    print "denmat_Chi_sh = "; denmat_Chi_sh[0].show_matrix()
    print "denmat_Phi_sh = "; denmat_Phi_sh[0].show_matrix()

    #f = open("_pop_Chi_se.txt","w");   f.close()
    #f = open("_pop_Chi_sh.txt","w");   f.close()
    f = open("_pop_Phi_se.txt","w");   f.close()
    f = open("_pop_Phi_sh.txt","w");   f.close()


    print "\nEntering SH Dynamics\n"
    print "istate = ", istate
    
    #=============== Now handle the electronic structure ==========================
    for i in xrange(init_time, init_time + len_traj):

        #============== TD-SE and surface hopping =================

        for tr in xrange(num_sh_traj):

            # Coherent evolution
            propagate_electronic(dt, Coeff_Phi[tr], H_vib[i])  # propagate in the diabatic basis (Chi)

            ksi = rnd.uniform(0.0, 1.0);
            ksi2 = rnd.uniform(0.0, 1.0)

            # Surface hopping
            istate[tr] = tsh.hopping(Coeff_Phi[tr], H_vib[i], istate[tr], sh_method, do_collapse, ksi, ksi2, dt, T)
    
  
        #============== Projections and analysis  =================

        # Update SE-derived density matrices
        denmat_Chi_se = tsh.amplitudes2denmat(Coeff_Chi)
        denmat_Phi_se = tsh.amplitudes2denmat(Coeff_Phi)

        # Update SH-derived density matrices
        for tr in xrange(num_sh_traj):
            denmat_Phi_sh[tr].set(istate[tr],istate[tr], 1.0, 0.0)   # see Eqs. 24-26
            denmat_Chi_sh[tr] = C2P.H() * denmat_Phi_sh[tr] * C2P     # see Eqs. 24-26


        ave_pop_Chi_sh, ave_pop_Chi_se = tsh.ave_pop(denmat_Chi_sh, denmat_Chi_se)
        ave_pop_Phi_sh, ave_pop_Phi_se = tsh.ave_pop(denmat_Phi_sh, denmat_Phi_se)


        #===================== Print out ==============================

        #add_printout(i, ave_pop_Chi_sh, "_pop_Chi_sh.txt")
        #add_printout(i, ave_pop_Chi_se, "_pop_Chi_se.txt")

        add_printout(i, ave_pop_Phi_sh, "_pop_Phi_sh.txt")
        add_printout(i, ave_pop_Phi_se, "_pop_Phi_se.txt")


