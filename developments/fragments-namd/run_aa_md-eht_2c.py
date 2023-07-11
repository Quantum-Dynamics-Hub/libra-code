#*********************************************************************************
#* Copyright (C) 2016 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
#
# Here, we switch back to orthogonal orbitals in the frontier subspace
# We use local adiabatic states!
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


def hop_py(initstate, g, ksi):

    nstates = g.num_of_cols
    finstate = initstate;

    left, right = 0.0, 0.0

    for i in xrange(nstates):
      if i==0:
          left = 0.0
          right = g.get(initstate,i)
      else:
          left = right
          right = right + g.get(initstate,i)
 
      if((left<ksi) and (ksi<=right)):
          finstate = i

    return finstate


def set_random_state(prob, ksi):

    nstates = len(prob)
    finstate = 0;

    left, right = 0.0, 0.0

    for i in xrange(nstates):
      if i==0:
          left = 0.0
          right = prob[i]
      else:
          left = right
          right = right + prob[i]
 
      if((left<ksi) and (ksi<=right)):
          finstate = i

    return finstate



def ida_py(Coeff, old_st, new_st, E_old, E_new, T, ksi, do_collapse):

    kb = 3.166811429e-6  # Hartree/K
    res = old_st;  C = CMATRIX(Coeff)
    dE = (E_new - E_old);   boltz_f = 1.0

    if dE>0.0:
        argg = dE/(kb*T)
        if argg > 50.0:
            boltz_f = 0.0
        else:
            boltz_f = math.exp(-argg)

        if ksi<boltz_f:
            res = new_st  # accepted hop
            
            # Collapse the wavefunction to the new state 
            if do_collapse:
                C *= 0.0; C.set(new_st, 1.0+0.0j)
        else:
            # Unsuccessful hop - collapse wfc back to the original state^M
            if do_collapse:
                C *= 0.0; C.set(old_st, 1.0+0.0j)
    else:
        res = new_st

    return res , C


def compute_S_H_el(opt, sub_ham, ham, active_orb, print_amt):
    ##
    # Compute S and H matrices for the chosen space of active fragment MOs (diabatic MOs)
    #
    # \param[in] opt - option controlling which type of calculations is to perform:
    #            opt = 0  - disregard electronic couplings
    #            opt = 1  - use the Martinez formula with overlaps - faster, but approximate
    #            opt = 2  - use AO Fock matrix transformation - slower, but more rigorous
    #
    # \param[in] el - A list of Electronic structure objects - one for each fragment  - includes also eigenvalues/eigenvectors
    # \param[in] el_all - A single Electronic structure object containing the overlap and Fock matrix for the entire system
    #                     we do not expect eigenvalues/eigenvectors (so don't compute them outside - that would save a lot of time)
    #            Set to "None", if the "opt" value chooses an algorithm that does not use this variable.
    # \param[in] basis - the list of AO lists for all fragments
    # \param[in] active_orb The list containing the nfrags lists of the orbital indices from each configuration    
    # \param[in] print_amt The flag controlling the amount of printout


    t = Timer();  t.start()

    nfrags = len(active_orb)     # the number of fragments

    nmo = 0                      # the total number of MOs in all fragments
    for fr in xrange(nfrags):
        nmo = nmo + len(active_orb[fr])


    # Result matrices
    S_el = MATRIX(nmo, nmo)
    H_el = MATRIX(nmo, nmo)
    
    dist2 = 10000.0



    if opt==0 or opt==1:

        start_f1 = 0
        for f1 in xrange(nfrags):        
            el1 = sub_ham[f1].get_electronic_structure()
            bas1 = sub_ham[f1].basis_ao
            nmo1 = len(active_orb[f1])  # the number of MOs in the fragment f1
            E1 = MATRIX(nmo1,nmo1)
            pop_submatrix( el1.get_E_alp(), E1, active_orb[f1], active_orb[f1] )  # Energies
        
            start_f2 = 0
            for f2 in xrange(nfrags):
                el2 = sub_ham[f2].get_electronic_structure()
                bas2 = sub_ham[f2].basis_ao
                nmo2 = len(active_orb[f2])  # the number of MOs in the fragment f2
                E2 = MATRIX(nmo2,nmo2)
                pop_submatrix( el2.get_E_alp(), E2, active_orb[f2], active_orb[f2] )  # Energies
        
                smo = MATRIX(nmo1, nmo2) # overlap of MOs of fragment f1 and f2
                hmo = MATRIX(nmo1, nmo2) # electronic coupling of MOs of fragment f1 and f2

                MO_overlap(smo, bas1, bas2, el1.get_C_alp(), el2.get_C_alp(), Py2Cpp_int(active_orb[f1]), Py2Cpp_int(active_orb[f2]), dist2)  # fr1 x fr2 block

                if opt==0:
                    if f2==f1:
                        for i1 in xrange(nmo1):
                            hmo.set(i1,i1, E1.get(i1,i1) )  
                elif opt==1:        
                    for i1 in xrange(nmo1):
                        for i2 in xrange(nmo2):                  
                            hmo.set(i1,i2, 0.5*(E1.get(i1,i1) + E2.get(i2,i2))*smo.get(i1,i2) )
        
                push_submatrix( S_el, smo, range(start_f1, start_f1 + nmo1), range(start_f2, start_f2 + nmo2))
                push_submatrix( H_el, hmo, range(start_f1, start_f1 + nmo1), range(start_f2, start_f2 + nmo2))

        
                del el2; del bas2; del E2; del smo; del hmo;
                start_f2 = start_f2 + nmo2

            del el1; del bas1; del E1;
            start_f1 = start_f1 + nmo1


    elif opt==2:

        S = ham.get_electronic_structure().get_Sao()
        H = ham.get_electronic_structure().get_Hao()

        start_f1 = 0
        start_f1a = 0
        for f1 in xrange(nfrags): 
            el1 = sub_ham[f1].get_electronic_structure()
            bas1 = sub_ham[f1].basis_ao
            nmo1 = len(active_orb[f1])  # the number of MOs in the fragment f1
            nao1 = el1.Norb             # the number of AOs in the fragment f2
            C1 = MATRIX(nao1,nmo1)
            pop_submatrix( el1.get_C_alp(), C1, range(0,nao1),  active_orb[f1] )  # MO-LCAO
        
            start_f2 = 0
            start_f2a = 0
            for f2 in xrange(nfrags):
                el2 = sub_ham[f2].get_electronic_structure()
                bas2 = sub_ham[f2].basis_ao
                nmo2 = len(active_orb[f2])  # the number of MOs in the fragment f2
                nao2 = el2.Norb             # the number of AOs in the fragment f2
                C2 = MATRIX(nao2,nmo2)
                pop_submatrix( el2.get_C_alp(), C2, range(0,nao2),  active_orb[f2] )  # MO-LCAO


                # Extract matrices in AO basis
                s = MATRIX(nao1,nao2)
                h = MATRIX(nao1,nao2)

                pop_submatrix( S, s, range(start_f1a, start_f1a + nao1), range(start_f2a, start_f2a + nao2))
                pop_submatrix( H, h, range(start_f1a, start_f1a + nao1), range(start_f2a, start_f2a + nao2))


                # Transformation                
                x = C1.T() * (s * C2)   # Identity matrix for the diagonal blocks  (nmo1 x nao1) x (nao1 x nao2) x (nao2 x nmo2) 
                y = C1.T() * (h * C2)   # Eigenvalues for the the diagonal blocks  (nmo1 x nao1) x (nao1 x nao2) x (nao2 x nmo2)


                # Push matrices in MO basis
                push_submatrix( S_el, x, range(start_f1, start_f1 + nmo1), range(start_f2, start_f2 + nmo2))
                push_submatrix( H_el, y, range(start_f1, start_f1 + nmo1), range(start_f2, start_f2 + nmo2))


                del el2; del bas2; del C2; del s; del h; del x; del y
        
                start_f2 = start_f2 + nmo2
                start_f2a = start_f2a + nao2

            del el1; del bas1; del C1

            start_f1 = start_f1 + nmo1 
            start_f1a = start_f1a + nao1

        del S; del H; 

    if print_amt==1:
        print "Final overlap";     S_el.show_matrix()
        print "Final Hamiltonian"; H_el.show_matrix()

    t.stop();  print "Overlap and electronic Hamiltonian = ", t.show()

    return S_el, H_el



def mo_fmo_overlap(sub_ham, active_orb, ham, active_orb_adi):
# This function computes the matrix:
#  S = <A-MO|A-FMO> 
    t = Timer();  t.start()

    print "In mo_fmo_overlap..."

    nfrags = len(active_orb)     # the number of fragments

    nmo = 0                      # the total number of MOs in all fragments
    for fr in xrange(nfrags):
        print fr, active_orb[fr]
        nmo = nmo + len(active_orb[fr])

    nmo_adi = len(active_orb_adi) # the number of adiabatic MOs
    print active_orb_adi


    el_str_adi = ham.get_electronic_structure()

    nao = el_str_adi.Norb
    C_mo = el_str_adi.get_C_alp()   #  nao x nao
    S_ao = el_str_adi.get_Sao()     #  nao x nao


    c_mo = MATRIX(nao, nmo_adi)
    pop_submatrix( C_mo, c_mo, range(0,nao),  active_orb_adi )  


    C_fmo = MATRIX(nao, nmo)  # all FMO orbitals in one place

    start_f1 = 0
    start_f1a = 0

    for f1 in xrange(nfrags): 
 
        el1 = sub_ham[f1].get_electronic_structure()
        nmo1 = len(active_orb[f1])  # the number of MOs in the fragment f1
        nao1 = el1.Norb             # the number of AOs in the fragment f1

        print f1, nmo1, nao1

        C1 = MATRIX(nao1,nmo1)
        pop_submatrix( el1.get_C_alp(), C1, range(0,nao1),  active_orb[f1] )  # MO-LCAO
        push_submatrix( C_fmo, C1, range(start_f1a, start_f1a + nao1), range(start_f1, start_f1 + nmo1) )

#        print "D-FMO = "; el1.get_C_alp().show_matrix()

        del el1; del C1

        start_f1 = start_f1 + nmo1 
        start_f1a = start_f1a + nao1

#    print "C_fmo = \n"; C_fmo.show_matrix()
#    print "c_mo = \n"; c_mo.show_matrix()

    # Result matrices
    S = c_mo.T() * S_ao * C_fmo  #  (nmo_adi x nao) x (nao x nao) x (nao x nmo) = nmo_adi x nmo

    if 0:
        print "\nVerification in mo_fmo_overlap\n"
        print "<A-MO|A-MO> = "; (c_mo.T() * S_ao * c_mo).show_matrix()
        print "<D-FMO|D-FMO> = "; (C_fmo.T() * S_ao * C_fmo).show_matrix()
        print "S.T() * S = "; (S.T()*S).show_matrix()

    del C_mo; del c_mo; del C_fmo; del S_ao;  del el_str_adi

    t.stop();  print "MO-FMO overlaps time = ", t.show()

    return CMATRIX(S)  # Complex



def compute_NAC(sub_ham_old, sub_ham_cur, active_orb):
    ##
    # Compute NAC matrix for the entire space of active fragment MOs
    #
    # \param[in] el_old - A list of Electronic structure objects - one for each fragment  - includes also eigenvalues/eigenvectors at t-dt
    # \param[in] el_cur - A list of Electronic structure objects - one for each fragment  - includes also eigenvalues/eigenvectors at t 
    # \param[in] dt - the timestep between the two configurations
    # \param[in] active_orb The list containing the nfrags lists of the orbital indices from each configuration    


    t = Timer();  t.start()

    nfrags = len(active_orb)     # the number of fragments

    nmo = 0                      # the total number of MOs in all fragments
    for fr in xrange(nfrags):
        nmo = nmo + len(active_orb[fr])


    # Result matrices
    D = MATRIX(nmo, nmo);    dist2 = 10000.0


    start_f1 = 0
    for f1 in xrange(nfrags):
        nmo1 = len(active_orb[f1])  # the number of MOs in the fragment f1

        start_f2 = 0
        for f2 in xrange(nfrags):   
            nmo2 = len(active_orb[f2])  # the number of MOs in the fragment f2


            Dmo = MATRIX(nmo1, nmo2)

            bas_old = sub_ham_old[f1].basis_ao
            bas_cur = sub_ham_cur[f2].basis_ao
            el_old = sub_ham_old[f1].get_electronic_structure()
            el_cur = sub_ham_cur[f2].get_electronic_structure()

            MO_overlap(Dmo, bas_old, bas_cur, el_old.get_C_alp(), el_cur.get_C_alp(), Py2Cpp_int(active_orb[f1]), Py2Cpp_int(active_orb[f2]),dist2)  # fr1 x fr2 block
            push_submatrix( D, Dmo, range(start_f1, start_f1 + nmo1), range(start_f2, start_f2 + nmo2))


            del Dmo; del bas_old; del bas_cur; del el_old; del el_cur;
            start_f2 = start_f2 + nmo2
        start_f1 = start_f1 + nmo1

#    D = (0.5/dt)*(D - D.T()) # the Smo_dot term will cancel out during solution of TD-SE !!! 
   
    t.stop();  print "Nonadiabatic coupligs time = ", t.show()

    return D




def update_Hvib_L(opt, sub_ham_old, sub_ham_cur, ham_old, ham_cur, active_orb, dt, print_amt):
    ##
    # updates vibronic Hamiltonian
    # \param[in, out] Hvib  vibronic Hamiltonian - the CMATRIX object N x N 
    # \param[in, out] Smo orbital overlap - the CMATRIX object N x N
    # \param[in] active_orb The list of the length nfrags, containing the lists of the orbital indices from each fragment
    # must be: N = len(active_orb[0]) + len(active_orb[1]) + ... + len(active_orb[nfrags-1])
    # \param[in] mo_alp_old The MO information at the time t-dt
    # \param[in] mo_alp The MO information at the time t
    # \param[in] print_amt The flag controlling whether to print out the information

    #.get_electronic_structure()

    t = Timer();  t.start()

    Smo_old, Hel_old = compute_S_H_el(opt, sub_ham_old, ham_old, active_orb, print_amt)
    Smo_cur, Hel_cur = compute_S_H_el(opt, sub_ham_cur, ham_cur, active_orb, print_amt)
    Dmo = compute_NAC(sub_ham_old, sub_ham_cur, active_orb)
    Dmo = (0.5/dt)*(Dmo - Dmo.T()) # the Smo_dot term will cancel out during solution of TD-SE !!! 

    Smo = 0.5*(Smo_cur + Smo_old) 
    Hel = 0.5*(Hel_cur + Hel_old) 

        
    Hvib = CMATRIX(Hel, -1.0*Dmo)  # so Hvib is now Hermitian

    sz = Smo.num_of_cols  
    S    = CMATRIX(Smo)
    S_half = CMATRIX(sz,sz)
    S_i_half = CMATRIX(sz,sz)

    sqrt_matrix(S, S_half, S_i_half)

    if 0:
        print "In update Hvib_L:\n"
        print "S_i_half * S_i_half * S = "; (S_i_half * S_i_half * S ).show_matrix()


    # Now transform the Hamiltonian
    Heff = CMATRIX(sz,sz) 
    Heff = S_i_half.H() * Hvib * S_i_half

    if print_amt==1:
        print "Original Hamiltonian = "; Hvib.show_matrix()
        print "Effective Hamiltonian = "; Heff.show_matrix()
        print "Transformation matrix S^{-1/2} = "; S_i_half.show_matrix()
   
    del Smo_old; del Hel_old; del Smo_cur; del Hel_cur; del Dmo; del Smo; del Hel; del S; del Hvib; 

    t.stop();  print "update_Hvib_L time = ", t.show()

    return Heff, S_half, S_i_half  # all are CMATRIX



def update_Hvib_A(opt, sub_ham_old, sub_ham_cur, ham_old, ham_cur, active_orb, dt, print_amt):
    ##
    # updates vibronic Hamiltonian
    # \param[in, out] Hvib  vibronic Hamiltonian - the CMATRIX object N x N 
    # \param[in, out] Smo orbital overlap - the CMATRIX object N x N
    # \param[in] active_orb The list of the length nfrags, containing the lists of the orbital indices from each fragment
    # must be: N = len(active_orb[0]) + len(active_orb[1]) + ... + len(active_orb[nfrags-1])
    # \param[in] mo_alp_old The MO information at the time t-dt
    # \param[in] mo_alp The MO information at the time t
    # \param[in] print_amt The flag controlling whether to print out the information

    #.get_electronic_structure()

    t = Timer();  t.start()

    # Diabatic FMO:
    Smo_old, Hel_old = compute_S_H_el(opt, sub_ham_old, ham_old, active_orb, print_amt)
    Smo_cur, Hel_cur = compute_S_H_el(opt, sub_ham_cur, ham_cur, active_orb, print_amt)


    # Let us first diagonalize the diabatic Hamiltonians in the FMO basis at two timesteps
    # at time "t-dt"
    sz = Hel_old.num_of_cols;
    C_old = MATRIX(sz, sz);  E_old = MATRIX(sz, sz) 
    #solve_eigen_gen(sz, Hel_old, Smo_old, E_old, C_old)  # H * C = S * C * E  
    symm = 0 # symmetrize Hel_old and Smo_old
    solve_eigen(Hel_old, Smo_old, E_old, C_old, symm) # H * C = S * C * E

    # at time "t"
    C_cur = MATRIX(sz, sz);  E_cur = MATRIX(sz, sz) 
    #solve_eigen_gen(sz, Hel_cur, Smo_cur, E_cur, C_cur)  # H * C = S * C * E  
    solve_eigen(Hel_cur, Smo_cur, E_cur, C_cur, symm) # H * C = S * C * E
    
    Smo = Smo_cur # 0.5*(Smo_cur + Smo_old)   # overlap is the identity matrix
    Hel = E_cur #0.5*(E_cur + E_old)   # adiabatic Hamiltonian


    # Now compute NACs between adiabatic FMO
    Dmo = compute_NAC(sub_ham_old, sub_ham_cur, active_orb)
    Dmo_adi = MATRIX(sz,sz)
    Dmo_adi = C_old.T() * Dmo * C_cur
    Dmo_adi = (0.5/dt)*(Dmo_adi - Dmo_adi.T()) # the Smo_dot term will cancel out during solution of TD-SE !!! 
        
    Hvib = CMATRIX(Hel, -1.0*Dmo_adi)  # so Hvib is now Hermitian


    if print_amt==1:
        print "Vibronic Hamiltonian = "; Hvib.show_matrix()
        print "ADI-FMO = LC-DIA-FMO coefficients = "; C_cur.show_matrix()
  

    del E_old; del E_cur;
    del C_old; #del C_cur;
    del Smo_old; del Smo_cur;
    del Hel_old; del Hel_cur;
    del Dmo; del Smo; del Hel; del Dmo_adi;

    t.stop();  print "update_Hvib_A time = ", t.show()

    return Hvib, CMATRIX(C_cur)   # all are CMATRIX




def compute_sh_statistics(nstates, istate):

    num_sh_traj = len(istate)
    f = 1.0/float(num_sh_traj)

    coeff_sh = MATRIX(nstates, 1)

    for i in xrange(num_sh_traj):
        st = istate[i]
        coeff_sh.set(st, coeff_sh.get(st) + f)
 
    return coeff_sh
    


def print_mo_frag(active_orb, sub_ham_cur, syst, prefix):

    nfrags = len(active_orb) 
    # For fragments
    for fr in xrange(nfrags):
        prms = Control_Parameters()
        prms.nx_grid, prms.ny_grid, prms.nz_grid = 40, 40, 40
        prms.charge_density_prefix = prefix+"_fmo_"
        prms.orbs = Py2Cpp_int(active_orb[fr]) 
        el_str_fr = sub_ham_cur[fr].get_electronic_structure()
        charge_density( el_str_fr, syst, sub_ham_cur[fr].basis_ao, prms)



def print_mo_adi(active_orb_adi, ham_cur, syst, prefix):

    prms = Control_Parameters()
    prms.nx_grid, prms.ny_grid, prms.nz_grid = 40, 40, 40
    prms.charge_density_prefix = prefix+"_mo_"
    prms.orbs = Py2Cpp_int(active_orb_adi) 
    el_str_adi = ham_cur.get_electronic_structure()
    charge_density( el_str_adi, syst, ham_cur.basis_ao, prms)



def dum():
    if system=="Pc-C60-deca":
       if frag_mode==0: # no fragmentation
           # Start as a superposition of L+3 + L+4
           nr = (1.0/math.sqrt(2.0)) * (1.0+0.0j)
           Coeff.set(3, nr);  Coeff.set(4, nr);
           istate = [4]*num_sh_traj  # is the initial state for the given trajectory
       else:
           # Start as a superposition of L + L+1
           nr = (1.0/math.sqrt(2.0)) * (1.0+0.0j)
           if method==0:
               Coeff.set(0, nr);  Coeff.set(1, nr) 
               istate = [0]*num_sh_traj  # is the initial state for the given trajectory
           elif method==1:
               Coeff.set(nstates-2, nr);  Coeff.set(nstates-1, nr)  # because the diagonalization of diabatic FMOs will change ordering
               istate = [nstates-1]*num_sh_traj  # 



def rep_transform(bastyp, S_half, S_i_half, U_cur, S_mo_fmo, coeffs):
# This function performs transformations from one basis to another
# The function will transform coefficients from the selected representation into all other, provided all needed
# matrices are here. The 
# bastyp - selects the initial representation, for which we assume we already have the coefficents
# coeffs - a container of coefficients in different representations: 
#             0         1         2       3
# coeffs = [Coeff_D, Coeff_L, Coeff_A, Coeff], where

# Coeff_D - D-FMO - diabatic FMO basis
# Coeff_L - L - Lowdin-transformed diabatic FMO basis
# Coeff_A - A-FMO - adiabatic FMO basis
# Coeff - A-MO (MO) - adiabatic MO (overall system) basis

#    print "In rep_transform..."
    
    if bastyp==0:   # D-FMO
        if S_half!=None:            
            coeffs[1] = S_half * coeffs[0]          # D-FMO -> L
        if U_cur!=None and S_half!=None:
            S = S_half * S_half
            coeffs[2] = U_cur.H() * S * coeffs[0]   # D-FMO -> A-FMO
        if S_mo_fmo!=None:
            coeffs[3] = S_mo_fmo * coeffs[0]        # D-FMO -> A-MO
   
    elif bastyp==1:  # L
        if S_i_half!=None:
            coeffs[0] = S_i_half * coeffs[1]                    #  L -> D-FMO
        if U_cur!=None and S_i_half!=None:
            coeffs[2] = U_cur.H() * S_i_half * coeffs[1]       #  L -> A-FMO
        if S_mo_fmo!=None and S_i_half!=None:                                        
#            print "S_mo_fmo matrix = "; S_mo_fmo.show_matrix()
#            print "Transformation matrix = "; (S_mo_fmo * S_i_half).show_matrix()
#            print "coeffs[1] = "; coeffs[1].show_matrix()

            coeffs[3] = S_mo_fmo * S_i_half * coeffs[1]        #  L -> A-MO
   
    elif bastyp==2:   # A-FMO
        if U_cur!=None:
            coeffs[0] = U_cur * coeffs[2]                       #  A-FMO -> D-FMO
        if U_cur!=None and S_half!=None:
            coeffs[1] = S_half * U_cur * coeffs[2]              #  A-FMO -> L
        if S_mo_fmo!=None and U_cur!=None:
            coeffs[3] = S_mo_fmo * U_cur * coeffs[2]            #  A-FMO -> A-MO
   
    elif bastyp==3:   # A-MO
        if S_mo_fmo!=None and S_i_half!=None:                                         
            S_i = S_i_half * S_i_half
            coeffs[0] = (S_i * S_mo_fmo.T()) * coeffs[3]        #  A-MO -> D-FMO
        if S_mo_fmo!=None and S_half!=None:                    
            coeffs[1] = (S_half * S_mo_fmo.T()) * coeffs[3]     #  A-MO -> L
        if S_mo_fmo!=None and U_cur!=None:
            coeffs[2] = U_cur.H() * S_mo_fmo.T() * coeffs[3]    #  A-MO -> A-FMO



def populations(bastyp, S_half, S_i_half, U_cur, S_mo_fmo, coeffs):
# This function performs transformations from one basis to another
# The function will transform coefficients from the selected representation into all other, provided all needed
# matrices are here. The 
# bastyp - selects the initial representation, for which we assume we already have the coefficents
# coeffs - a container of coefficients in different representations: 
#             0         1         2       3
# coeffs = [Coeff_D, Coeff_L, Coeff_A, Coeff], where

# Coeff_D - D-FMO - diabatic FMO basis
# Coeff_L - L - Lowdin-transformed diabatic FMO basis
# Coeff_A - A-FMO - adiabatic FMO basis
# Coeff - A-MO (MO) - adiabatic MO (overall system) basis

    # First, transform the propagated variables into other representations
    rep_transform(bastyp, S_half, S_i_half, U_cur, S_mo_fmo, coeffs)    

    # Compute density matrix in all representations
    denmat = [None, None, None, None]

    if S_half!=None:
        S = S_half * S_half
        denmat[0] = S * coeffs[0] * coeffs[0].H() * S
    denmat[1] = coeffs[1] * coeffs[1].H()
    denmat[2] = coeffs[2] * coeffs[2].H()
    denmat[3] = coeffs[3] * coeffs[3].H()

    # Compute populations
    pops = []
    for rep in [0,1,2,3]:  # all representations
        sz = coeffs[rep].num_of_elems

        pop_r = [0.0]*sz
        for i in xrange(sz):
            pop_r[i] = denmat[rep].get(i,i).real

        pops.append(pop_r)

    return pops


def ave_Pop(nstates, nstates_adi, num_sh_traj, Pop_se, istate):

    ave_Pop_sh = [[0.0]*nstates, [0.0]*nstates, [0.0]*nstates, [0.0]*nstates_adi] 
    ave_Pop_se = [[0.0]*nstates, [0.0]*nstates, [0.0]*nstates, [0.0]*nstates_adi] 

    den = 1.0/float(num_sh_traj)

    for rep in xrange(0,4):

        for i in xrange(num_sh_traj):
            #=========== SE ===================
            for st in xrange(len(ave_Pop_se[rep])):
                ave_Pop_se[rep][st] = ave_Pop_se[rep][st] + den * Pop_se[i][rep][st]

            #=========== SH ===================
            ave_Pop_sh[rep][ istate[i][rep] ] = ave_Pop_sh[rep][ istate[i][rep] ] + den


    return ave_Pop_sh, ave_Pop_se



def init_coeffs(system, ist, bastyp, num_sh_traj, S_half, S_i_half, U_cur, S_mo_fmo, rnd, active_orb, active_orb_adi):
### NOTE: This function is not generic - only specific to considered systems and processes

# system (string: either "Pc-C60" or "Pc-C60-deca") - selection the molecular system
# ist (integer)- initial state index in selected basis
# bastyp (string, options: "D-FMO", "L", "A-FMO", "A-MO")
# num_sh_traj (integer): the number of SH trajectories
# S_half = S^{1/2}, where S_ij = <D-FMO_i|D-FMO_j>
# S_i_half = S^{-1/2}, where S_ij = <D-FMO_i|D-FMO_j>
# U_cur = matrix that diagonalizes: <D_FMO_i | H_el | D_FMO_j >
# U_cur_adi = matrix that diagonalizes <AO_i | H_el | AO_j >
# S_mo_fmo = overlap matrix between MO and FMOs:  <MO_i | FMO_j> = C_MO^T() * S_ao * { C_FMO } 
# S_mo_fmo_i - pseudoinverse of S_mo_fmo

    print "In init_coeffs"

    # Total number of orbitals
    nfrags = len(active_orb)
    nstates = 0 
    for fr in xrange(nfrags):
        nstates = nstates + len(active_orb[fr])

    nstates_adi = len(active_orb_adi[0])     


    if 1:   
        print "nfrags = ", nfrags
        print "nstates = ", nstates
        print "nstates_adi = ", nstates_adi
        print "S_half = \n"; S_half.show_matrix()
        print "S_i_half = \n"; S_i_half.show_matrix()
        print "U_cur = \n"; U_cur.show_matrix()
        if S_mo_fmo!=None:
            print "S_mo_fmo = \n"; S_mo_fmo.show_matrix()
     


    # Setup the coefficients of the coherent superposition in all representations
    #  |PSI> = sum_i { c_i  * |i> }
    Coeff = []   
    for i in xrange(num_sh_traj):

        coeff_i = [CMATRIX(nstates, 1), CMATRIX(nstates, 1), CMATRIX(nstates, 1), CMATRIX(nstates_adi, 1)]

        # Setup coefficients in the selected basis
        if bastyp in [0,1,2]: 
            if ist>nstates-1:
                print "Error in init_coeffs: ist = ", ist, " max allowed state index = ", nstates-1
                sys.exit(0)
            coeff_i[bastyp].set(ist, 1.0, 0.0)                # D-FMO, L, or A-FMO

        elif bastyp==3:
            if ist>nstates_adi-1:
                print "Error in init_coeffs: ist = ", ist, " max allowed state index = ", nstates_adi-1
                sys.exit(0)
            coeff_i[3].set(ist, 1.0, 0.0)                # A-MO

        # Transform to all other representations
        rep_transform(bastyp, S_half, S_i_half, U_cur, S_mo_fmo, coeff_i)

        if 0:
            print i
            print "coeff[0] = "; coeff_i[0].show_matrix()
            print "coeff[1] = "; coeff_i[1].show_matrix()
            print "coeff[2] = "; coeff_i[2].show_matrix()
            print "coeff[3] = "; coeff_i[3].show_matrix()

        Coeff.append(coeff_i)

    # Setup the populations in all representations
    # P_i = <X|PSI><PSI|X>.get(i,i)  <-- diagonal element of the density matrix
    Pop_se = []
    istate = []

#    sys.exit(0)

    for i in xrange(num_sh_traj):
        # TS-SE
        pop_i = populations(bastyp, S_half, S_i_half, U_cur, S_mo_fmo, Coeff[i])
        Pop_se.append(pop_i)

        # And initialize states
        istate_i = [0, 0, 0, 0]

        for rep in [0,1,2,3]:
            ksi = rnd.uniform(0.0, 1.0)
            istate_i[rep] = set_random_state(pop_i[rep], ksi)

        istate.append(istate_i)


    ave_Pop_sh = [[0.0]*nstates, [0.0]*nstates, [0.0]*nstates, [0.0]*nstates_adi] 
    ave_Pop_se = [[0.0]*nstates, [0.0]*nstates, [0.0]*nstates, [0.0]*nstates_adi] 

    den = 1.0/float(num_sh_traj)

    for rep in xrange(0,4):

        for i in xrange(num_sh_traj):
            #=========== SE ===================
            for st in xrange(len(ave_Pop_se[rep])):
                ave_Pop_se[rep][st] = ave_Pop_se[rep][st] + den * Pop_se[i][rep][st]

            #=========== SH ===================
            ave_Pop_sh[rep][ istate[i][rep] ] = ave_Pop_sh[rep][ istate[i][rep] ] + den



    return nstates, nstates_adi, Coeff, Pop_se, istate, ave_Pop_se, ave_Pop_sh




        

def main():

    opt_ham = 2      # 0 - neglect electronic couplings, 1 - MO-EHT formula, 2 - transformation formula
    do_collapse = 0  # 0 - no decoherence, 1 - decoherence
    compute_adi = 1  # 0 - do not compute (no related info will be available), 1 - do compute
                     # this is needed for computing projections (even if the fragmentation is sed)
    sh_method = 0    # 0 - MSSH,  1 - FSSH


#    print_initial_mo = False
#    print_starting_mo = False
#    print_basis_mo = False

    print_initial_mo = True   # flag controlling whether we want to compute and print MOs (cube files) in the starting configuration
    print_starting_mo = True  # -- for the very first step of NA-MD (so after cooling and thermalization)
    print_basis_mo = True     # -- for the actual basis states that will be used (Lowdin orbitals, adiabatic-FMO or full MO)

    num_sh_traj = 1000;

#    system = "2benz";  print_amt = 1
    system = "Pc-C60";  print_amt = 1;  N_C60 = 1
#    system = "Pc-C60-3";  print_amt = 0;   N_C60 = 3
#    system = "Pc-C60-5";  print_amt = 0;   N_C60 = 5
#    system = "Pc-C60-7";  print_amt = 0;   N_C60 = 7
#    system = "Pc-C60-9";  print_amt = 0;   N_C60 = 9


#    init_bastyp = 2    # initial basis type: 0 = D-FMO, 1 = L, 2 = A-FMO, and 3 = A-MO - the one in which we start with initial state
#    ist = 3*N_C60 + 1  # 3 states from each C60 + 2 states from SubPc -1 since indexing starts with 0
#    prop_bastyp = 2    # propagation basis type

    init_bastyp = 2
    ist =  4  #7    # if work in MO basis
    prop_bastyp = 3

    if prop_bastyp==0:
        print "prop_bastyp must be any of the orthogonal bases [1,2,3], not 0"    

        



    rnd = Random()
    #--------------------- Initialization ----------------------

    # Create Universe and populate it
    U = Universe(); LoadPT.Load_PT(U, "elements.dat")

    # Create force field
    uff = ForceField({"bond_functional":"Harmonic", "angle_functional":"Fourier",
                      "dihedral_functional":"General0", "oop_functional":"Fourier",
                      "mb_functional":"LJ_Coulomb","R_vdw_on":40.0,"R_vdw_off":55.0 })
    LoadUFF.Load_UFF(uff,"uff.dat")

    # Create molecular system and initialize the properties
    syst = System()

    if system=="Pc-C60":
        LoadMolecule.Load_Molecule(U, syst, "Pc-C60.ent", "pdb")
    elif system=="2benz":
        LoadMolecule.Load_Molecule(U, syst, "2benz_aa.ent", "pdb")
    else: #system=="Pc-C60-deca":
        LoadMolecule.Load_Molecule(U, syst, "Pc-C60_deca1.ent", "pdb")

    syst.determine_functional_groups(0)  # do not assign rings
    syst.init_fragments()
    print "Number of atoms in the system = ", syst.Number_of_atoms
    print "Number of bonds in the system = ", syst.Number_of_bonds
    print "Number of angles in the system = ", syst.Number_of_angles
    print "Number of dihedrals in the system = ", syst.Number_of_dihedrals
    print "Number of impropers in the system = ", syst.Number_of_impropers
    atlst1 = range(1,syst.Number_of_atoms+1)


    # Creating Hamiltonian and initialize it
    ham = Hamiltonian_Atomistic(1, 3*syst.Number_of_atoms)
    ham.set_Hamiltonian_type("MM")
    ham.set_interactions_for_atoms(syst, atlst1, atlst1, uff, 1, 0)  # 0 - verb, 0 - assign_rings
    ham.show_interactions_statistics()

    # Bind Hamiltonian and the system   
    ham.set_system(syst);   ham.compute();   print "Energy = ", ham.H(0,0), " a.u."

    # Electronic DOFs
    el = Electronic(1,0)

    # Nuclear DOFs
    mol = Nuclear(3*syst.Number_of_atoms)

    # Initialize MD variables
    nve_md.nve_md_init(syst, mol, el, ham)

    # Initialize EHT control options
    prms = Control_Parameters(); 

    #control_filename = "control_parameters_eht.dat"
    #if system=="Pc-C60" or system=="Pc-C60-deca" or system=="2benz":

    control_filename = "control_parameters_eht_pc_c60.dat"

    # Define fragments and fragment-specific properties. Use the atom IDs (start at 1)
    frag = []
    frag_adi = []

    if system=="Pc-C60":
        Pc = range(1,45);    print "Pc = ", Pc
        C60 = range(45,105);   print "C60 = ", C60    
        frag = [ Pc , C60 ]
        frag_adi = [ Pc + C60 ]


    elif system=="2benz":
        Bz1 = range(1,13);    print "Bz1 = ", Bz1
        Bz2 = range(13,25);   print "Bz2 = ", Bz2
 
        frag = [Bz1, Bz2]
        frag_adi = [ Bz1 + Bz2 ]


    elif system=="Pc-C60-3":
        Pc = range(1,45);    print "Pc = ", Pc
        C60_1 = range(45,105);   print "C60_1 = ", C60_1
        C60_2 = range(105,165);  print "C60_2 = ", C60_2
        C60_3 = range(165,225);  print "C60_3 = ", C60_3
        C60_4 = range(225,285);  print "C60_4 = ", C60_4
        C60_5 = range(285,345);  print "C60_5 = ", C60_5
        C60_6 = range(345,405);  print "C60_6 = ", C60_6
        C60_7 = range(405,465);  print "C60_7 = ", C60_7
        C60_8 = range(465,525);  print "C60_8 = ", C60_8
        C60_9 = range(525,585);  print "C60_9 = ", C60_9
        C60_10= range(585,645);  print "C60_10 = ", C60_10

        frag = [ Pc , C60_1 , C60_2 , C60_3 ] #, C60_4 , C60_5 , C60_7 , C60_9 ]
        frag_adi = [ Pc + C60_1 + C60_2 + C60_3 ] #+ C60_4 + C60_5 + C60_7 + C60_9 ]

    elif system=="Pc-C60-5":
        Pc = range(1,45);    print "Pc = ", Pc
        C60_1 = range(45,105);   print "C60_1 = ", C60_1
        C60_2 = range(105,165);  print "C60_2 = ", C60_2
        C60_3 = range(165,225);  print "C60_3 = ", C60_3
        C60_4 = range(225,285);  print "C60_4 = ", C60_4
        C60_5 = range(285,345);  print "C60_5 = ", C60_5
        C60_6 = range(345,405);  print "C60_6 = ", C60_6
        C60_7 = range(405,465);  print "C60_7 = ", C60_7
        C60_8 = range(465,525);  print "C60_8 = ", C60_8
        C60_9 = range(525,585);  print "C60_9 = ", C60_9
        C60_10= range(585,645);  print "C60_10 = ", C60_10

        frag = [ Pc , C60_1 , C60_2 , C60_3 , C60_4 , C60_5 ] #, C60_7 , C60_9 ]
        frag_adi = [ Pc + C60_1 + C60_2 + C60_3 + C60_4 + C60_5] # + C60_7 + C60_9 ]


    elif system=="Pc-C60-7":
        Pc = range(1,45);    print "Pc = ", Pc
        C60_1 = range(45,105);   print "C60_1 = ", C60_1
        C60_2 = range(105,165);  print "C60_2 = ", C60_2
        C60_3 = range(165,225);  print "C60_3 = ", C60_3
        C60_4 = range(225,285);  print "C60_4 = ", C60_4
        C60_5 = range(285,345);  print "C60_5 = ", C60_5
        C60_6 = range(345,405);  print "C60_6 = ", C60_6
        C60_7 = range(405,465);  print "C60_7 = ", C60_7
        C60_8 = range(465,525);  print "C60_8 = ", C60_8
        C60_9 = range(525,585);  print "C60_9 = ", C60_9
        C60_10= range(585,645);  print "C60_10 = ", C60_10

        frag = [ Pc , C60_1 , C60_2 , C60_3 , C60_4 , C60_5 , C60_6 , C60_7 ]
        frag_adi = [ Pc + C60_1 + C60_2 + C60_3 + C60_4 + C60_5 + C60_6 + C60_7 ]

    elif system=="Pc-C60-9":
        Pc = range(1,45);    print "Pc = ", Pc
        C60_1 = range(45,105);   print "C60_1 = ", C60_1
        C60_2 = range(105,165);  print "C60_2 = ", C60_2
        C60_3 = range(165,225);  print "C60_3 = ", C60_3
        C60_4 = range(225,285);  print "C60_4 = ", C60_4
        C60_5 = range(285,345);  print "C60_5 = ", C60_5
        C60_6 = range(345,405);  print "C60_6 = ", C60_6
        C60_7 = range(405,465);  print "C60_7 = ", C60_7
        C60_8 = range(465,525);  print "C60_8 = ", C60_8
        C60_9 = range(525,585);  print "C60_9 = ", C60_9
        C60_10= range(585,645);  print "C60_10 = ", C60_10

        frag = [ Pc , C60_1 , C60_2 , C60_3 , C60_4 , C60_5 , C60_6 , C60_7 , C60_8, C60_9 ] #, C60_10]
        frag_adi = [ Pc + C60_1 + C60_2 + C60_3 + C60_4 + C60_5 + C60_6 + C60_7 + C60_8 + C60_9 ] # + C60_10]



    # Convert to indices (start from 0)
    # Fragments
    nfrags = len(frag)
    for i in xrange(nfrags):
        for j in xrange(len(frag[i])):
            frag[i][j] = frag[i][j] - 1

    # Adiabatic sub-system
    for j in xrange(len(frag_adi[0])):
        frag_adi[0][j] = frag_adi[0][j] - 1



    #======= Generate subsystems =========
    sub_syst, sub_syst_adi = [], None

    # Fragments
    for fr in xrange(nfrags):
        # Atom types and coordinates of the atoms belonging to given fragment        
        grad = []; mol_at_types_fr = StringList(); R_fr = VECTORList()

        for i in frag[fr]:
            print fr, i
            mol_at_types_fr.append(syst.Atoms[i].Atom_element)
            R_fr.append(syst.Atoms[i].Atom_RB.rb_cm)
            grad.append(-1.0*syst.Atoms[i].Atom_RB.rb_force)

        # Generate sub-systems
        sub_syst.append( init_system.init_system(mol_at_types_fr, R_fr, grad, rnd, 0.0, 0.0, 0, os.getcwd()+"/elements.dat") )

    # Adiabatic sub-system
    grad = []; mol_at_types_fr = StringList(); R_fr = VECTORList()
    for i in frag_adi[0]:
        mol_at_types_fr.append(syst.Atoms[i].Atom_element)
        R_fr.append(syst.Atoms[i].Atom_RB.rb_cm)
        grad.append(-1.0*syst.Atoms[i].Atom_RB.rb_force)
    sub_syst_adi = [ init_system.init_system(mol_at_types_fr, R_fr, grad, rnd, 0.0, 0.0, 0, os.getcwd()+"/elements.dat")  ]



    #======= Now handle Hamiltonians =========
    E, homo, homo_adi = [0.0]*nfrags, [], 0
    sub_ham_old, sub_ham_cur = [], []  # fragment

 
    # Adiabatic sub-system.     
    ham_old = [ listHamiltonian_QM(control_filename, sub_syst_adi[0] ) ] # this includes the overlap calculation
    if compute_adi==0:
        ham_old[0].compute_core_Hamiltonian(sub_syst_adi[0]) 
    elif compute_adi==1:
        ham_old[0].compute_scf(sub_syst_adi[0])     # this includes the core Hamiltonian calculations
    ham_cur = [ listHamiltonian_QM(ham_old[0]) ]

    el_str_adi = ham_cur[0].get_electronic_structure()
    homo_adi = el_str_adi.Nocc_alp - 1  # index of the HOMO orbital in the adiabatic sub-system
    print "Adiabatic sub-system: number of electrons (alpha) = ", el_str_adi.Nocc_alp, " homo = ", homo_adi

    if compute_adi==1:
        print "Electronic structure of the adiabatic sub-system"
        print "Index  Bands(alp)    Occupations(alp)       Bands(bet)    Occupations(bet)"
        for j in xrange(el_str_adi.Norb):
            print "%5i  %12.8f   %12.8f  %12.8f   %12.8f" %(j, el_str_adi.get_bands_alp(j), el_str_adi.get_occ_alp(j), el_str_adi.get_bands_bet(j), el_str_adi.get_occ_bet(j) )
    

#    syst.exit(0)

    # Fragments
    for fr in xrange(nfrags):
        lstHamQM_fr = listHamiltonian_QM(control_filename, sub_syst[fr]) # this includes the overlap calculation
        E[fr] = lstHamQM_fr.compute_scf(sub_syst[fr])     # this includes the core Hamiltonian calculation
        sub_ham_old.append(lstHamQM_fr)
        sub_ham_cur.append(listHamiltonian_QM(lstHamQM_fr))

        el_str_fr = sub_ham_cur[fr].get_electronic_structure()
        homo_fr = el_str_fr.Nocc_alp - 1  # index of the HOMO orbital
        print "fragment ", fr, " number of electrons (alpha) = ", el_str_fr.Nocc_alp, " homo = ", homo_fr
        homo.append(homo_fr)

        print "Electronic structure of the fragment ",fr
        print "Index  Bands(alp)    Occupations(alp)       Bands(bet)    Occupations(bet)"
        for j in xrange(el_str_fr.Norb):
            print "%5i  %12.8f   %12.8f  %12.8f   %12.8f" %(j, el_str_fr.get_bands_alp(j), el_str_fr.get_occ_alp(j), el_str_fr.get_bands_bet(j), el_str_fr.get_occ_bet(j) )



    #======= Define spaces of active states (considered in NA-MD) =========
    active_orb = []
    active_orb_adi = []

    if system=="2benz":
        active_orb_adi = [ range(homo_adi+1,homo_adi+3) ]   # L, L+1
                                                            # 0   1  
        for fr in xrange(nfrags):
            active_orb.append([ homo[fr]+1 ] )  # L



    else: # system=="Pc-C60" or system=="Pc-C60-deca" or system=="2benz":
        active_orb_adi = [ range(homo_adi+1,homo_adi+10) ]   # L, L+1, L+2  L+3, L+4
                                                             # 0   1    2    3   
        for fr in xrange(nfrags):
            if fr==0:
                active_orb.append(range(homo[fr]+1, homo[fr]+3)) #   L, L+1    SubPc
                                                                 #   0   1
            else:
                active_orb.append(range(homo[fr]+1,homo[fr]+4))  #   L  L+1, L+2,   C60
                                                                 #   0   1    2    (of that fragment)


    #========= Print fragment orbitals here ===========
    if print_initial_mo:
        print_mo_frag(active_orb, sub_ham_cur, syst, "_initial_dia_")

        if compute_adi==1: 
            print_mo_adi(active_orb_adi[0], ham_cur[0], syst, "_initial_adia_")




    #=================== Propagation ====================
    ########################## Cooling #################################

    md = MD({"max_step":100,"ensemble":"NVE","integrator":"DLML","terec_exp_size":10,"dt":40.0,"n_medium":1,"n_fast":1,"n_outer":1}); 
    md.show_info()

    # Thermostat
    therm = Thermostat({"Temperature":278.0,"Q":100.0,"thermostat_type":"Nose-Hoover","nu_therm":0.001,"NHC_size":5});  
    therm.show_info()

    # State: system + thermostat + (optional barostat)
    ST = State()
    ST.set_system(syst); ST.set_thermostat(therm); ST.set_md(md)
    ST.init_md(mol, el, ham, rnd)    

#    anneal_schedule = [ {"dt":20.0, "nsteps":5, "ncycles":20}] # for test (very fast)
#    anneal_schedule = [ {"dt":20.0, "nsteps":50, "ncycles":20}] # for test
    anneal_schedule = [ {"dt":20.0, "nsteps":500, "ncycles":20}]
    f = open("_en_cooling.txt","w"); f.close()

    for item in anneal_schedule:
        md.dt = item["dt"]; md.max_step = item["nsteps"] 

        for i in xrange(item["ncycles"]):
            syst.set_atomic_q(mol.q)
            syst.print_xyz("_mol_cooling.xyz",i)
            ST.run_md(mol, el, ham)
            ekin = ST.E_kin; epot = ST.E_pot
            ST.cool()

            f = open("_en_cooling.txt","a")
            f.write("i= %3i ekin= %8.5f  epot= %8.5f  etot= %8.5f  H_NP= %8.5f  curr_T= %8.5f\n" % (i, ekin, epot, ST.E_tot, ST.H_NP, ST.curr_T ))
            f.close()

    ########################## Production MD - thermalization #################################

    syst.init_atom_velocities(300.0, rnd)  # must be this!!!
    # syst.init_fragment_velocities(300.0) not this!

    md.max_step = 10;  md.ensemble = "NVT"; md.dt = 20.0;
#    md.dt = 0.0

    f = open("_en_md-therm.txt","w"); f.close()

#    for i in xrange(5):  # for test (very fast)
#    for i in xrange(50):  # for test
    for i in xrange(500):
        syst.set_atomic_q(mol.q)
        syst.print_xyz("_mol_md-therm.xyz",i)
        ST.run_md(mol, el, ham)

        f = open("_en_md-therm.txt","a")
        f.write("i= %3i ekin= %8.5f  epot= %8.5f  etot= %8.5f  H_NP= %8.5f  curr_T= %8.5f\n" % (i, ST.E_kin, ST.E_pot, ST.E_tot, ST.H_NP, ST.curr_T ))
        f.close()

    ########################## Production MD - NA-MD #################################
    f = open("_en_md.txt","w"); f.close()
    f_ham = open("_Hvib.txt","w"); f_ham.close()

    md.max_step = 1;  md.ensemble = "NVE";  md.dt = 20.0;




    nstates, nstates_adi, Coeff, Pop_se, istate, ave_Pop_se, ave_Pop_sh = None, None, None, None, None, None, None

    for i in xrange(10000):
        syst.set_atomic_q(mol.q)
        syst.print_xyz("_mol_md.xyz",i)
        ST.run_md(mol, el, ham)



        #=============== Now handle the electronic structure ==========================

        t = Timer();  t.start()

        #== Adiabatic sub-system ==
        # Update coordinates of the sub-system atoms
        sz_fr = len(frag_adi[0])   # the # of atoms in adiabatic sub-system
        for i_fr in xrange(sz_fr): 
            sub_syst_adi[0].Atoms[i_fr].Atom_RB.rb_cm = syst.Atoms[frag_adi[0][i_fr]].Atom_RB.rb_cm
        
        # Update positions of the basis AOs!!!
        print ham_cur[0].Norb
        for o in xrange(ham_cur[0].Norb):
            at_indx_loc = ham_cur[0].ao_to_atom_map[o]  # indices of the atoms in the system bound to this Hamiltonin object, that
                                                        # is sub_syst_adi !!!  so we will need another mapping of the indices of those atoms to
                                                        # the indices to the overall system atoms
            at_indx_glob = frag_adi[0][at_indx_loc]
            ham_cur[0].basis_ao[o].set_position(syst.Atoms[at_indx_glob].Atom_RB.rb_cm)

        # Update overlaps and core Hamiltonian
        ham_cur[0].compute_overlap(sub_syst_adi[0])        
        if compute_adi==0:
            ham_cur[0].compute_core_Hamiltonian(sub_syst_adi[0]) 
        elif compute_adi==1:
            ham_cur[0].compute_scf(sub_syst_adi[0])     # this includes the core Hamiltonian calculations

        t.stop();  print "Time for electronic structure calculation in the adiabatic sub-system = ", t.show()


        t = Timer();  t.start()
        #== Fragments ==
        for fr in xrange(nfrags):
            # Update coordinates of the sub-system atoms
            sz_fr = len(frag[fr])
            for i_fr in xrange(sz_fr):                           
                sub_syst[fr].Atoms[i_fr].Atom_RB.rb_cm = syst.Atoms[frag[fr][i_fr]].Atom_RB.rb_cm

            # Update positions of the basis AOs!!!
            print sub_ham_cur[fr].Norb
            for o in xrange(sub_ham_cur[fr].Norb):
                at_indx_loc = sub_ham_cur[fr].ao_to_atom_map[o]
                at_indx_glob = frag[fr][at_indx_loc]
                sub_ham_cur[fr].basis_ao[o].set_position(syst.Atoms[at_indx_glob].Atom_RB.rb_cm)

            # Update overlap, MOs and energies
            sub_ham_cur[fr].compute_overlap(sub_syst[fr])
            E[fr] = sub_ham_cur[fr].compute_scf(sub_syst[fr])   # this includes the core Hamiltonian calculations
            
        t.stop();  print "Time for electronic structure calculation on all fragments = ", t.show()


        #======= Update vibronic Hamiltonians (different options) =========== 

        t = Timer(); t.start()
        # Compute all the options (at least for now)
        # Lowdin-orthogonalized basis
        print "Compute Hvib_L"
        Hvib_low, S_half, S_i_half = update_Hvib_L(opt_ham, sub_ham_old, sub_ham_cur, ham_old[0], ham_cur[0], active_orb, md.dt, print_amt)        
        # Adiabatic-FMO basis
        print "Compute Hvib_A"
        Hvib_fmo, U_cur = update_Hvib_A(opt_ham, sub_ham_old, sub_ham_cur, ham_old[0], ham_cur[0], active_orb, md.dt, print_amt)           

 
        # Adiabatic MO basis (so U_cur_adi will be an identity matrix)
        Hvib_adi, S_mo_fmo, U_cur_adi = None, None, None
        if compute_adi==1:
            print "Compute Hvib_A-MO"
            Hvib_adi, U_cur_adi = update_Hvib_A(opt_ham, ham_old, ham_cur, ham_old[0], ham_cur[0], active_orb_adi, md.dt, print_amt)            
            S_mo_fmo = mo_fmo_overlap(sub_ham_cur, active_orb, ham_cur[0], active_orb_adi[0])  # nmo_adi x nmo
            print "S_mo_fmo = "; S_mo_fmo.show_matrix()


        Hvib = None
        if prop_bastyp==1:
            Hvib = Hvib_low
        elif prop_bastyp==2:
            Hvib = Hvib_fmo
        elif prop_bastyp==3:
            Hvib = Hvib_adi
            
        t.stop();   print "Time for update_Hvib = ", t.show()


        if i==0:  # initialization
            nstates, nstates_adi, Coeff, Pop_se, istate, ave_Pop_se, ave_Pop_sh = init_coeffs(syst, ist, init_bastyp, num_sh_traj, S_half, S_i_half, U_cur, S_mo_fmo, rnd, active_orb, active_orb_adi)


        if i==0:
            print "nstates = ", nstates
            print "nstates_adi = ", nstates_adi
            print "Coeff = ", Coeff
            print "Pop_se = ", Pop_se
            print "istate = ", istate
            print "ave_Pop_se = ", ave_Pop_se
            print "ave_Pop_sh = ", ave_Pop_sh


        #=========== Here we will print all the basis orbitals (Lowdin, adiabatic-FMO, or full MO)
        if i==0 and print_starting_mo==True:
            print_mo_frag(active_orb, sub_ham_cur, syst, "_thermalized_dia_")
            
            if compute_adi==1: 
                print_mo_adi(active_orb_adi[0], ham_cur[0], syst, "_thermalized_adi_")


        if i==0 and print_basis_mo==True:            
            prms1 = Control_Parameters()
            prms1.nx_grid, prms1.ny_grid, prms1.nz_grid = 40, 40, 40

            # Lowdin
            prms1.orbs = Py2Cpp_int(range(0,nstates)) 
            prms1.charge_density_prefix = "_thermalized_lowdin_fmo_"
            charge_density(S_i_half.real(), sub_ham_cur, syst, active_orb, prms1);

            # Adiabatic FMO
            prms1.orbs = Py2Cpp_int(range(0,nstates)) 
            prms1.charge_density_prefix = "_thermalized_adi_fmo_"
            charge_density(U_cur.real(), sub_ham_cur, syst, active_orb, prms1);

            # Adiabatic MO
            if compute_adi==1:
                prms1.orbs = Py2Cpp_int(range(0,len(active_orb_adi[0]))) 
                prms1.charge_density_prefix = "_thermalized_adi_mo_"
                charge_density(U_cur_adi.real(), ham_cur, syst, active_orb_adi, prms1);



        f_ham = open("_Hvib.txt","a"); 
        line = "i= %3i " % i

        if prop_bastyp in [0,1,2]:
            start_f = 0
            for fr in xrange(nfrags):
                print "Active orbitals of the fragment ", fr
                sz_f = len(active_orb[fr])
         
                for act_o in xrange(sz_f):
                    ff = active_orb[fr][act_o]
                    e_el  = sub_ham_cur[fr].get_electronic_structure().get_E_alp().get(ff, ff)
                    e_vib = Hvib.get(start_f + act_o, start_f + act_o).real
                   
                    print act_o, ff, e_el, e_vib
                    line = line + " %8.5f  %8.5f " % (e_el, e_vib)
         
                start_f = start_f + sz_f
        line = line + "\n"
        f_ham.write(line)
        f_ham.close()


        f_ham = open("_Hvib_mo.txt","a"); 
        line = "i= %3i " % i

        if prop_bastyp in [0,1,2,3]:
            start_f = 0
            for fr in [0]:
                print "Active orbitals of the fragment ", fr
                sz_f = len(active_orb_adi[fr])
         
                for act_o in xrange(sz_f):
                    ff = active_orb_adi[fr][act_o]
                    e_el  = ham_cur[fr].get_electronic_structure().get_E_alp().get(ff, ff)
                    e_vib = Hvib_adi.get(start_f + act_o, start_f + act_o).real
                   
                    print act_o, ff, e_el, e_vib
                    line = line + " %8.5f  %8.5f " % (e_el, e_vib)
         
                start_f = start_f + sz_f

        line = line + "\n"
        f_ham.write(line)
        f_ham.close()



        #============== TD-SE solution and surface hopping =================
        for tr in xrange(num_sh_traj):
            # Coherent dynamics
            propagate_electronic(md.dt, Coeff[tr][prop_bastyp], Hvib)  # propagate only in the basis selected


            g = None
            if sh_method==0:
                #g = compute_hopping_probabilities_mssh(Coeff[tr][prop_bastyp])
                g = compute_hopping_probabilities_mssh(Coeff[tr][prop_bastyp],Hvib,1,therm.Temperature)
            elif sh_method==1:
                #g = compute_hopping_probabilities_fssh(Coeff[tr][prop_bastyp], Hvib, md.dt)
                g = compute_hopping_probabilities_fssh(Coeff[tr][prop_bastyp], Hvib, md.dt,1,therm.Temperature)
            #if tr==0:
            #    print "Hopping matrix for the first trajectory is: "; g.show_matrix()

            ksi = rnd.uniform(0.0, 1.0)
            ksi2 = rnd.uniform(0.0, 1.0)
            
            old_st = istate[tr][prop_bastyp]
            new_st = hop(istate[tr][prop_bastyp], g, ksi)


            if new_st != old_st:
                E_old = Hvib.get(old_st,old_st).real
                E_new = Hvib.get(new_st,new_st).real

                # ID-A decoherence                
                istate[tr][prop_bastyp],Coeff[tr][prop_bastyp] = ida_py(Coeff[tr][prop_bastyp], old_st, new_st, E_old, E_new, 300.0, ksi2, do_collapse)


            Pop_se[tr] = populations(prop_bastyp, S_half, S_i_half, U_cur, S_mo_fmo, Coeff[tr])


        ave_Pop_sh, ave_Pop_se = ave_Pop(nstates, nstates_adi, num_sh_traj, Pop_se, istate)




        #===================== Start printing out ==============================

        #========= Fragment states ================

        for rep in [0,1,2]:

            f = open("_populations_rep%i.txt" % rep,"a")
            line = "i= %3i " % i    
            tot_pop = 0.0
            for st in xrange(nstates):
                pop_o = ave_Pop_se[rep][st] 
                tot_pop = tot_pop + pop_o
                line = line + " %8.5f " % pop_o
            line = line + "%8.5f\n" % tot_pop
            f.write(line)
            f.close()


            f = open("_populations_sh_rep%i.txt" % rep,"a")
            line = "i= %3i " % i
            tot_pop = 0.0
            for st in xrange(nstates):
                pop_o = ave_Pop_sh[rep][st] 
                tot_pop = tot_pop + pop_o
                line = line + " %8.5f " % pop_o
            line = line + "%8.5f\n" % tot_pop
            f.write(line)
            f.close()


        #=============== Adiabatic MO states ========================
        rep = 3
        if 1:
            f = open("_populations_rep%i.txt" % rep,"a")
            line = "i= %3i " % i    
            tot_pop = 0.0
            for st in xrange(nstates_adi):
                pop_o = ave_Pop_se[rep][st]
                tot_pop = tot_pop + pop_o
                line = line + " %8.5f " % pop_o
            line = line + "%8.5f\n" % tot_pop
            f.write(line)
            f.close()


            f = open("_populations_sh_rep%i.txt" % rep,"a")
            line = "i= %3i " % i
            tot_pop = 0.0
            for st in xrange(nstates_adi):
                pop_o = ave_Pop_sh[rep][st]
                tot_pop = tot_pop + pop_o
                line = line + " %8.5f " % pop_o
            line = line + "%8.5f\n" % tot_pop
            f.write(line)
            f.close()



        for fr in xrange(nfrags):
            f = open("_fragment_"+str(fr)+"elec_struct.txt","a")
            line = "i= %3i E_elec= %8.5f " % (i, E[fr])
            el_str = sub_ham_cur[fr].get_electronic_structure()
            for o in xrange(el_str.Norb):
                line = line + "%8.5f  "  %  el_str.get_E_alp().get(o,o) #mo_alp[fr][3][o][1]
            line = line + "\n"
            f.write(line)
            f.close()

        f = open("_en_md.txt","a")
        f.write("i= %3i ekin= %8.5f  epot= %8.5f  etot= %8.5f  H_NP= %8.5f  curr_T= %8.5f\n" % (i, ST.E_kin, ST.E_pot, ST.E_tot, ST.H_NP, ST.curr_T ))
        f.close()


        #============== Update: new becomes old ==========================
        if opt_ham==2:
            #del ham_old[0]
            ham_old[0] = listHamiltonian_QM(ham_cur[0])
            

        for fr in xrange(nfrags):
            #del sub_ham_old[fr]
            sub_ham_old[fr] = listHamiltonian_QM(sub_ham_cur[fr])



main()

