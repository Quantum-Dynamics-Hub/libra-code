#***********************************************************
# * Copyright (C) 2017-2018 Wei Li and Alexey V. Akimov
# * This file is distributed under the terms of the
# * GNU General Public License as published by the
# * Free Software Foundation; either version 3 of the
# * License, or (at your option) any later version.
# * http://www.gnu.org/copyleft/gpl.txt
#***********************************************************/

import os
import sys

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *


def get_value(params,key,default,typ):
# Function to extract parameter from the dictionary
    # Try to get value from the dictionary
    str_val = "None"
    if key in params:
        if params[key]!=None:
            str_val = params[key]

    # If nothing found - use default value
    if str_val!="None":
        pass  # All is ok
    else: 
        str_val = default
        print "Warning: Parameter with key = %s does not exist in dictionary" % key
        print "Using the default value of %s" % default

    # Convert string to desired data type
    if typ=="s":
        return str_val
    elif typ=="f":
        return float(str_val)
    elif typ=="i":
        return int(float(str_val))




def split_orbitals_energies(C, E):
    """ In SOC, non-collinear case, the orbitals are 2-component spinors:
             | psi_i^alp |          | E_i_alp          |
     psi_i = |           |,  so E = |                  |
             | psi_i^bet |          |          E_i_bet |
  
    So, the wfc we read from the QE calculations, in this case is composed of such
    pairs, going in this order, so:

    psi_i^alp, psi_i^bet, psi_{i+1}^alp, psi_{i+1}^bet, ...

    Thus, there are 2*N_mo_adi columns, so we need to extract spin-components 
    into 2 matrices

    """

    N_pw = C.num_of_rows
    N_mo_adi = C.num_of_cols/2   # C.num_of_cols has to be even

    C_alp = CMATRIX(N_pw, N_mo_adi)
    C_bet = CMATRIX(N_pw, N_mo_adi)
    E_alp = CMATRIX(N_mo_adi, N_mo_adi)
    E_bet = CMATRIX(N_mo_adi, N_mo_adi)


    stenc1, stenc2 = [], []

    for i in xrange(N_mo_adi):
        stenc1.append(2*i)
        stenc2.append(2*i+1)

    pop_submatrix(C, C_alp, range(0, N_pw), stenc1 )
    pop_submatrix(C, C_bet, range(0, N_pw), stenc2 )

    pop_submatrix(E, E_alp, stenc1, stenc1)
    pop_submatrix(E, E_bet, stenc2, stenc2)

    return C_alp, C_bet, E_alp, E_bet

    



def merge_orbitals(Ca, Cb):
    """
    This function puts two matrices together into a single matrix
    """

    npw_a = Ca.num_of_rows
    N_mo_a = Ca.num_of_cols
 
    npw_b = Cb.num_of_rows
    N_mo_b = Cb.num_of_cols

    if npw_a != npw_b:
        print "Error: The number of rows of the two matrices should be equal"
        sys.exit(0)

    C = CMATRIX(npw_a, N_mo_a + N_mo_b)

    push_submatrix(C, Ca, range(0,npw_a), range(0,N_mo_a));
    push_submatrix(C, Cb, range(0,npw_a), range(N_mo_a, N_mo_a + N_mo_b));

    return C



def post_process(coeff, ene, issoc):

    """

    issoc = 0 - no SOC, = 1 - SOC

    In SOC case (spinor):

        coeff[0] - a N_pw x 2*N matrix of type:   (psi_0^alp, psi_0^bet, ... psi_{N-1}^alp, psi_{N-1}^bet) 
                   where each psi_ is a colum of the PW coefficients for given KS orbital (spatial component 
                   for each spin)

        ene[0] -   a 2*N x 2*N matrix of energies coming in pairs:  e_{2*i} = e_{2*i+1}, because these are the 
                   energies of the same orbital, just its different spin components

   
        We then split the coeff[0] matrix into a pair of N_pw x N matrices that represent alpha and beta
        spatial components of the wavefunction, separately

        However, the number of <b>spin<\b>-orbitals will be twice that of the N,
        so we need to construct new matrices with spin-orbitals.

        So that we have both psi_i = (psi_i_alp, psi_i_bet) and psi_{i+N} = (psi_i_bet, psi_i_alp)
        pairs of spin-orbitals (this is needed to represent the indistinguishable nature of electrons)
        i = 0,...N-1, where N - is the number of pairs of read spinors  = N_adi_ks_orb


    In non-SOC case (spin-polarized):
    We directly get a pair of N_pw x N_mo_dia matrices that represent alpha and beta spatial components of the wavefunction

        coeff[0] - a N_pw x N matrix of type:   (psi_0^alp, ... psi_{N-1}^alp) 
        coeff[1] - a N_pw x N matrix of type:   (psi_0^bet, ... psi_{N-1}^bet) 

                   where each psi_ is a colum of the PW coefficients for given KS orbital (spatial component 
                   for each spin)
        ene[0] -   a N x N matrix of alpha KS orbital energies
        ene[1] -   a N x N matrix of beta KS orbital energies



    Eventually:

                "alpha-block"        "beta-block"

    C_adi[0] = (psi_0^alp,... psi_{N-1}^alp,  psi_N^alp, ... psi_{2N-1}^alp)     alpha-components of spinors 
    C_adi[1] = (psi_0^bet,... psi_{N-1}^bet,  psi_N^bet, ... psi_{2N-1}^bet)     beta-components of spinors

                  ^                              ^
                  |______________________________|
                                |
        same energies in SOC case or non-polarized non-SOC
        different energies in spin-polarized non-SOC
        

    Also:  psi_0^alp = psi_N^bet  and psi_0^bet = psi_N^alp, and so on

             | E^alpha-block          0         |
    E_adi =  |                                  |
             |      0             E^beta-block  |

    E^alpha-block  = E^beta-block (SOC or non-polarized non-SOC)
    E^alpha-block != E^beta-block (spin-polarized non-SOC)

    """

    C, E = None, None

    if issoc==1:  # SOC, spinor case

        c_a, c_b, e_a, e_b = split_orbitals_energies(coeff[0], ene[0])

        N_ks_orb = c_a.num_of_cols  # should be equal to Cb.num_of_cols

        C_a = merge_orbitals(c_a, c_b)
        C_b = merge_orbitals(c_b, c_a)
        C = [C_a, C_b]  # "2-component spinor" format

        # Same with energies:
        E = CMATRIX(2*N_ks_orb, 2*N_ks_orb)
        push_submatrix(E, e_a, range(0,N_ks_orb), range(0,N_ks_orb) )
        push_submatrix(E, e_b, range(N_ks_orb,2*N_ks_orb), range(N_ks_orb,2*N_ks_orb) )


    elif issoc==0:  # no SOC

        if len(coeff)==1:  # spin-unpolarized (1-k point)

            N_ks_orb = coeff[0].num_of_cols  
        
            zero = CMATRIX(coeff[0].num_of_rows, N_ks_orb)

            C_a = merge_orbitals(coeff[0], zero)
            C_b = merge_orbitals(zero, coeff[0])
            C = [C_a, C_b]  # "2-component spinor" format

            # Same with energies:
            E = CMATRIX(2*N_ks_orb, 2*N_ks_orb)
            push_submatrix(E, ene[0], range(0,N_ks_orb), range(0,N_ks_orb) )
            push_submatrix(E, ene[0], range(N_ks_orb,2*N_ks_orb), range(N_ks_orb,2*N_ks_orb) )


        elif len(coeff)==2:  # spin-polarized (or 2-k points, non-polarized case, beware!)

            N_ks_orb = coeff[0].num_of_cols  
        
            zero = CMATRIX(coeff[0].num_of_rows, N_ks_orb)

            C_a = merge_orbitals(coeff[0], zero)
            C_b = merge_orbitals(zero, coeff[1])
            C = [C_a, C_b]  # "2-component spinor" format

            # Same with energies:
            E = CMATRIX(2*N_ks_orb, 2*N_ks_orb)
            push_submatrix(E, ene[0], range(0,N_ks_orb), range(0,N_ks_orb) )
            push_submatrix(E, ene[1], range(N_ks_orb,2*N_ks_orb), range(N_ks_orb,2*N_ks_orb) )


    return C, E




def orthogonalize_orbitals(C):
    """
    C = N_pw x N_mo   (just alpha or beta orbitals)
    flag == 2:   C = N_pw x 2*N_mo (both alpha and beta orbitals)

    This function takes an input of orbitals (C), which may not
    be rigorously orthogonal, finds a suitable transformation (U)
    and converts them into rigorously orthogonal orbitals (C_tilda)

    C_tilda = C * U, so if you want

    C_tilda^+  * C_tilda = I, you'll find that

    U = S^{-1/2}, where S = C^+ * C

    """

    S = C.H() * C  # overlap matrix

    # Test is S is invertabile
    print "\nTesting if S is invertabile\n"
    print FullPivLU_rank_invertible(S)
    print "Det = ", FullPivLU_det(S)
    is_inv = FullPivLU_rank_invertible(S)
    if is_inv[1] != 1:
        print "Error, S is not invertible, Exiting Program" 
        sys.exit(0)

    S_half = CMATRIX(S.num_of_rows, S.num_of_cols)
    S_i_half = CMATRIX(S.num_of_rows, S.num_of_cols)

    sqrt_matrix(S, S_half, S_i_half)

    C_tilda = C * S_i_half

    return C_tilda


def orthogonalize_orbitals2(Ca,Cb):
    """
    Ca and Cb = N_pw x N_mo   - represent the spin-components
    of the adiabatic states

    This function takes an input of orbitals (C), which may not
    be rigorously orthogonal, finds a suitable transformation (U)
    and converts them into rigorously orthogonal orbitals (C_tilda)

    For each channel:

    C_tilda = C * U, so if you want

    C_tilda^+  * C_tilda = I, you'll find that

    U = S^{-1/2}, where S = Ca^+ * Ca + Cb^+ * Cb

    """

    S = Ca.H() * Ca + Cb.H() * Cb  # overlap matrix

    S_half = CMATRIX(S.num_of_rows, S.num_of_cols)
    S_i_half = CMATRIX(S.num_of_rows, S.num_of_cols)

    sqrt_matrix(S, S_half, S_i_half)

    Ca_tilda = Ca * S_i_half
    Cb_tilda = Cb * S_i_half

    return Ca_tilda, Cb_tilda


