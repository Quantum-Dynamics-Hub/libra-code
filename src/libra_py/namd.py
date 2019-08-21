#*********************************************************************************                     
#* Copyright (C) 2018-2019 Alexey V. Akimov                                                   
#*                                                                                                     
#* This file is distributed under the terms of the GNU General Public License                          
#* as published by the Free Software Foundation, either version 2 of                                   
#* the License, or (at your option) any later version.                                                 
#* See the file LICENSE in the root directory of this distribution   
#* or <http://www.gnu.org/licenses/>.          
#***********************************************************************************
"""
.. module:: namd
   :platform: Unix, Windows
   :synopsis: This module implements functions for running NA-MD
.. moduleauthor:: Alexey V. Akimov

"""

import os
import sys
import math
import copy

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *



def compute_Hvib(ham_old, ham_cur, orb, dt):
    """

    Computes the vibronic Hamiltonian - in the MO basis - this is equivalent to the general version 
    below with the 1-electron "Slater determinant" basis functions.
    This function is going to be deprecated, but let keep it here, for backward compatibility

    Args:     
        ham_old ( Hamiltonian ): Hamiltonian at time t-dt [units: a.u. of energy]
        ham_cur ( Hamiltonian ): Hamiltonian at time t [units: a.u. of energy]
        orb ( list of ints ): indices of the orbitals included in the active space. The Hvib dimensions will 
            be determined by the N_act = len(orb) and the elements of Hvib will reflect only the orbitals included 
            in this active state. Indexing starts with 0.
        dt ( float ): Time step [units: a.u. of time]

    Returns:
        CMATRIX(N_act,N_act): vibronic Hamiltonian matrix in the MO basis:  Hvib = Hel - i*hbar*d_ij

    """

    N_act = len(orb)  # number of orbitals in the active space

    es_old = ham_old.get_electronic_structure()
    Nocc = es_old.Nocc_alp
    Fao_old = CMATRIX(es_old.get_Fao_alp() )
    Sao_old = CMATRIX(es_old.get_Sao() )
    res_old = Fock_to_P(Fao_old, Sao_old, Nocc, 1.0, 0.001, 1e-8, 0) 
    E_old = res_old[0] # MO energies
    C_old = res_old[1] # MO-LCAO
    e_old = CMATRIX(N_act, N_act)
    pop_submatrix(E_old, e_old, orb, orb)


    e_old.show_matrix()

    es_cur = ham_cur.get_electronic_structure()
    Nocc = es_cur.Nocc_alp
    Fao_cur = CMATRIX(es_cur.get_Fao_alp() )
    Sao_cur = CMATRIX(es_cur.get_Sao() )
    res_cur = Fock_to_P(Fao_cur, Sao_cur, Nocc, 1.0, 0.001, 1e-8, 0) 
    E_cur = res_cur[0] # MO energies
    C_cur = res_cur[1] # MO-LCAO
    e_cur = CMATRIX(N_act, N_act)
    pop_submatrix(E_cur, e_cur, orb, orb)


    e_cur.show_matrix()    


    # <MO(t)|MO(t-dt)>
    dist2 = 10000.0

    D = CMATRIX(N_act,N_act)
    MO_overlap(D, ham_old.basis_ao, ham_cur.basis_ao, C_old, C_cur, Py2Cpp_int(orb), Py2Cpp_int(orb), dist2 )  # 
    D = (0.5/dt)*(D - D.H())


    Eadi = 0.5*(e_old + e_cur)        
    Hvib = Eadi - 1.0j*D       # direct adiabatic

    
    return Hvib



def compute_Hvib_sd(ham_old, ham_cur, orb, SD_basis, dt):
    """ Computes the vibronic Hamiltonian - in the MO basis

    Args:     
        ham_old ( Hamiltonian ): Hamiltonian at time t-dt [units: a.u. of energy]
        ham_cur ( Hamiltonian ): Hamiltonian at time t [units: a.u. of energy]
        orb ( list of ints ): indices of the orbitals included in the active space. The Hvib dimensions will 
            be determined by the norb = len(orb) and the elements of Hvib will reflect only the orbitals included 
            in this active state. Indexing starts with 0.
        SD_basis ( list of lists of ints ): Slater determinant (SD) basis - are defined in terms of the indices 
            withing the active space, with the numbering starting from 0. In no-SOC case, the numbers smaller
            than N_act = len(orb) will correspond to alpha orbitals other - to beta orbitals
        dt ( float ): Time step [units: a.u. of time]

    Returns:
        CMATRIX(N,N): vibronic Hamiltonian matrix in the SD basis:  Hvib = Hel - i*hbar*d_ij. Here, N = len(SD_basis)
        
    
    Example:
        If we have a system with HOMO being the orbital with index 20, LUMO - with 21, then 
        if we want to consider a minimal 2-state system, we can choose orb = [20, 21]
        Mapping::
            0       1      2      3     
            H(alp) L(alp) H(bet)  L(bet) 
    
                     GS = |H(alp),H(beta)|    S1 = |H(alp),L(beta)|
     
            and SD_basis = [ [0,    2],             [0,      3]   ]
                             /\    /\                /\      /\
                             |     |                 |       |  
                          H alp   H bet             H alp   L beta

    Example:
        If we have a system with HOMO being the orbital with index 20, LUMO - with 21, then 
        if we want to consider a minimal 2-state system with spin-flip, we can choose orb = [20, 21]
        Mapping::

            0       1      2      3     
            H(alp) L(alp) H(bet)  L(bet) 
       
                    GS = |H(alp),H(beta)|     S1 =|H(alp),L(beta)|    T0 = |H(alp),L(alp)|    T1 = |H(bet),L(bet)|
                                                                                  
            and SD_basis = [ [0,    2],             [0,      3]  ,            [0,   1] ,               [ 2,    3 ] ]
                             /\    /\                /\      /\               /\    /\                  /\     /\
                             |     |                 |       |                |     |                   |      |   
                          H alp   H bet             H alp   L beta           H alp  L alp              H bet   L bet

    Example:
        If we have a radiacl in a system with 3 electrons, occuplying  HOMO-1 (doubly) and HOMO (singly).
        Say, the HOMO is the orbital with index 20, HOMO-1 - with 19, then 
        if we may want to consider the active space of 3 orbitals: [19, 20, 21]
        Mapping::
            0       1      2      3        4      5
            H-1(alp) H(alp) L(alp) H-1(bet) H(bet) L(bet)
    
                    GS = |H-1(alp),H-1(beta), H(alp)| 
                 
            and SD_basis = [  [0,       3,        1],    [0,       3,        4],       ...     ]
                              /\        /\       /\      /\        /\       /\  
                              |         |        |       |         |        |   
                             H-1 alp  H-1 bet   H alp   H-1 alp  H-1 bet   H bet

    Example:
        We may consider a transition between certain orbitals, say we are interested in transition between 
        LUMO+10 (index 30) and LUMO+15 (index 35), then we can choose orb = [30, 35]
        Note, the definition of the SD basis used in Example #1 can be re-used. 
        Mapping::
            0       1            2            3     
            L+10(alp) L+15(alp)  L+10(bet)   L+15(bet) 
    
                    GS = |L+10(alp),L+10(beta)|    S1 = |L+10(alp),L+15(beta)|
    
            and SD_basis = [ [0,      2],             [0,      3]   ]
                             /\      /\                /\      /\
                             |       |                 |       |    
                          L+10 alp  L+15 bet         L+10 alp L+15 beta

    """


    N_act = len(orb)  # number of orbitals in the active space

    es_old = ham_old.get_electronic_structure()
    Nocc = es_old.Nocc_alp
    Fao_old = CMATRIX(es_old.get_Fao_alp() )
    Sao_old = CMATRIX(es_old.get_Sao() )
    res_old = Fock_to_P(Fao_old, Sao_old, Nocc, 1.0, 0.001, 1e-8, 0) 
    E_old = res_old[0] # MO energies
    C_old = res_old[1] # MO-LCAO
    e_old = CMATRIX(N_act, N_act)
    pop_submatrix(E_old, e_old, orb, orb)


    es_cur = ham_cur.get_electronic_structure()
    Nocc = es_cur.Nocc_alp
    Fao_cur = CMATRIX(es_cur.get_Fao_alp() )
    Sao_cur = CMATRIX(es_cur.get_Sao() )
    res_cur = Fock_to_P(Fao_cur, Sao_cur, Nocc, 1.0, 0.001, 1e-8, 0) 
    E_cur = res_cur[0] # MO energies
    C_cur = res_cur[1] # MO-LCAO
    e_cur = CMATRIX(N_act, N_act)
    pop_submatrix(E_cur, e_cur, orb, orb)



    # <MO(t)|MO(t-dt)>
    dist2 = 10000.0

    # construct the overlaps in the MO (spatial orbitals) space
    Smo = CMATRIX(N_act,N_act)     
    MO_overlap(Smo, ham_old.basis_ao, ham_cur.basis_ao, C_old, C_cur, Py2Cpp_int(orb), Py2Cpp_int(orb) , dist2 )  # 


    # construct the overlaps in the MO (spin-orbitals) space    
    #
    #  No SOC:
    #
    # the orb x orb block is for alpha orbitals and orb + N_act x orb + N_act  - for beta orbitals
    Sso = CMATRIX(2*N_act,2*N_act) 

    orb_alp = Py2Cpp_int(range(0, N_act))
    orb_bet = Py2Cpp_int(range(N_act, 2*N_act))

    push_submatrix(Sso, Smo, orb_alp, orb_alp)
    push_submatrix(Sso, Smo, orb_bet, orb_bet)


    # construct the overlaps in the SD space
    # convert to internal data type:
    sd_basis = intList2()
    for sd in SD_basis:
        sd_basis.append( Py2Cpp_int(sd) )

    Ssd = overlap_sd(Sso, sd_basis);  
    Ssd = (0.5/dt)*(Ssd - Ssd.H())   # scalar NAC


    # Compute energies in the SD basis:
    sd_bas_size = len(SD_basis)
    Esd = CMATRIX(sd_bas_size, sd_bas_size)
    for I in range(0,sd_bas_size):
        e = 0.0
        
        for i in SD_basis[I]:
            if i<N_act:
                e = e + 0.5*(e_cur.get(i,i) + e_old.get(i,i))
            else:
                j = i - N_act
                e = e + 0.5*(e_cur.get(j,j) + e_old.get(j,j))

        Esd.set(I,I, e)


    # Finally, compute the Hvib
    Hvib = Esd - 1.0j*Ssd    

    
    return Hvib



