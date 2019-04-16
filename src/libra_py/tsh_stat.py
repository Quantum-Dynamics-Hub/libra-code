#*********************************************************************************                     
#* Copyright (C) 2019 Alexey V. Akimov                                                   
#* Copyright (C) 2016-2019 Kosuke Sato, Alexey V. Akimov                                                   
#*                                                                                                     
#* This file is distributed under the terms of the GNU General Public License                          
#* as published by the Free Software Foundation, either version 2 of                                   
#* the License, or (at your option) any later version.                                                 
#* See the file LICENSE in the root directory of this distribution   
#* or <http://www.gnu.org/licenses/>.          
#***********************************************************************************
"""
.. module:: tsh_stat
   :platform: Unix, Windows
   :synopsis: This module implements various functions for analysis of the statistics 
       in the TSH calculations
.. moduleauthor:: Alexey V. Akimov

"""

__author__ = "Alexey V. Akimov, Kosuke Sato"
__copyright__ = "Copyright 2016-2019 Kosuke Sato, Alexey V. Akimov"
__credits__ = ["Alexey V. Akimov", "Kosuke Sato"]
__license__ = "GNU-3"
__version__ = "1.0"
__maintainer__ = "Alexey V. Akimov"
__email__ = "alexvakimov@gmail.com"
__url__ = "https://quantum-dynamics-hub.github.io/libra/index.html"



import os
import sys
import math
import copy

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

import units
import probabilities


def compute_etot(ham, p, Cdia, Cadi, iM, rep):
    """Computes the Ehrenfest potential energy

    This function computes the average kinetic, potential, and total
    energies for an ensemble of trajectories - according to the Ehrenfest recipe

    Args: 
        ham ( nHamiltonian ): object that handles Hamiltonian-related calculations with many trajectories
        p ( MATRIX(ndof, ntraj) ): nuclear momenta of multiple trajectories
        Cdia ( CMATRIX(ndia, ntraj) ): amplitudes of diabatic states in the TD wavefunction expansion
        Cadi ( CMATRIX(ndia, ntraj) ): amplitudes of adiabatic states in the TD wavefunction expansion
        iM ( MATRIX(ndof, 1) ): inverse masses for all nuclear DOFs
        rep ( int ): The selector of the representation that is of current interest.

            - 0: diabatic
            - 1: adiabatic

    Returns: 
        tuple: ( Ekin, Epot, Etot, dEkin, dEpot, dEtot ): here

            * Ekin ( double ): average kinetic energy of the ensemble
            * Epot ( double ): average potential energy of the ensemble
            * Etot ( double ): average total energy of the ensemble
            * dEkin ( double ): standard deviation of the kinetic energy in the ensemble
            * dEpot ( double ): standard deviation of the potential energy in the ensemble
            * dEtot ( double ): standard deviation of the total energy in the ensemble

    """

    ntraj = p.num_of_cols
    ndof = p.num_of_rows

    epot, ekin = [], []    
    Epot, Ekin = 0.0, 0.0

    nst = 1
    if rep==0:
        nst = Cdia.num_of_rows
    elif rep==1:
        nst = Cadi.num_of_rows


    C = CMATRIX(nst, 1)

    for traj in xrange(ntraj):

        if rep==0:
            pop_submatrix(Cdia, C, Py2Cpp_int(range(0,nst)), Py2Cpp_int([traj]))    
            epot.append( ham.Ehrenfest_energy_dia(C, Py2Cpp_int([0,traj])).real )
            Epot = Epot + epot[traj]
        elif rep==1:
            pop_submatrix(Cadi, C, Py2Cpp_int(range(0,nst)), Py2Cpp_int([traj]))    
            epot.append( ham.Ehrenfest_energy_adi(C, Py2Cpp_int([0,traj])).real )
            Epot = Epot + epot[traj]

        tmp = 0.0
        for dof in xrange(ndof):
            tmp = tmp + 0.5 * iM.get(dof, 0) * (p.get(dof, traj) ** 2)
        ekin.append(tmp)
        Ekin = Ekin + ekin[traj]

    Ekin = Ekin / float(ntraj)
    Epot = Epot / float(ntraj)
    Etot = Ekin + Epot

    # Variances:
    dEkin, dEpot = 0.0, 0.0
    for traj in xrange(ntraj):
        dEkin = dEkin + (ekin[traj] - Ekin)**2
        dEpot = dEpot + (epot[traj] - Epot)**2

    dEtot = dEkin + dEpot

    dEkin = math.sqrt(dEkin/ float(ntraj))
    dEpot = math.sqrt(dEpot/ float(ntraj))
    dEtot = math.sqrt(dEtot/ float(ntraj))
    

    return Ekin, Epot, Etot, dEkin, dEpot, dEtot




def compute_etot_tsh(ham, p, Cdia, Cadi, act_states, iM, rep):
    """Compute the adiabatic potential energy

    This function computes the average kinetic, potential, and total
    energies for an ensemble of trajectories - according to the TSH recipe

    Args: 
        ham ( nHamiltonian ): object that handles Hamiltonian-related calculations with many trajectories
        p ( MATRIX(ndof, ntraj) ): nuclear momenta of multiple trajectories
        Cdia ( CMATRIX(ndia, ntraj) ): amplitudes of diabatic states in the TD wavefunction expansion
        Cadi ( CMATRIX(ndia, ntraj) ): amplitudes of adiabatic states in the TD wavefunction expansion
        iM ( MATRIX(ndof, 1) ): inverse masses for all nuclear DOFs
        rep ( int ): The selector of the representation that is of current interest.

            - 0: diabatic
            - 1: adiabatic

    Returns: 
        tuple: ( Ekin, Epot, Etot, dEkin, dEpot, dEtot ): here

            * Ekin ( double ): average kinetic energy of the ensemble
            * Epot ( double ): average potential energy of the ensemble
            * Etot ( double ): average total energy of the ensemble
            * dEkin ( double ): standard deviation of the kinetic energy in the ensemble
            * dEpot ( double ): standard deviation of the potential energy in the ensemble
            * dEtot ( double ): standard deviation of the total energy in the ensemble

    """

    ntraj = p.num_of_cols
    ndof = p.num_of_rows

    epot, ekin = [], []    
    Epot, Ekin = 0.0, 0.0

    nst = 1
    if rep==0:
        nst = Cdia.num_of_rows
    elif rep==1:
        nst = Cadi.num_of_rows


    C = CMATRIX(nst, 1)
    states = CMATRIX(nst, ntraj)

    tsh_indx2vec(ham, states, act_states)

    for traj in xrange(ntraj):

        pop_submatrix(states, C, Py2Cpp_int(range(0,nst)), Py2Cpp_int([traj]))      

        if rep==0:
            epot.append( ham.Ehrenfest_energy_dia(C, Py2Cpp_int([0,traj])).real )
            Epot = Epot + epot[traj]
        elif rep==1:
            epot.append( ham.Ehrenfest_energy_adi(C, Py2Cpp_int([0,traj])).real )
            Epot = Epot + epot[traj]

        tmp = 0.0
        for dof in xrange(ndof):
            tmp = tmp + 0.5 * iM.get(dof, 0) * (p.get(dof, traj) ** 2)
        ekin.append(tmp)
        Ekin = Ekin + ekin[traj]

    Ekin = Ekin / float(ntraj)
    Epot = Epot / float(ntraj)
    Etot = Ekin + Epot

    # Variances:
    dEkin, dEpot, dEtot = 0.0, 0.0, 0.0
    for traj in xrange(ntraj):
        dEkin = dEkin + (ekin[traj] - Ekin)**2
        dEpot = dEpot + (epot[traj] - Epot)**2
        dEtot = dEtot + (ekin[traj] + epot[traj] - Etot)**2

    dEkin = math.sqrt(dEkin/ float(ntraj))
    dEpot = math.sqrt(dEpot/ float(ntraj))
    dEtot = math.sqrt(dEtot/ float(ntraj))
    

    return Ekin, Epot, Etot, dEkin, dEpot, dEtot




def compute_dm(ham, Cdia, Cadi, rep, lvl):
    """

    Compute the trajectory-averaged density matrices in diabatic
    or adiabatic representations

    Args: 
        ham ( nHamiltonian ): object that handles Hamiltonian-related calculations with many trajectories
        Cdia ( CMATRIX(ndia, ntraj) ): amplitudes of diabatic states in the TD wavefunction expansion
        Cadi ( CMATRIX(ndia, ntraj) ): amplitudes of adiabatic states in the TD wavefunction expansion
        rep ( int ): a selector of which representation is considered main (being propagated)
            E.g. if rep = 0 - that means we propagate the diabatic coefficients, that is the calculation 
            of the diabatic density matrix is straightforward, but we need to involve some transformations 
            to compute the adiabatic density matrix and vice versa, if rep = 1, the propagation is done according
            to the adiabatic properties and we'd need to convert to the diabatic representation in the end
             
            - 0: diabatic
            - 1: adiabatic

        lvl ( int ): The level of the Hamiltonian that treats the transformations:
            - 0: ham is the actual Hamiltonian to use (use with single trajectory),
            - 1: ham is the parent of the Hamiltonians to use (use with multiple trajectories)

    Returns:
        tuple: ( dm_dia, dm_adi ):

            * dm_dia ( CMATRIX(ndia, ndia) ): the trajectory-averaged density matrix in
                the diabatic representation. Here, ndia - is the number of diabatic basis
                states
            * dm_adi ( CMATRIX(nadi, nadi) ): the trajectory-averaged density matrix in
                the adiabatic representation. Here, nadi - is the number of adiabatic basis
                states

    """

    ntraj = Cdia.num_of_cols
    ndia = Cdia.num_of_rows
    nadi = Cadi.num_of_rows

   
    dm_dia, dm_adi = CMATRIX(ndia, ndia), CMATRIX(nadi, nadi)


    for traj in xrange(ntraj):
        indx = None
        if lvl==0:
            indx = Py2Cpp_int([0])
        elif lvl==1:
            indx = Py2Cpp_int([0,traj])

    
        if rep==0:
            S = ham.get_ovlp_dia(indx)
            U = ham.get_basis_transform(indx) 
            #correct_phase(U)
    
            dm_tmp = S * Cdia.col(traj) * Cdia.col(traj).H() * S
            dm_dia = dm_dia + dm_tmp
            dm_adi = dm_adi + U.H() * dm_tmp * U
       
    
        elif rep==1:
            c = Cadi.col(traj)
            M = ham.get_ordering_adi(Py2Cpp_int([0, traj]))
            iM = inverse_permutation(M)

            c.permute_rows(iM)
            dm_tmp = c * c.H()
            dm_adi = dm_adi + dm_tmp

            S = ham.get_ovlp_dia(indx)
            U = ham.get_basis_transform(indx)     
            correct_phase(U)
            su = S * U
            dm_dia = dm_dia + su * dm_tmp * su.H()
    
    dm_dia = dm_dia / float(ntraj)        
    dm_adi = dm_adi / float(ntraj)

    return dm_dia, dm_adi



def compute_sh_statistics(nstates, istate):
    """

    This function computes the SH statistics for an ensemble of trajectories

    Args: 
        nstates ( int ): The number of considered quantum states
        istate ( list of integers ): The list containing the info about the index
            of a quantum state in which each trajectory is found. 
            The length of the list is equal to the number of trajectories. 
            Each element of the list is the state index for that trajectory. 
            In other words, istate[0] is the quantum state for a trajectory 0, 
            istate[1] is the quantum state for a trajectory 1, etc.

    Returns: 
        MATRIX(nstates, 1): coeff_sh: The list containing the average 
            population of each quantum state. The length of the list is equal to the 
            total number of quantum states considered, ```nstates```
    """

    ntraj = len(istate)
    f = 1.0/float(ntraj)

    coeff_sh = MATRIX(nstates, 1)

    for i in xrange(ntraj):
        st = istate[i]
        coeff_sh.add(st, 0, f)
 
    return coeff_sh


def update_sh_pop(istate, nstates):  
    """

    Args: 
        istate ( list of integers ): The list containing the info about the index
            of a quantum state in which each trajectory is found. 
            The length of the list is equal to the number of trajectories. 
            Each element of the list is the state index for that trajectory. 
            In other words, istate[0] is the quantum state for a trajectory 0, 
            istate[1] is the quantum state for a trajectory 1, etc.
        nstates ( int ): The number of considered quantum states

    Returns: 
        ( list of ```nstates``` ints ): pops: The list containing the average 
            SH-based population of each quantum state. The length of the list is equal to the 
            total number of quantum states considered, ```nstates```

    Note: 
        The functionality is the same as of ```compute_sh_statistics```, just a different 
        format of the output

    """

    pops = [0.0] * nstates
    ntraj = len(istate)

    incr = 1.0/float(ntraj)

    for j in xrange(ntraj): # for all trajectories
        pops[ istate[j] ] += incr

    return pops

    

def avarage_populations(el):
    """

    This function computes the SH statistics for an ensemble of trajectories

    Args:
        el ( list of Electronic ): The list containing electronic DOF variables
            for all trajectories in ensemble. The length of the list determines
            the number of trajectories in ensemble

    Returns:
        tuple: ( sh_pops, se_pops, rho ):

            * sh_pops ( list of N float ): The list containing the average population
                of each quantum state based on the statistics of the discrete states
                in which each trajectory resides. Here, N is the number of states
            * se_pops ( list of N float ): The list containing the average population
                of each quantum state based on the amplitudes of all quantum states
                as obtained from the TD-SE solution. Here, N is the number of states
            * rho ( CMATRIX(N,N) ): The matrix containing the trajectory-averaged 
                SE populations and coherences. Here, N is the number of states

    """

    ntraj = len(el)        # the total number of trajectories
    nstat = el[0].nstates  # the number of quantum states, assume that all objects in the "el" list 
                           # are similar

    sh_pops = [0.0] * nstat   # average SH populations of all states
    se_pops = [0.0] * nstat   # average SE populations of all states
    rho = CMATRIX(nstat, nstat) # trajectory-averaged density matrix

    f = 1.0/float(ntraj)
    for traj in xrange(ntraj): # for all trajectories
        sh_pops[ el[traj].istate ] += f

        for st1 in xrange(nstat):
            se_pops[ st1 ] += f * el[traj].rho(st1,st1).real

            for st2 in xrange(nstat):
                rho.set(st1, st2, rho.get(st1,st2) + f * el[traj].rho(st1,st2) )

    return sh_pops, se_pops, rho


def ave_pop(denmat_sh, denmat_se):
    """

    Compute the ensemble-averaged SH and SE density matrices

    Args: 
        denmat_sh ( list of CMATRIX(N, N) ): SH density matrix (diagonal)
            for each trajectory. Such that ```denmat_sh[itraj]``` corresponds to the trajectory ```itraj```
        denmat_se ( list of CMATRIX(N, N) ): SE density matrix
            for each trajectory. Such that ```denmat_se[itraj]``` corresponds to the trajectory ```itraj```

    Returns: 
        tuple: ( ave_pop_sh, ave_pop_se ):

            * ave_pop_sh ( CMATRIX(N, N) ): the ensemble-averaged SH density matrix
            * ave_pop_se ( CMATRIX(N, N) ): the ensemble-averaged SE density matrix

    """

    ntraj = len(denmat_sh)
    nst_out = denmat_sh[0].num_of_cols

    ave_pop_sh = CMATRIX(nst_out, nst_out)
    ave_pop_se = CMATRIX(nst_out, nst_out)
    den = 1.0/float(ntraj)

    for i in xrange(ntraj):
        ave_pop_se = ave_pop_se + den * denmat_se[i]   # SE
        ave_pop_sh = ave_pop_sh + den * denmat_sh[i]   # SH


    return ave_pop_sh, ave_pop_se



def ave_en(denmat_sh, denmat_se, Hvib):
    """Computes ensemble averaged SH and SE energies

    Args:
        denmat_sh ( list of CMATRIX(nst_in, nst_in) ): SE density matrices (diagonal) for each trajectory
        denmat_se ( list of CMATRIX(nst_in,nst_in) ): SE density matrices for each trajectory 
        Hvib ( list of CMATRIX(nst_in,nst_in) ): Hvib for each trajectory [units: arbitrary]

    Returns: 
        (double, double): ave_en_sh, ave_en_se, where:

            * ave_en_sh ( double ): SH-averaged energy [ units: same as Hvib ]
            * ave_en_se ( double ): SE-averaged energy [ units: same as Hvib ]
    """

    ntraj = len(denmat_sh)
    nst_out = denmat_sh[0].num_of_cols

    ave_en_sh = 0.0
    ave_en_se = 0.0
    den = 1.0/float(ntraj)

    for i in xrange(ntraj):
        ave_en_se =  ave_en_se + den * (denmat_se[i].real() * Hvib[i].real() ).tr()  # SE
        ave_en_sh =  ave_en_sh + den * (denmat_sh[i].real() * Hvib[i].real() ).tr()  # SH

    return ave_en_sh, ave_en_se





def amplitudes2denmat(coeffs):
    """
 
    Converts the wavefunction amplitudes for all trajectories to the corresponding 
        density matrices

    Args:
        coeffs ( list of CMATRIX(nstates, 1) ): wavefunction amplitudes for all trajectories

    Returns:
        ( list of CMATRIX(nstates, nstate) ): density matrices for all trajectory

    """

    ntraj = len(coeffs)
    denmat = []

    for tr in xrange(ntraj):
        denmat.append( coeffs[tr] * coeffs[tr].H() )

    return denmat



def denmat2prob(P):
    """Convert the density matrix to the populations

    Args:
        P ( CMATRIX(N, N) ): Density matrix

    Returns:
        ( list of doubles ): populations of all states

    """
    nst = P.num_of_cols
    prob = [0.0] * nst

    for i in xrange(nst):
        prob[i] = P.get(i,i).real

    return prob


