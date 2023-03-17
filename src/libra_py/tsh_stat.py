#*********************************************************************************                     
#* Copyright (C) 2019-2022 Alexey V. Akimov                                                   
#* Copyright (C) 2016-2019 Kosuke Sato, Alexey V. Akimov                                                   
#*                                                                                                     
#* This file is distributed under the terms of the GNU General Public License                          
#* as published by the Free Software Foundation, either version 3 of                                   
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
import numpy as np

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

import util.libutil as comn

from . import units
from . import probabilities
from libra_py import data_conv

def compute_etot(ham, p, Cdia, Cadi, projectors, iM, rep):
    """Computes the Ehrenfest potential energy

    This function computes the average kinetic, potential, and total
    energies for an ensemble of trajectories - according to the Ehrenfest recipe

    Args: 
        ham ( nHamiltonian ): object that handles Hamiltonian-related calculations with many trajectories
        p ( MATRIX(ndof, ntraj) ): nuclear momenta of multiple trajectories
        Cdia ( CMATRIX(ndia, ntraj) ): amplitudes of diabatic states in the TD wavefunction expansion
        Cadi ( CMATRIX(ndia, ntraj) ): amplitudes of adiabatic states in the TD wavefunction expansion
        projectors ( list of CMATRIX(nst, nst)): dynamically-consistent corrections
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

    for traj in range(0,ntraj):

        if rep==0:
            pop_submatrix(Cdia, C, Py2Cpp_int(list(range(0,nst))), Py2Cpp_int([traj]))    
            epot.append( ham.Ehrenfest_energy_dia(C, Py2Cpp_int([0,traj])).real )
            Epot = Epot + epot[traj]
        elif rep==1:

            pop_submatrix(Cadi, C, Py2Cpp_int(list(range(0,nst))), Py2Cpp_int([traj]))    
            C = projectors[traj] * C  # dyn-const -> raw
            epot.append( ham.Ehrenfest_energy_adi(C, Py2Cpp_int([0,traj])).real )
            Epot = Epot + epot[traj]

        tmp = 0.0
        for dof in range(0,ndof):
            tmp = tmp + 0.5 * iM.get(dof, 0) * (p.get(dof, traj) ** 2)
        ekin.append(tmp)
        Ekin = Ekin + ekin[traj]

    Ekin = Ekin / float(ntraj)
    Epot = Epot / float(ntraj)
    Etot = Ekin + Epot

    # Variances:
    dEkin, dEpot = 0.0, 0.0
    for traj in range(0,ntraj):
        dEkin = dEkin + (ekin[traj] - Ekin)**2
        dEpot = dEpot + (epot[traj] - Epot)**2

    dEtot = dEkin + dEpot

    dEkin = math.sqrt(dEkin/ float(ntraj))
    dEpot = math.sqrt(dEpot/ float(ntraj))
    dEtot = math.sqrt(dEtot/ float(ntraj))
    

    return Ekin, Epot, Etot, dEkin, dEpot, dEtot




def compute_etot_tsh(ham, p, Cdia, Cadi, projectors, act_states, iM, rep):
    """Compute the adiabatic potential energy

    This function computes the average kinetic, potential, and total
    energies for an ensemble of trajectories - according to the TSH recipe

    Args: 
        ham ( nHamiltonian ): object that handles Hamiltonian-related calculations with many trajectories
        p ( MATRIX(ndof, ntraj) ): nuclear momenta of multiple trajectories
        Cdia ( CMATRIX(ndia, ntraj) ): amplitudes of diabatic states in the TD wavefunction expansion
        Cadi ( CMATRIX(nadi, ntraj) ): amplitudes of adiabatic states in the TD wavefunction expansion
        projectors ( list of CMATRIX(nst, nst)): dynamically-consistent corrections
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

    #tsh_indx2vec(ham, states, act_states)
    states = tsh_indx2ampl(act_states, nst)

    for traj in range(0,ntraj):

        pop_submatrix(states, C, Py2Cpp_int(list(range(0,nst))), Py2Cpp_int([traj]))      

        if rep==0:
            epot.append( ham.Ehrenfest_energy_dia(C, Py2Cpp_int([0,traj])).real )
            Epot = Epot + epot[traj]
        elif rep==1:
            # C is supposed to be a dyn-consistent amplitudes, so we'd need to convert
            # the Hamiltonian into it, but since the Ehrenfest energy is invariant w.r.t. 
            # the choice of raw/dynamically-consystent rep, we better convert C
            C = projectors[traj] * C  # dyn-const -> raw
            epot.append( ham.Ehrenfest_energy_adi(C, Py2Cpp_int([0,traj])).real )
            Epot = Epot + epot[traj]

        tmp = 0.0
        for dof in range(0,ndof):
            tmp = tmp + 0.5 * iM.get(dof, 0) * (p.get(dof, traj) ** 2)
        ekin.append(tmp)
        Ekin = Ekin + ekin[traj]

    Ekin = Ekin / float(ntraj)
    Epot = Epot / float(ntraj)
    Etot = Ekin + Epot

    # Variances:
    dEkin, dEpot, dEtot = 0.0, 0.0, 0.0
    for traj in range(0,ntraj):
        dEkin = dEkin + (ekin[traj] - Ekin)**2
        dEpot = dEpot + (epot[traj] - Epot)**2
        dEtot = dEtot + (ekin[traj] + epot[traj] - Etot)**2

    dEkin = math.sqrt(dEkin/ float(ntraj))
    dEpot = math.sqrt(dEpot/ float(ntraj))
    dEtot = math.sqrt(dEtot/ float(ntraj))
    

    return Ekin, Epot, Etot, dEkin, dEpot, dEtot




def compute_dm(ham, Cdia, Cadi, projectors, rep, lvl, isNBRA=0):
    """

    Compute the trajectory-averaged density matrices in diabatic
    or adiabatic representations

    Args: 
        ham ( nHamiltonian ): object that handles Hamiltonian-related calculations with many trajectories
        Cdia ( CMATRIX(ndia, ntraj) ): amplitudes of diabatic states in the TD wavefunction expansion
        Cadi ( CMATRIX(ndia, ntraj) ): amplitudes of adiabatic states in the TD wavefunction expansion
        projectors ( list of CMATRIX(nst, nst)): dynamically-consistent corrections   
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

        isNBRA ( int ): The flag for NBRA type calculations:
            - 0: the Hamiltonian related properties are computed for all of the trajectories [default]
            - 1: the Hamiltonian related properties are computed only for one trajectory 

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

    # Dynamically-consistent
    dm_dia, dm_adi = CMATRIX(ndia, ndia), CMATRIX(nadi, nadi)

    # Raw
    dm_dia_raw, dm_adi_raw = CMATRIX(ndia, ndia), CMATRIX(nadi, nadi)


    for traj in range(0,ntraj):

        # In the NBRA case - we are going to use the S and U only for the
        # first Hamiltonian child, since other children have not Hamiltonians         
        ham_indx = traj
        if isNBRA==1: 
            ham_indx = 0
       
        indx = None
        if lvl==0:
            indx = Py2Cpp_int([0])
        elif lvl==1:
            indx = Py2Cpp_int([0,ham_indx])
    
        S = ham.get_ovlp_dia(indx)
        U = ham.get_basis_transform(indx) 

        # However, in both NBRA and non-NBRA cases, we have all the trajectories        
        if rep==0:
            c = Cdia.col(traj)
            tmp = S * (c * c.H()) * S
    
            # Dia dyn-consistent and dia raw
            dm_dia = dm_dia + tmp
            dm_dia_raw = dm_dia_raw + tmp
    
            # Adi raw                     
            tmp = U.H() * tmp * U
            dm_adi_raw = dm_adi_raw + tmp

            # Adi dyn-consistent    
            dm_adi = dm_adi + tmp
                   
    
        elif rep==1:
            # Raw
            c = Cadi.col(traj)            
            tmp = c * c.H()
                        
            # Adi raw
            dm_adi_raw = dm_adi_raw + tmp
            
            # Adi dyn-consistent
            dm_adi = dm_adi + tmp
            #dm_tmp_raw = projectors[traj] * dm_tmp * projectors[traj].H()
            #dm_adi_raw = dm_adi_raw + dm_tmp_raw
    
            # Dia raw and dyn-consistent
            su = S * U
            tmp =  su * tmp * su.H()
            dm_dia = dm_dia + tmp 
            dm_dia_raw = dm_dia_raw + tmp  

    
    dm_dia = dm_dia / float(ntraj)        
    dm_adi = dm_adi / float(ntraj)
    dm_dia_raw = dm_dia_raw / float(ntraj)        
    dm_adi_raw = dm_adi_raw / float(ntraj)

    return dm_dia, dm_adi, dm_dia_raw, dm_adi_raw


def compute_sh_statistics(nstates, istate, projectors, isNBRA=0):
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
        isNBRA ( int ): The flag for NBRA type calculations:
            - 0: the Hamiltonian related properties are computed for all of the trajectories [ default ]
            - 1: the Hamiltonian related properties are computed only for one trajectory 

    Returns: 
        MATRIX(nstates, 1): coeff_sh: The list containing the average 
            population of each quantum state. The length of the list is equal to the 
            total number of quantum states considered, ```nstates```
    """

    ntraj = len(istate)    
    coeff_sh = MATRIX(nstates, 1)        
    coeff_sh_raw = MATRIX(nstates, 1)        

    if isNBRA==1:        

        pops_numpy = np.zeros((nstates, ntraj))
        for i in range(ntraj):
            pops_numpy[istate[i], i] = 1

        coeff_sh_numpy = np.average(pops_numpy, axis=1).reshape(nstates,1)
        coeff_sh = data_conv.nparray2MATRIX(coeff_sh_numpy)        
        coeff_sh_raw = coeff_sh

    else:

        for i in range(0,ntraj):

            st = istate[i]
            #print(st)
            pop = CMATRIX(nstates, nstates)
            pop.set(st, st, 1.0+0.0j)
            #projectors[i].show_matrix()
            pop_raw = projectors[i] * pop * projectors[i].H()
            
            for j in range(nstates):
                coeff_sh.add(j, 0, pop.get(j,j).real)
                coeff_sh_raw.add(j, 0, pop_raw.get(j,j).real)
                
        coeff_sh.scale(-1,0, 1.0/float(ntraj))
        coeff_sh_raw.scale(-1,0, 1.0/float(ntraj))
 
    return coeff_sh, coeff_sh_raw



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

    for j in range(0,ntraj): # for all trajectories
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
    for traj in range(0,ntraj): # for all trajectories
        sh_pops[ el[traj].istate ] += f

        for st1 in range(0,nstat):
            se_pops[ st1 ] += f * el[traj].rho(st1,st1).real

            for st2 in range(0,nstat):
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

    for i in range(0,ntraj):
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

    for i in range(0,ntraj):
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

    for tr in range(0,ntraj):
        denmat.append( coeffs[tr] * coeffs[tr].H() )

    return denmat


def pops2denmat(pops):
    """
 
    Converts the populations (vector) of all states for all trajectories to the corresponding 
        density matrices (matrix). This is just a convenience function

    Args:
        pops ( list of CMATRIX(nstates, 1) ): state populations for all trajectories

    Returns:
        ( list of CMATRIX(nstates, nstate) ): density matrices for all trajectories

    """

    ntraj = len(pops)
    nstates = pops[0].num_of_rows
    denmat = []

    for tr in range(0,ntraj):
        denmat.append( CMATRIX(nstates, nstates) )
        for i in range(0,nstates):
            denmat[tr].set(i,i, pops[tr].get(i,0))

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

    for i in range(0,nst):
        prob[i] = P.get(i,i).real

    return prob





def probabilities_1D_scattering(q, states, nst, params):
    """Computes the scattering probabilities in 1D

    Args:
        _q ( MATRIX(nnucl, ntraj) ): coordinates of the "classical" particles [units: Bohr]
        states ( intList, or list of ntraj ints ): the quantum state of each trajectory
        nst ( int ): the number of possible quantum states in the problem
        params ( dictionary ): parameters of the simulation, should contain
 
            * **params["act_dof"]** ( int ): index of the nuclear DOF that is considered active (scattering coord)
            * **params["left_boundary"] ( double ): the beginning of the reflected particles counter [units: Bohr]
            * **params["right_boundary"] ( double ): the beginning of the transmitted particles counter [units: Bohr]

    Returns:
        tuple: ( pop_refl, pop_transm ): where

            * pop_refl ( MATRIX(nst, 1) ): probabilities of reflection on each state
            * pop_transm ( MATRIX(nst, 1) ): probabilities of transmission on each state

    """

    critical_params = [  ] 
    default_params = {"act_dof":0, "left_boundary":-10.0, "right_boundary":10.0 }
    comn.check_input(params, default_params, critical_params)


    act_dof = params["act_dof"]
    left_boundary = params["left_boundary"]
    right_boundary = params["right_boundary"]


    ntraj = len(states)

    pop_transm = MATRIX(nst, 1)  # transmitted
    pop_refl = MATRIX(nst, 1)    # reflected

    ntransm, nrefl = 0.0, 0.0
    for traj in range(0,ntraj):

        if q.get(act_dof, traj) < left_boundary:
            pop_refl.add(states[traj], 0, 1.0)
            nrefl += 1.0

        if q.get(act_dof, traj) > right_boundary:
            pop_transm.add(states[traj], 0, 1.0)
            ntransm += 1.0         

    ntot = ntransm + nrefl 
    if ntot > 0.0:
        pop_transm = pop_transm / ntot
        pop_refl = pop_refl / ntot

    return pop_refl, pop_transm
