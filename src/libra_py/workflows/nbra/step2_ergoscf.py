#***********************************************************
# * Copyright (C) 2019 Alexey V. Akimov
# * This file is distributed under the terms of the
# * GNU General Public License as published by the
# * Free Software Foundation; either version 3 of the
# * License, or (at your option) any later version.
# * http://www.gnu.org/copyleft/gpl.txt
#***********************************************************/

"""
.. module:: step2_ergoscf
   :platform: Unix, Windows
   :synopsis: This module implements functions for doing NAC calculaitons in the 1-electron
       basis of MOs, with the ErgoSCF package

.. moduleauthor:: Alexey V. Akimov

"""

import os
import sys

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

import util.libutil as comn
import libra_py.packages.ergo.methods as ERGO_methods
from libra_py import units




def do_step(i, params, run):
    """

    Runs a single-point SCF calculation for a given geometry along a trajectory

    Args:
        i ( int ): index of the time step to be used from the trajectory file
        params ( dictionary ): the control parameters of the simulation
        
            * **params["EXE"]** ( string ): path to the ErgoSCF executable [ default: ergo ]
            * **params["md_file"]** ( string ): the name of the xyz file containing the trajectory - the 
                file should be in the general xyz format. [default: "md.xyz"]
            * **params["mo_active_space"]** ( list of ints or None ): indices of the MOs we care about 
                The indexing starts from 0, not 1! If set to None - all MOs will be returned. [default: None]
            * **params["mo_indexing_convention"]** ( string ): orbital indexing convention:
                - "rel" : relative, with respect to the HOMO position: HOMO = 0, HOMO-1 = -1, LUMO = 1, etc. [default]
                - "abs" : absolute, with respect to all possible MOs of the system 
            * **params["direct_MO"]** ( int ): a flag to decide to use orbitals as read from 
                the ErgoSCF output:
                - 0 : read in the Fock matrices and diagonalize them internally
                - 1 : read the MOs produced by ErgoSCF, this is the option to really go 
                    with the linear-scaling costs [ default ]
            * **params["spinpolarized"]** ( int ): whether the ErgoSCF calculations (and hence the
                following post-processing) are done withing the restricted or unrestricted formulation:
                - 0 : non-spin-polarized (restricted) [ default ] 
                - 1 : spin-polarized (unrestricted)

        run (Python function ): the function that defines the ErgoSCF input generation - the user 
            has to define all the control parameters in it, to be able to run the ErgoSCF calculations

            Example:

                In the example below, the outermost quotes should be tripled 

                Note: the function should follw the signature shown here

                def run(EXE, COORDS):
                    inp = "#!bin/sh
                %s << EOINPUT > /dev/null
                spin_polarization = 0
                molecule_inline
                %sEOF
                basis = "STO-3G"
                use_simple_starting_guess=1
                scf.create_mtx_files_F = 1
                scf.create_mtx_file_S = 1
                XC.sparse_mode = 1
                run "LDA"
                EOINPUT
                " % (EXE, COORDS)
                    return inp


    Returns:
        (list, list): (E_sub, MO_sub), where:
        
            * E_sub ( list of CMATRIX(M, M) ): the matrix of the converged Hamiltonian eigenvalues (at given geometry)
                Here, M = len(params["mo_active_space"]) - we output only the MO energies that are of interest to us.
                E_sub[0] is the alpha component (for both restricted and unrestricted calculations)
                E_sub[1] is the beta component (only for the unrestricted calculations)
                       
            * MO_sub ( list of CMATRIX(A, M) ): the matrix of the converged Hamiltonian eigenvalues (at given geometry)
                Here, M = len(params["mo_active_space"]) - we output only the MO energies that are of interest to us.
                A - is the number of AOs in this calculation
                MO_sub[0] is the alpha component (for both restricted and unrestricted calculations)
                MO_sub[1] is the beta component (only for the unrestricted calculations)


    """

    # Now try to get parameters from the input
    critical_params = [ ] 
    default_params = { "EXE":"ergo",  "md_file":"md.xyz",
                       "mo_active_space":None, "mo_indexing_convention":"rel",
                       "direct_MO":1,  "spinpolarized":0
                     }
    comn.check_input(params, default_params, critical_params)
    
    # Get the parameters
    EXE = params["EXE"]
    md_file = params["md_file"]
    mo_indexing_convention = params["mo_indexing_convention"]
    direct_MO = params["direct_MO"]
    spinpolarized = params["spinpolarized"]

    if direct_MO==1 and mo_indexing_convention=="abs":
        print("WARNING in step2_ergoscf.do_step: \
               Cannot use the direct_MO = 1 with the mo_indexing_convention=\"abs\"\
               Resetting mo_indexing_convention to =\"rel\" ")
        mo_indexing_convention = "rel"

                
    # Make an input file for SP calculations     
    R = ERGO_methods.xyz_traj2gen_sp(md_file, i)

    # Run SCF calculations
    command = run(EXE, R)
    os.system("%s" % (command))

    # Read in the common info
    S = ERGO_methods.get_mtx_matrices("S_matrix.mtx")    

    
    E_sub, MO_sub = [], []

    if direct_MO == 0:

        # Get the dimensions
        ao_sz = S.num_of_cols        
        ao_act_sp = list(range(0, ao_sz))
      
        mo_sz = ao_sz
        mo_act_sp = list(range(0, mo_sz))
    
        if params["mo_active_space"] != None:
            mo_sz = len(params["mo_active_space"])        
            mo_act_sp = list(params["mo_active_space"])

        if spinpolarized==0:

            # Get the last Fock matrix 
            last_indx, last_filename = ERGO_methods.find_last_file("F_matrix_", ".mtx")
            F = ERGO_methods.get_mtx_matrices(last_filename)    
                
            # Solve the eigenvalue problem with the converged Fock matrix
            # get the converged MOs
            E = CMATRIX(ao_sz, ao_sz)
            MO = CMATRIX(ao_sz, ao_sz)
            solve_eigen(F, S, E, MO, 0)  

            # Extract the E sub-matrix
            e_sub = CMATRIX(mo_sz, mo_sz)
            pop_submatrix(E, e_sub, mo_act_sp, mo_act_sp)           
            E_sub.append(e_sub)
    
            # Extract the MO sub-matrix
            mo_sub = CMATRIX(ao_sz, mo_sz)
            pop_submatrix(MO, mo_sub, ao_act_sp, mo_act_sp)  
            MO_sub.append(mo_sub)

        elif spinpolarized==1:

            for suff in ["alpha", "beta"]:

                # Get the last Fock matrix 
                last_indx, last_filename = ERGO_methods.find_last_file(F"F_matrix_{suff}_", ".mtx")
                F = ERGO_methods.get_mtx_matrices(last_filename)    
                
                # Solve the eigenvalue problem with the converged Fock matrix
                # get the converged MOs
                E = CMATRIX(ao_sz, ao_sz)
                MO = CMATRIX(ao_sz, ao_sz)
                solve_eigen(F, S, E, MO, 0)  

                # Extract the E sub-matrix
                e_sub = CMATRIX(mo_sz, mo_sz)
                pop_submatrix(E, e_sub, mo_act_sp, mo_act_sp)           
                E_sub.append(e_sub)
    
                # Extract the MO sub-matrix
                mo_sub = CMATRIX(ao_sz, mo_sz)
                pop_submatrix(MO, mo_sub, ao_act_sp, mo_act_sp)  
                MO_sub.append(mo_sub)



    elif direct_MO == 1:

        if spinpolarized==0:
            [e_a], [nocc, nvirt] = ERGO_methods.read_spectrum_restricted() 
            E = ERGO_methods.energies(e_a, nocc, nvirt, params["mo_active_space"])
            E_sub.append(E)

            [mo_a] = ERGO_methods.read_mo_restricted(nocc, nvirt, params["mo_active_space"]) 
            MO_sub.append(mo_a)

        elif spinpolarized==1:
            [e_a, e_b], [nocc, nvirt] = ERGO_methods.read_spectrum_unrestricted() 
            Ea = ERGO_methods.energies(e_a, nocc, nvirt, params["mo_active_space"])
            Eb = ERGO_methods.energies(e_b, nocc, nvirt, params["mo_active_space"])
            E_sub.append(Ea)
            E_sub.append(Eb)

            [mo_a, mo_b] = ERGO_methods.read_mo_unrestricted(nocc, nvirt, params["mo_active_space"]) 
            MO_sub.append(mo_a)
            MO_sub.append(mo_b)

    return E_sub, MO_sub



def do_ovlp(i, params, run):
    """

    Compute the overlap matrix in the AO basis for two geometries, i and i+1

    Args:
        i ( int ): index of the time step to be used from the trajectory file
        params ( dictionary ): the control parameters of the simulation
        
            * **params["EXE"]** ( string ): path to the ERGOSCF executable [ default: ergo ]
            * **params["md_file"]** ( string ): the name of the xyz file containing the trajectory - the 
                file should be in the general xyz format. [default: "md.xyz"]

        run (Python function ): the function that defines the ErgoSCF input generation - the user 
            has to define all the control parameters in it, to be able to run the ErgoSCF calculations

            Example:

                In the example below, the outermost quotes should be tripled 

                Note: the function should follw the signature shown here

                def run(EXE, COORDS):
                    inp = "#!bin/sh
                %s << EOINPUT > /dev/null
                spin_polarization = 0
                molecule_inline
                %sEOF
                basis = "STO-3G"
                scf.create_mtx_file_S = 1
                scf.create_mtx_files_S_and_quit = 1
                XC.sparse_mode = 1
                run "LDA"
                EOINPUT
                " % (EXE, COORDS)
                    return inp

    Returns:
        (CMATRIX(A, A), CMATRIX(A, A), CMATRIX(A, A)): (S11, S22, S12), where: 

            * S11 - the matrices of the AO overlaps for the first geometry with itself (normal AO overlaps)
            * S22 - the matrices of the AO overlaps for the second geometry with itself (normal AO overlaps)
            * S12 - the matrices of the AO overlaps for the two geometries at the adjacent geometries (for TDM)

            where A - is the size of the AO basis

    """

    # Now try to get parameters from the input
    critical_params = [ ] 
    default_params = { "EXE":"ergo", "md_file":"md.xyz"   }
    comn.check_input(params, default_params, critical_params)
    
    # Get the parameters
    EXE = params["EXE"]
    md_file = params["md_file"]
    
        
    # Make an input file for the overlap calculations 
    R = ERGO_methods.xyz_traj2gen_ovlp(md_file, i, i+1)
    
    # Run SCF calculations
    command = run(EXE, R)
    os.system("%s" % (command))
    
    # Get the overlap matrix
    S = ERGO_methods.get_mtx_matrices("S_matrix.mtx")

    # In this function we run the calculation for a super-system 
    # (a dimer of the actual chemical one), so the number of orbitals is 
    # defined like this:
    norbs = int(S.num_of_cols/2)

    """
    BEWARE: The following orbital ordering works only with the 
    following changes in the ErgoSCF code:

    ergo_scripted.cc
    const int skip_sort_shells = 1;  /** AVA on 4/28/2019 */

    Great thanks to Elias Rudberg for this suggestion!
    """

    act_sp1 = list(range(0,norbs))
    act_sp2 = list(range(norbs,2*norbs))

    S11 = CMATRIX(norbs, norbs)
    S12 = CMATRIX(norbs, norbs)
    S22 = CMATRIX(norbs, norbs)

    pop_submatrix(S, S11, act_sp1, act_sp1)
    pop_submatrix(S, S12, act_sp1, act_sp2)
    pop_submatrix(S, S22, act_sp2, act_sp2)
    
    return S11, S12, S22


def clean():
    os.system("rm ergoscf.out")
    os.system("rm density.bin")
    os.system("rm F_matrix_*.mtx")

    
def run_step2(params, run1, run2):
    """
    
    Calculate the overlaps, transition dipole moments, and vibronic Hamiltonian matrix elements in the AO basis

    Args:
        params ( dictionary ): the control parameters of the simulation
        
        * **params["dt"]** ( double ): nuclear dynamics timestep - as encoded in the trajectory [ units: a.u., default: 41.0 ]
        * **params["isnap"]** ( int ): initial frame  [ default: 0 ]
        * **params["fsnap"]** ( int ): final frame  [ default: 1 ]
        * **params["out_dir"]** ( string ): the path to the directory that will collect all the results
            If the directory doesn't exist, it will be created  [ default: "res" ]

        * **params["EXE"]** ( string ): path to the ErgoSCF executable [ default: ergo ]
        * **params["md_file"]** ( string ): the name of the xyz file containing the trajectory - the 
            file should be in the general xyz format. [default: "md.xyz"]
        * **params["mo_active_space"]** ( list of ints or None ): indices of the MOs we care about 
            The indexing starts from 0, not 1! If set to None - all MOs will be returned. [default: None]
        * **params["mo_indexing_convention"]** ( string ): orbital indexing convention:
            - "rel" : relative, with respect to the HOMO position: HOMO = 0, HOMO-1 = -1, LUMO = 1, etc. [default]
            - "abs" : absolute, with respect to all possible MOs of the system 
        * **params["direct_MO"]** ( int ): a flag to decide to use orbitals as read from 
            the ErgoSCF output:
            - 0 : read in the Fock matrices and diagonalize them internally
            - 1 : read the MOs produced by ErgoSCF, this is the option to really go 
                with the linear-scaling costs [ default ]
        * **params["spinpolarized"]** ( int ): whether the ErgoSCF calculations (and hence the
            following post-processing) are done withing the restricted or unrestricted formulation:
            - 0 : non-spin-polarized (restricted) [ default ] 
            - 1 : spin-polarized (unrestricted)

        SeeAlso:  the description of the parameters in ```do_ovlp(i, params)``` and in ```do_step(i, params)```

    Returns:
        None: but generates the files (S, St, and hvib) indexed in the range [isnap, fsnap), for 
        instance, if isnap = 0, fsnap = 3, we will have files hvib_0, hvib_1, hvib_2

    """

    critical_params = [ ] 
    default_params = { "dt":1.0*units.fs2au,  "isnap":0, "fsnap":1, "out_dir":"res",
                       "EXE":"ergo",  "md_file":"md.xyz",
                       "mo_active_space":None, "mo_indexing_convention":"rel",
                       "direct_MO":1,  "spinpolarized":0
                     }
    comn.check_input(params, default_params, critical_params)


    # Get the parameters
    dt = params["dt"]
    isnap = params["isnap"]
    fsnap = params["fsnap"]
    out_dir = params["out_dir"]
    spinpolarized = params["spinpolarized"]


    # Create <out_dir> directory if it does not exist yet
    if os.path.isdir(out_dir):
        pass
    else:
        os.system("mkdir %s" % out_dir)


    # Compute
    clean()
    E_prev, U_prev = do_step(isnap, params, run1)
    
    for i in range(isnap+1, fsnap):         

        clean()
        E_curr, U_curr = do_step(i, params, run1)

        clean()
        S11, S22, S12 = do_ovlp(i-1, params, run2)

        spins = []
        if spinpolarized==0:
            spins = [0]
        elif spinpolarized==1:
            spins = [0, 1]


        # Calculations and print out
        for spin in spins:

            TDM = U_prev[spin].H() * S12 * U_curr[spin]
            Hvib = 0.5*(E_prev[spin] + E_curr[spin]) - (0.5j/dt) * ( TDM - TDM.H() )

            # Overlaps
            s = 0.5 * (U_prev[spin].H() * S11 * U_prev[spin]  +  U_curr[spin].H() * S22 * U_curr[spin])
            s.real().show_matrix(F"{out_dir}/S_{spin}{spin}_{i-1}_re" )
            s.imag().show_matrix(F"{out_dir}/S_{spin}{spin}_{i-1}_im" )

            # Time-overlaps
            TDM.real().show_matrix(F"{out_dir}/St_{spin}{spin}_{i-1}_re" )
            TDM.imag().show_matrix(F"{out_dir}/St_{spin}{spin}_{i-1}_im" )
            
            # Vibronic Hamiltonians
            Hvib.real().show_matrix(F"{out_dir}/hvib_{spin}{spin}_{i-1}_re" )
            Hvib.imag().show_matrix(F"{out_dir}/hvib_{spin}{spin}_{i-1}_im" )

        
        # Current becomes the old 
        if spinpolarized==0:
            E_prev[0] = CMATRIX(E_curr[0])
            U_prev[0] = CMATRIX(U_curr[0])

        elif spinpolarized==1:
            E_prev[0] = CMATRIX(E_curr[0])
            U_prev[0] = CMATRIX(U_curr[0])

            E_prev[1] = CMATRIX(E_curr[1])
            U_prev[1] = CMATRIX(U_curr[1])


