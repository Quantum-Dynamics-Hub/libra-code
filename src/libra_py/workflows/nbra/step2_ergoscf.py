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

from libra_py import ERGO_methods
from libra_py import units
import util.libutil as comn



def do_step(i, params, run):
    """

    Runs a single-point SCF calculation for a given geometry along a trajectory

    Args:
        i ( int ): index of the time step to be used from the trajectory file
        params ( dictionary ): the control parameters of the simulation
        
            * **params["EXE"]** ( string ): path to the ERGOSCF executable [ default: ergo ]
            * **params["mo_active_space"]** ( list of ints or None ): indices of the MOs we care about 
                The indexing starts from 0, not 1! If set to None - all MOs will be returned. [default: None]
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
                use_simple_starting_guess=1
                scf.create_mtx_files_F = 1
                scf.create_mtx_file_S = 1
                XC.sparse_mode = 1
                run "LDA"
                EOINPUT
                " % (EXE, COORDS)
                    return inp


    Returns:
        tuple: (E, MO), where:
        
            * E ( CMATRIX(M, M) ), the matrix of the converged Hamiltonian eigenvalues (at given geometry)
                Here, M = len(params["mo_active_space"]) - we output only the MO energies that are of interest to us

            * MO ( CMATRIX(A, M) ), the matrix of the converged Hamiltonian eigenvalues (at given geometry)
                Here, M = len(params["mo_active_space"]) - we output only the MO energies that are of interest to us.
                A - is the number of AOs in this calculation


    """

    # Now try to get parameters from the input
    critical_params = [ ] 
    default_params = { "EXE":"ergo",  "mo_active_space":None,  "md_file":"md.xyz"    }
    comn.check_input(params, default_params, critical_params)
    
    # Get the parameters
    EXE = params["EXE"]
    md_file = params["md_file"]

        
    # Make an input file for SP calculations     
    R = ERGO_methods.xyz_traj2gen_sp(md_file, i)

    # Run SCF calculations
    command = run(EXE, R)
    os.system("%s" % (command))
    
    # Get the last Fock matrix 
    last_indx, last_filename = ERGO_methods.find_last_file("F_matrix_", ".mtx")
    F = ERGO_methods.get_mtx_matrices(last_filename)    
    S = ERGO_methods.get_mtx_matrices("S_matrix.mtx")
    
        
    # Get the dimensions
    ao_sz = F.num_of_cols        
    ao_act_sp = range(0, ao_sz)
    
    mo_sz = ao_sz
    mo_act_sp = range(0, mo_sz)
    
    if params["mo_active_space"] != None:
        mo_sz = len(params["mo_active_space"])        
        mo_act_sp = list(params["mo_active_space"])
    
    
    # Solve the eigenvalue problem with the converged Fock matrix
    # get the converged MOs
    E = CMATRIX(ao_sz, ao_sz)
    MO = CMATRIX(ao_sz, mo_sz)
    solve_eigen(F, S, E, MO, 0)  

    # Extract the E sub-matrix
    E_sub = CMATRIX(mo_sz, mo_sz)
    pop_submatrix(E, E_sub, mo_act_sp, mo_act_sp)  
    
    # Extract the MO sub-matrix
    MO_sub = CMATRIX(ao_sz, mo_sz)
    pop_submatrix(MO, MO_sub, ao_act_sp, mo_act_sp)  
    
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
                use_simple_starting_guess=1
                scf.create_mtx_files_F = 1
                scf.create_mtx_file_S = 1
                XC.sparse_mode = 1
                run "LDA"
                EOINPUT
                " % (EXE, COORDS)
                    return inp

    Returns:
        CMATRIX(A, A): the matrix of the AO overlaps for two geometries, where A - is the size of the AO basis

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
        
    norbs = S.num_of_cols/2    
    act_sp1 = range(0,norbs)
    act_sp2 = range(norbs,2*norbs)
    S_sub = CMATRIX(norbs, norbs)
    pop_submatrix(S, S_sub, act_sp1, act_sp2)
    
    return S_sub    


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

        SeeAlso:  the description of the parameters in ```do_ovlp(i, params)``` and in ```do_step(i, params)```

    """

    critical_params = [ ] 
    default_params = { "dt":1.0*units.fs2au,  "isnap":0, "fsnap":1  }
    comn.check_input(params, default_params, critical_params)


    # Get the parameters
    dt = params["dt"]
    isnap = params["isnap"]
    fsnap = params["fsnap"]


    # Compute
    clean()
    E_curr, U_curr = do_step(isnap, params, run1)
    
    for i in xrange(isnap+1, fsnap-1):         

        clean()
        E_next, U_next = do_step(i, params, run1)

        clean()
        S = do_ovlp(i, params, run2)

        S.real().show_matrix("res/AOS_%i_re" % (i) )
        
        TDM = U_curr.H() * S * U_next
        Hvib = 0.5*(E_curr + E_next) - (0.5j/dt) * ( TDM - TDM.H() )

        # Overlaps
        #s = 0.5 * (U_curr.H() * Sao_curr[0] * U_curr  +  U_next.H() * Sao_next[0] * U_next)
        #s.real().show_matrix("res/S_%i_re" % (i) )
        #s.imag().show_matrix("res/S_%i_im" % (i) )


        # Time-overlaps
        TDM.real().show_matrix("res/St_%i_re" % (i) )
        TDM.imag().show_matrix("res/St_%i_im" % (i) )
                
        # Vibronic Hamiltonians
        Hvib.real().show_matrix("res/hvib_%i_re" % (i) )
        Hvib.imag().show_matrix("res/hvib_%i_im" % (i) )

        
        # Current becomes the old 
        E_curr = CMATRIX(E_next)
        U_curr = CMATRIX(U_next)
        

