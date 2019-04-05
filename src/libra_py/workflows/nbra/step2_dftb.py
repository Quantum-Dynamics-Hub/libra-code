#***********************************************************
# * Copyright (C) 2019 Alexey V. Akimov
# * This file is distributed under the terms of the
# * GNU General Public License as published by the
# * Free Software Foundation; either version 3 of the
# * License, or (at your option) any later version.
# * http://www.gnu.org/copyleft/gpl.txt
#***********************************************************/

"""
.. module:: step2_dftb
   :platform: Unix, Windows
   :synopsis: This module implements functions for doing NAC calculaitons in the 1-electron
       basis of KS orbitals, with the DFTB+ package

.. moduleauthor:: Alexey V. Akimov

"""

import os
import sys

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from libra_py import DFTB_methods
from libra_py import units
#import libra_py.common_utils as comn
import util.libutil as comn




def do_step(i, params):
    """

    Runs a single-point SCF calculation for a given geometry along a trajectory

    Args:
        i ( int ): index of the time step to be used from the trajectory file
        params ( dictionary ): the control parameters of the simulation
        
        * **params["EXE"]** ( string ): path to the DFTB+ executable [ default: dftb+ ]
        * **params["mo_active_space"]** ( list of ints or None ): indices of the MOs we care about 
            The indexing starts from 0, not 1! If set to None - all MOs will be returned. [default: None]
        * **params["md_file"]** ( string ): the name of the xyz file containing the trajectory - the 
            file should be in the xyz format produced by the DFTB+ program. [default: "md.xyz"]
        * **params["sp_gen_file"]** ( string ): the name of the .gen file that is listed in the 
            DFTB+ input file and contains the geometry of the system (the content of this file will be 
            updated for every i value). [default: "x1.gen" ]
        * **params["syst_spec"]** ( string ): the string that is a part of the DFTB+ .gen file and defines
            whether the system is a non-periodic/cluster ("C") or periodic ("S"). [default: "C"]
        * **params["scf_in_file"]** ( string ): the name of the file containing the template for running regular
            SCF calculations for a single-point DFTB+ calculations. 
 
            - It should use the geometry file defined by params["sp_gen_file"]. 

            [default: "dftb_in_ham1.hsd"]

        * **params["hs_in_file"]** ( string ): the name of the file containing the template for running 
            calculations that construct H and S matrices and print them out by the DFTB+ calculations.  

            - It should use the geometry file defined by params["sp_gen_file"]. 
            - It should have the section: "ReadInitialCharges = Yes" to use the previously-converged charge density
            - It should have the section: "WriteHS = Yes" to initialize the writing of the H and S matrices
 
            [default: "dftb_in_ham2.hsd"]

    Returns:
        tuple: (Ei, MOi, Hi, Si), where:
        
            * Ei ( CMATRIX(M, M) ), the matrix of the converged Hamiltonian eigenvalues (at given geometry)
                Here, M = len(params["mo_active_space"]) - we output only the MO energies that are of interest to us

            * MOi ( CMATRIX(A, M) ), the matrix of the converged Hamiltonian eigenvalues (at given geometry)
                Here, M = len(params["mo_active_space"]) - we output only the MO energies that are of interest to us.
                A - is the number of AOs in this calculation

            * Hi (list of CMATRIX(A, A) ): the Hamiltonian matrices in the AO basis, for each k-point

            * Si (list of CMATRIX(A, A) ): the overlap matrices in the AO basis, for each k-point

    """


    # Now try to get parameters from the input
    critical_params = [ ] 
    default_params = { "EXE":"dftb+",
                       "mo_active_space":None,
                       "md_file":"md.xyz", "sp_gen_file":"x1.gen", "syst_spec":"C" ,
                       "scf_in_file":"dftb_in_ham1.hsd",
                       "hs_in_file":"dftb_in_ham2.hsd"
                     }
    comn.check_input(params, default_params, critical_params)

    
    # Get the parameters
    EXE = params["EXE"]
    md_file = params["md_file"]
    sp_gen_file = params["sp_gen_file"]
    syst_spec = params["syst_spec"]
    scf_in_file = params["scf_in_file"]
    hs_in_file = params["hs_in_file"]

    
    # Make an input file for SP calculations 
    DFTB_methods.xyz_traj2gen_sp(md_file, sp_gen_file, i, syst_spec)

    # Run SCF calculations and generate the charge density for a converged calculations
    # The file x1.gen is used as a geometry
    os.system("cp %s dftb_in.hsd" % scf_in_file )
    os.system( "%s" % EXE )

    # Just generate the Hamiltonian corresponding to the converged density matrix
    os.system("cp %s dftb_in.hsd" % hs_in_file )
    os.system( "%s" % EXE )

    # Get the Hamiltonian    
    Hi = DFTB_methods.get_dftb_matrices("hamsqr1.dat")
    Si = DFTB_methods.get_dftb_matrices("oversqr.dat")


    # Get the dimensions
    ao_sz = Hi[0].num_of_cols
    mo_sz = ao_sz
    ao_act_sp = range(0, ao_sz)
    mo_act_sp = range(0, mo_sz)

    if params["mo_active_space"] != None:
        mo_sz = len(params["mo_active_space"])
        mo_act_sp = list(params["mo_active_space"])


    # Extract the sub-matrix of interest
    H_sub = CMATRIX(ao_sz, mo_sz)
    pop_submatrix(Hi[0], H_sub, ao_act_sp, mo_act_sp)  # list element #0 = gamma-point

    # Get the orbitals
    Ei = CMATRIX(mo_sz, mo_sz)
    MOi = CMATRIX(ao_sz, mo_sz)
    solve_eigen(H_sub, Si[0], Ei, MOi, 0)     # list element #0 = gamma-point
    
    return Ei, MOi, Hi, Si



def do_ovlp(i, params):
    """

    Compute the overlap matrix in the AO basis for two geometries, i and i+1

    Args:
        i ( int ): index of the time step to be used from the trajectory file
        params ( dictionary ): the control parameters of the simulation
        
        * **params["EXE"]** ( string ): path to the DFTB+ executable [ default: dftb+ ]
        * **params["md_file"]** ( string ): the name of the xyz file containing the trajectory - the 
            file should be in the xyz format produced by the DFTB+ program. [default: "md.xyz"]
        * **params["ovlp_gen_file"]** ( string ): the name of the .gen file that is listed in the 
            DFTB+ input file and contains the geometry of the system at two time steps
            (the content of this file will be updated for every i value). [default: "x2.gen" ]
        * **params["syst_spec"]** ( string ): the string that is a part of the DFTB+ .gen file and defines
            whether the system is a non-periodic/cluster ("C") or periodic ("S"). [default: "C"]
        * **params["ovlp_in_file"]** ( string ): the name of the file containing the template for running 
            calculations that construct H and S matrices and print them out by the DFTB+ calculations.  

            - It should use the geometry file defined by params["ovlp_gen_file"]. 
            - It should have the section: "WriteHS = Yes" to initialize the writing of the H and S matrices
 
            [default: "dftb_in_overlaps.hsd"]


    Returns:
        CMATRIX(A, A): the matrix of the AO overlaps for two geometries, where A - is the size of the AO basis

    """

    # Now try to get parameters from the input
    critical_params = [ ] 
    default_params = { "EXE":"dftb+",
                       "md_file":"md.xyz", "ovlp_gen_file":"x2.gen", "syst_spec":"C" ,
                       "ovlp_in_file":"dftb_in_overlaps.hsd"
                     }
    comn.check_input(params, default_params, critical_params)
    
    # Get the parameters
    EXE = params["EXE"]
    md_file = params["md_file"]
    ovlp_gen_file = params["ovlp_gen_file"]
    syst_spec = params["syst_spec"]
    ovlp_in_file = params["ovlp_in_file"]

        
    # Make an input file for the overlap calculations 
    DFTB_methods.xyz_traj2gen_ovlp(md_file, ovlp_gen_file, i, i+1, syst_spec)
    
    # Run SCF calculations and generate the charge density for a converged calculations
    # The file x2.gen is used as a geometry
    os.system("cp %s dftb_in.hsd" % ovlp_in_file)
    os.system( "%s" % EXE )

    # Get the Hamiltonian    
    Sbig = DFTB_methods.get_dftb_matrices("oversqr.dat")
    norbs = Sbig[0].num_of_cols/2
    
    act_sp1 = range(0,norbs)
    act_sp2 = range(norbs,2*norbs)
    S = CMATRIX(norbs, norbs)
    pop_submatrix(Sbig[0], S, act_sp1, act_sp2)
    
    return S
    
    

    
def run_step2(params):
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
    E_curr, U_curr, Hao_curr, Sao_curr = do_step(isnap, params)
    
    for i in xrange(isnap+1, fsnap-1):
        E_next, U_next, Hao_next, Sao_next = do_step(i, params)
        S = do_ovlp(i, params)

        S.real().show_matrix("res/AOS_%i_re" % (i) )
        
        TDM = U_curr.H() * S * U_next
        Hvib = 0.5*(E_curr + E_next) - (0.5j/dt) * ( TDM - TDM.H() )

        # Overlaps
        s = 0.5 * (U_curr.H() * Sao_curr[0] * U_curr  +  U_next.H() * Sao_next[0] * U_next)
        s.real().show_matrix("res/S_%i_re" % (i) )
        s.imag().show_matrix("res/S_%i_im" % (i) )


        # Time-overlaps
        TDM.real().show_matrix("res/St_%i_re" % (i) )
        TDM.imag().show_matrix("res/St_%i_im" % (i) )
                
        # Vibronic Hamiltonians
        Hvib.real().show_matrix("res/hvib_%i_re" % (i) )
        Hvib.imag().show_matrix("res/hvib_%i_im" % (i) )

        
        # Current becomes the old 
        E_curr = CMATRIX(E_next)
        U_curr = CMATRIX(U_next)
        

