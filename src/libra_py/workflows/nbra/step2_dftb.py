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
import libra_py.common_utils as comn




def do_step(i, EXE, norbs):
    """
    Runs a single-point SCF calculation for a given geometry along a trajectory

    Args:
        i ( int ): index of the time step to be used from the trajectory file
        params ( dictionary ): the control parameters of the simulation
        
        * **params["EXE"]** ( string ): path to the DFTB+ executable [ Required! ]
        * **params["active_space"]** ( list of ints ): indices of the MOs to be included in 
            the active space 
    """


    # Now try to get parameters from the input
    critical_params = [ "EXE", "active_space" ] 
    default_params = { "md_file":"md.xyz", "sp_gen_file":"x1.gen", "syst_spec":"C"  }
    comn.check_input(params, default_params, critical_params)

    
    # normal size of the active space (AO basis)
    act_sp = params["active_space"]
    md_file = params["md_file"]
    sp_gen_file = params["sp_gen_file"]
    syst_spec = params["syst_spec"]

    
    # Make an input file for SP calculations 
    DFTB_methods.xyz_traj2gen_sp(md_file, sp_gen_file, i, syst_spec)

    # Run SCF calculations and generate the charge density for a converged calculations
    # The file x1.gen is used as a geometry
    os.system("cp dftb_in_ham1.hsd dftb_in.hsd")
    os.system( "%s" % EXE )

    # Just generate the Hamiltonian corresponding to the converged density matrix
    os.system("cp dftb_in_ham2.hsd dftb_in.hsd")
    os.system( "%s" % EXE )

    # Get the Hamiltonian    
    Hi = DFTB_methods.get_dftb_matrices("hamsqr1.dat", act_sp, act_sp)    
    Si = DFTB_methods.get_dftb_matrices("oversqr.dat", act_sp, act_sp) 
    Ei = CMATRIX(norbs, norbs)
    MOi = CMATRIX(norbs, norbs)
    solve_eigen(Hi[0], Si[0], Ei, MOi, 0)  # gamma-point only at time i
    
    return Ei, MOi


def do_ovlp(i, EXE, norbs):
        
    # Make an input file for the overlap calculations 
    DFTB_methods.xyz_traj2gen_ovlp("md.xyz", "x2.gen", i, i+1, "C")
    
    # Run SCF calculations and generate the charge density for a converged calculations
    # The file x2.gen is used as a geometry
    os.system("cp dftb_in_overlaps.hsd dftb_in.hsd")
    os.system( "%s" % EXE )

    # Get the Hamiltonian    
    act_sp = range(0,2*norbs)  
    Sbig = DFTB_methods.get_dftb_matrices("oversqr.dat", act_sp, act_sp)    
    
    act_sp1 = range(0,norbs)
    act_sp2 = range(norbs,2*norbs)
    S = CMATRIX(norbs, norbs)
    pop_submatrix(Sbig[0], S, act_sp1, act_sp2)
    
    return S
    
    

    
def run_step2(nsteps, EXE, norbs, dt):

    U_curr = CMATRIX(norbs, norbs)
    U_next = CMATRIX(norbs, norbs)
    E_curr = CMATRIX(norbs, norbs)
    E_next = CMATRIX(norbs, norbs)
        
    act_sp1 = range(0,norbs)  # normal size of the active space (AO basis)
    
    E_curr, U_curr = do_step(0, EXE, norbs)
    
    for i in xrange(1,nsteps-1):
        E_next, U_next = do_step(i, EXE, norbs)
        S = do_ovlp(i, EXE, norbs)
        
        TDM = U_curr.H() * S * U_next
        Hvib = 0.5*(E_curr + E_next) - (0.5j/dt) * ( TDM - TDM.H() )
        
        
        Hvib.real().show_matrix("res/hvib_%i_re" % (i) )
        Hvib.imag().show_matrix("res/hvib_%i_im" % (i) )
        
        # Current becomes the old 
        E_curr = CMATRIX(E_next)
        U_curr = CMATRIX(U_next)
        
