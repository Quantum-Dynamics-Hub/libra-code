#*********************************************************************************
#* Copyright (C) 2020 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 3 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
"""
.. module:: psi4_methods
   :platform: Unix, Windows
   :synopsis: This module implements functions for interfacing Psi4 and Libra codes

.. moduleauthor:: 
       Alexey V. Akimov 
  
"""


import os
import sys
import math
import cmath
import re

# Uncomment if you have psi4 installed
#from psi4 import energy, gradient
import numpy as np

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
import util.libutil as comn

from . import units
from . import scan


class tmp:
    pass

def run_psi4_adi(q, params_, full_id):
    """
   
    This function executes the Psi4 quantum chemistry calculations and 
    returns the key properties needed for dynamical calculations.

    Args: 
        q ( MATRIX(ndof, ntraj) ): coordinates of the particles [ units: Bohr ]
        params ( dictionary ): model parameters
 
            * **params["labels"]** ( list of strings ): the labels of atomic symbolc - for all atoms,
                and in a order that is consistent with the coordinates (in triples) stored in `q`.
                The number of this labels is `natoms`, such that `ndof` = 3 * `natoms`. [ Required ]
            * **params["nstates"]** ( int ): the total number of electronic states 
                in this model [ default: 1 - just the ground state ]
            * **params["grad_method_gs"]** ( string ):  the name of the methodology to compute the 
                energy and gradients on the ground state [ defaut: "ccsd/sto-3g" ]
                Examples: 
                  "pbe/sto-3g", "mp2/aug-cc-pVDZ", "ccsd/aug-cc-pVDZ" # ground state energies, gradients                   
            * **params["grad_method_ex"]** ( string ):  the name of the methodology to compute the 
                energy and gradients on the excited states [ defaut: "eom-ccsd/sto-3g" ]
                Examples:                                     
                  "eom-ccsd/aug-cc-pVDZ", # excited state energies, gradients                  
                If you need just excited state energies (but not gradients), consider: 
                "cisd/aug-cc-pVDZ", adc/aug-cc-pVDZ                
            * **params["charge"]** ( int ): the total charge of the system [ default: 0 ]
            * **params["spin_multiplicity"]** ( int ): the total spin multiplicity [ default: 1 - singlet ]
            * **params["options"]** ( dictionary ): additional parameters of calculations [ default: empty ]
                Examples: 
                  - {} - noting extra
                  - {'reference':'rohf'}, 
                  - {'roots_per_irrep':[3, 0, 0, 0], 'prop_root':1, 'reference':'rohf'}
                  - {'num_roots':3, 'follow_root':2, 'reference':'rohf'} - for state-resolved gradients
            * **params["verbosity"]** ( int ): the level of output of the execution-related 
                information [ default : 0]
                  
        full_id ( intList ): the "path" to the Hamiltonian in the Hamiltonian's hierarchy. Usually, 
            this is Py2Cpp_int([0, itraj]) - the index of the trajectory in a swarm of trajectories
            
    Returns:       
        PyObject: obj, with the members:

            * obj.ham_adi ( CMATRIX(nstates,nstates) ): adiabatic Hamiltonian
            * obj.hvib_adi ( CMATRIX(nstates,nstates) ): vibronic Hamiltonian in the adiabatic basis
            * obj.d1ham_adi ( list of ndof CMATRIX(nstates, nstates) objects ): 
                derivatives of the adiabatic Hamiltonian w.r.t. the nuclear coordinate            
 
    """

    # Make a copy of the input parameters dictionary
    params = dict(params_)
    
    # Defaults    
    critical_params = [ "labels" ] 
    default_params = { "nstates":1,
                       "grad_method_gs":"ccsd/sto-3g", 
                       "grad_method_ex":"eom-ccsd/sto-3g", 
                       "charge":0, "spin_multiplicity":1,                       
                       "options":{},
                       "verbosity":0
                     }
    comn.check_input(params, default_params, critical_params)
    
    # Extract the key variables
    grad_method_gs = params["grad_method_gs"]
    grad_method_ex = params["grad_method_ex"]
    charge = params["charge"]
    spin_multiplicity = params["spin_multiplicity"]
    labels = params["labels"]    
    nstates = params["nstates"]
    options = params["options"]
    verbosity = params["verbosity"]
    natoms = len(labels)
    ndof = 3 * natoms
    
    
    obj = tmp()
    obj.ham_adi = CMATRIX(nstates, nstates)    
    obj.hvib_adi = CMATRIX(nstates, nstates)            
    obj.d1ham_adi = CMATRIXList();            
    for idof in range(ndof):        
        obj.d1ham_adi.append( CMATRIX(nstates, nstates) )
  

    Id = Cpp2Py(full_id)
    indx = Id[-1]

    # Setup and execute the PSI4 calculations
    psi4.core.set_output_file('tmp.dat',False)    
    coords_str = scan.coords2xyz(labels, q, indx)

    mol = psi4.geometry(F"""
    {charge} {spin_multiplicity}
    {coords_str }

    units bohr
    """)

        
    for istate in range(nstates):

        E, grad = None, None
        if istate == 0:
                        
            # Compute the energy 
            if verbosity >0:
                print("Doing state 0 energy")
                
            E, wfc = energy(grad_method_gs, molecule=mol, return_wfn=True)

            # Compute force at the converged density
            if verbosity >0:
                print("Doing state 0 gradient")
            grad = np.asarray(gradient(grad_method_gs, ref_wfn=wfc))            
            
        else:
    
            opt = dict(options)
            opt.update({ 'prop_root':istate, 'roots_per_irrep':[3, 0, 0, 0],  'reference':'rohf'  })
            psi4.set_options(opt)

            # Compute the energy 
            if verbosity >0:
                print(F"Doing state {istate} energy")
                
            E, wfc = energy(grad_method_ex, molecule=mol, return_wfn=True)

            # Compute force at the converged density
            if verbosity >0:
                print(F"Doing state {istate} gradient")
            
            grad = np.asarray(gradient(grad_method_ex, ref_wfn=wfc))

        
        obj.ham_adi.set(istate, istate, E * (1.0+0.0j) )
        obj.hvib_adi.set(istate, istate, E * (1.0+0.0j) )        
        for iatom in range(natoms):        
            obj.d1ham_adi[3 * iatom + 0].set(istate, istate, grad[iatom, 0] * (1.0+0.0j) )
            obj.d1ham_adi[3 * iatom + 1].set(istate, istate, grad[iatom, 1] * (1.0+0.0j) )
            obj.d1ham_adi[3 * iatom + 2].set(istate, istate, grad[iatom, 2] * (1.0+0.0j) )
        
                        
    return obj
