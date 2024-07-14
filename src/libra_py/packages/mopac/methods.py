#*********************************************************************************
#* Copyright (C) 2024  Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 3 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
"""
.. module:: methods
   :platform: Unix, Windows
   :synopsis: This module implements functions for interfacing Libra to MOPAC package

.. moduleauthor:: 
       Alexey V. Akimov
  
"""


import os
import sys
import math
import re
import numpy as np

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
import util.libutil as comn

from libra_py import units
from libra_py import scan
from libra_py import regexlib as rgl

import libra_py.packages.cp2k.methods as CP2K_methods


def make_mopac_input(mopac_input_filename, mopac_run_params, labels, coords):
    """
    This function creates an input file for MOPAC package using the 
    parameters passed in the `mopac_input_params` dictionary
    
    Args: 
        * mopac_input_filename ( string ): the name of the input file to create

        * mopac_run_params ( string ): the string containing the specification for the MOPAC run. 
        E.g. one can use: "INDO C.I.=(6,3) CHARGE=0 RELSCF=0.000001 ALLVEC  WRTCONF=0.00"

        * labels (list of stings): element symbols for atoms in the system (N items), e.g.
         ["C", "H", "H", "H", "H"] for methane

        * coords ( MATRIX(3N, 1) ): Cartesian coordinates of all atoms ordered in triples x, y, z [ units: Bohr ]

    Returns:
        None :  just creates the files
        
    """
    
        
    # Create the actual output file
    mopac_input = open(mopac_input_filename,"w");

    mopac_input.write(F"{mopac_run_params}\n\n")

    nat = len(labels) # how many atoms 
    for i in range(nat):
        x = coords.get(3*i+0, 0)/units.Angst
        y = coords.get(3*i+0, 1)/units.Angst
        z = coords.get(3*i+0, 2)/units.Angst
        mopac_input.write(F"{labels[i]}   {x}  1  {y} 1   {z}  1\n")
     
    mopac_input.close()


class tmp:
    pass

def run_mopac_adi(q, params_, full_id):
    """
   
    This function executes the MOPAC quantum chemistry calculations and 
    returns the key properties needed for dynamical calculations.

    Args: 
        q ( MATRIX(ndof, ntraj) ): coordinates of the particle [ units: Bohr ]
        params ( dictionary ): model parameters

            * **params_["labels"]** ( list of strings ): the labels of atomic symbolc - for all atoms,
                and in a order that is consistent with the coordinates (in triples) stored in `q`.
                The number of this labels is `natoms`, such that `ndof` = 3 * `natoms`. [ Required ]                
            * **params["nstates"]** ( int ): the total number of electronic states 
                in this model [ default: 1 - just the ground state]                     
            * **params_["mopac_exe"]** ( string ):  the full path to `the mopac` executable [ defaut: "mopac" ]
            * **params_["mopac_run_params"]** ( string ): the control string to define the MOPAC job
                [default: "INDO C.I.=(6,3) CHARGE=0 RELSCF=0.000001 ALLVEC  WRTCONF=0.00"]
            * **params_["mopac_working_directory"]** ( string ) [ default: "mopac_wd"]
            * **params_["mopac_jobid"]** ( string ) [ default: "job_0000" ]
            * **params_["mopac_input_prefix"]** ( string ) [ default: "input_" ]
            * **params_["mopac_output_prefix"]** ( string ) [ default: "output_" ]
                            
    Returns:       
        PyObject: obj, with the members:

            * obj.ham_adi ( CMATRIX(nstates,nstates) ): adiabatic Hamiltonian             
            * obj.d1ham_adi ( list of ndof CMATRIX(nstates, nstates) objects ): 
                derivatives of the adiabatic Hamiltonian w.r.t. the nuclear coordinate            


        string method_id="HAMILTONIAN:";
        string time_step_id="TIME STEPS:";
        string trajectory_id="TRAJECTORY FILE:";
        string trajectory_directory_id="TRAJECTORY DIRECTORY:";
        string charge_id="CHARGE:";
        string convergence_id="SCF CONVERGENCE CRITERIA:";
        string mopac_prefix_id="MOPAC PREFIX:";
 
    """

    params = dict(params_)
    
    critical_params = [ "labels" ] 
    default_params = { "nstates":1,
                       "mopac_exe":"mopac", 
                       "mopac_run_params":"INDO C.I.=(6,3) CHARGE=0 RELSCF=0.000001 ALLVEC  WRTCONF=0.00",
                       "mopac_working_directory":"mopac_wd",
                       "mopac_jobid":"job_0000",
                       "mopac_input_prefix":"input_", "mopac_output_prefix":"output_"
                     }
    comn.check_input(params, default_params, critical_params)
        
    labels = params["labels"]      
    nstates = params["nstates"]    
    mopac_exe = params["mopac_exe"]
    mopac_run_params = params["mopac_run_params"]
    mopac_wd = params["mopac_working_directory"]
    mopac_jobid = params["mopac_jobid"]
    mopac_input_prefix = params["mopac_input_prefix"]
    mopac_output_prefix = params["mopac_output_prefix"]
    
    natoms = len(labels)
    ndof = 3 * natoms
        
    obj = tmp()
    obj.ham_adi = CMATRIX(nstates, nstates)    
    obj.nac_adi = CMATRIX(nstates, nstates)    
    obj.hvib_adi = CMATRIX(nstates, nstates)            
    obj.basis_transform = CMATRIX(nstates, nstates) 
    obj.time_overlap_adi = CMATRIX(nstates, nstates)
    obj.d1ham_adi = CMATRIXList();
    obj.dc1_adi = CMATRIXList();          
    for idof in range(ndof):
        obj.d1ham_adi.append( CMATRIX(nstates, nstates) )
        obj.dc1_adi.append( CMATRIX(nstates, nstates) )
  

    Id = Cpp2Py(full_id)
    indx = Id[-1]

    coords = q.col(indx)

    # Create working directory, if doesn't exist
    if not os.path.exists(mopac_wd):
        os.mkdir(mopac_wd)

    # Go into that directory
    os.chdir(mopac_wd)
 
    # Create input
    mopac_input_filename = F"{mopac_input_prefix}{mopac_jobid}"
    make_mopac_input(mopac_input_filename, mopac_run_params, labels, coords)

    # Run the MOPAC job
    mopac_output_filename = F"{mopac_output_prefix}{mopac_jobid}"
    os.system( F"{mopac_exe} {mopac_input_filename} > {mopac_output_filename}")

    # Go back to the original directory
    os.chdir("../")
    
    """    
        # At this point, we should have the "detailed.out" file created, so lets read it        
        E, grad = read_dftb_output(natoms, istate)

        # Now, populate the allocated matrices                
        obj.ham_adi.set(istate, istate, E * (1.0+0.0j) )
        obj.hvib_adi.set(istate, istate, E * (1.0+0.0j) )        
        obj.basis_transform.set(istate, istate, 1.0+0.0j )        
        obj.time_overlap_adi.set(istate, istate, 1.0+0.0j )
        for idof in range(ndof):        
            obj.d1ham_adi[idof].set(istate, istate, grad.get(idof, 0) * (1.0+0.0j) )                
    """
                    
    return obj


def mopac_nbra_workflow(params_):
    pass
    #labels, q = cp2k.read_trajectory_xyz_file("1_ring-pos-1.xyz", 0)
    #run_mopac_adi(q, params_, full_id)



def read_mopac_output(natoms, istate):
    """
    TO DO: This function is yet to be implemented 
    This file reads in the total energy (essentially the reference, ground state energy), 
    the excitation energies, and the corresponding forces on all atoms and return them in 
    a digital format for further processing.
    
    Args:
    
        natoms ( int ): the number of atoms in the systemm
        istate ( int ): the index of the calculation - just for archiving purposes
        
    Returns: 
        double, MATRIX(ndof, 1): the total energy of a given electronic state,         
            and the corresponding gradients
    """
    
    # Check the successful completion of the calculations like this:
    if os.path.isfile('detailed.out'):
        #print("Found detailed.out")
        os.system(F"cp detailed.out detailed_{istate}.out")
    else:
        print("\nCannot find the file detailed.out")
        print(F"Hint: Current working directory is: {os.getcwd()}")
        print("Is this where you expect the file detailed.out to be found?")
        print("Exiting program...\n")
        sys.exit(0)
            
    
    ndof = 3 * natoms
    grad = MATRIX(ndof, 1)
    
    
    f = open("detailed.out")
    output = f.readlines()
    f.close()
    
    E = 0.0
    nlines = len(output)
            
    for i in range( nlines ):

        output_line = output[i].split()

        if len( output_line ) >= 2:
            
            if output_line[0] == "Total" and output_line[1] == "Forces":

                for j in range( natoms ):

                    next_lines = output[i+j+1].split()  

                    grad.set( 3*j+0, 0, -float(next_lines[1]) )
                    grad.set( 3*j+1, 0, -float(next_lines[2]) )
                    grad.set( 3*j+2, 0, -float(next_lines[3]) )

            if output_line[0] == "Excitation" and output_line[1] == "Energy:" :
            
                E += float(output_line[2])  # energy in a.u.
                #print(output_line[2])
                
            if output_line[0] == "Total" and output_line[1] == "energy:" :
                
                E += float(output_line[2])  # energy in a.u.
                #print(output_line[2])
    
    #print(F"Final energy = {E}")
    return E, grad          

