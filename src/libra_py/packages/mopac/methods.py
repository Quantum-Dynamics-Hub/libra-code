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

    mopac_input.write(F"{mopac_run_params}\n\n\n")

    nat = len(labels) # how many atoms 
    for i in range(nat):
        x = coords.get(3*i+0, 0)/units.Angst
        y = coords.get(3*i+0, 1)/units.Angst
        z = coords.get(3*i+0, 2)/units.Angst
        mopac_input.write(F"{labels[i]}   {x}  1  {y} 1   {z}  1\n")
    mopac_input.write("\n")
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



def read_mopac_orbital_info(params_):
    """
    This function reads the MOs, configurations, and CI from the output files
    
    Args:
    
        params_ ( dict ): the dictionary containing key parameters

            * **params_["filename"]** ( string ) : the name of the file to read
        
    Returns: 
        (Es, MOs, E_CI, CI, configs): 
     
            * Es - MATRIX(nmo, nmo): the matrix of the MO energies
            * MOs - MATRIX(nao, nmo): the matrix of MO-LCAO coefficients
            * E_CI - MATRIX(nci, nci): the matrix of CI energies
            * CI - MATRIX(nconf, nci): the matrix of CI coefficients in the basis of spin-adapted configurations
            * configs - (list of lists): information about the configurations, the first list is [-1, -1] - it corresponds 
                 to the ground state configurations, other lists are of the kinds [i, j] which corresponds to i->j excitation
                 The indices i and j are enumerated starting from 0 (Python/C++ convention), which is different from the 
                 raw MOPAC output, where MO indexing starts from 1
    Notes:
            * nao - the number of AOs, it is the same as nmo
            * nmo - the number of MOs
            * nconf  - the number of configurations (spin-adapted Slater determinants)
            * nci - the number of CI states
                              
 
    """
    params = dict(params_)

    critical_params = [ ]
    default_params = { "filename":"output"  }
    comn.check_input(params, default_params, critical_params)

    out_file = params["filename"]

    
    # Check the successful completion of the calculations like this:
    if os.path.isfile(out_file):
        pass
    else:
        print(F"Cannot find the file {out_file}")
        print(F"Hint: Current working directory is: {os.getcwd()}")
        print("Is this where you expect the file detailed.out to be found?")
        print("Exiting program...\n")
        sys.exit(0)
                 
    f = open(out_file)
    output = f.readlines()
    f.close()
    nlines = len(output)

    # First, let's find where the MOs are and count how many of them we have
    ibeg, iend, nmo = 0, nlines-1, 0
    for i in range(nlines):
        line = output[i]

        if line.find("MOLECULAR ORBITALS") != -1:
            ibeg = i
        if line.find("Reference determinate nber") != -1:
            iend = i

        if line.find("ROOT NO.") != -1:
            nmo = int(float(line.split()[-1]))

    if False: # Make True for debugging
        print("nmo = ", nmo)
        print(output[ibeg:iend])

    nao = nmo 

    # Find the line indices that contain "ROOT NO." keyword
    # the last one will be the `iend`
    break_lines = []
    for i in range(ibeg, iend):
        line = output[i]
        if line.find("ROOT NO.") != -1:
            break_lines.append(i)
    break_lines.append(iend)
    nblocks = len(break_lines)

    # Now read energies:
    E = []
    Es = MATRIX(nmo, nmo)
    for j in range(nblocks-1):
        i = break_lines[j]
        tmp = output[i+2].split()
        for e in tmp:
            ener = float(e)
            E.append( ener )
    for i in range( nmo ):
        Es.set(i,i, E[i])

    if False: # Make True for debugging
        print(E)

    # Now read MOs:
    mo_indx = 0
    MOs = MATRIX(nao, nmo)
    for j in range(nblocks-1):
        i = break_lines[j]
        tmp = output[i+2].split()
        ncols = len(tmp)
        ao_indx = 0
        for i in range( break_lines[j]+5, break_lines[j+1]):
            tmp = output[i].split()
            sz = len(tmp)
            if( sz == ncols + 3 ):
                for a in range(ncols):
                    coeff = float(tmp[3+a])
                    MOs.set(ao_indx, mo_indx + a, coeff)
                ao_indx += 1
        mo_indx += ncols

    if False: # Make True for debugging
        MOs.show_matrix("MO.txt")


    # Find the configurations
    nconfig = 0
    configs = []
    for i in range(nlines):
        line = output[i]
        if line.find("The lowest") != -1 and line.find("spin-adapted configurations of multiplicity=  1"):
            tmp = line.split()
            if( len(tmp) > 3): 
                nconfig = int( float(tmp[2]) )
                #========= Found the configurations info, now read the configurations ========
                for iconfig in range(nconfig):
                    tmp = output[i + 4 + iconfig].split()
                    if len(tmp) == 12: 
                        configs.append( [-1, -1])  # this is the ground state determinants
                    if len(tmp) == 13:
                        i_orb = int(float(tmp[10].split(")->(")[0]))
                        j_orb = int(float(tmp[11]))
                        configs.append( [ i_orb-1, j_orb-1])  # orbital indexing from 0


    if False: # Make True for debugging
        print(F"The number of spin-adapted configurations = {nconfig}")
        print(F"configs are = {configs}")


    # Find the line indices that contain the beginning and end of the CI information
    ci_beg, ci_end = [], []
    for i in range(nlines):
        line = output[i]
        if line.find("State") != -1 and line.find("CI coeff") != -1 and line.find("CI percent") != -1:
            ci_beg.append(i)
        if line.find("Total coeff printed") != -1:
            ci_end.append(i)

    nci = len(ci_beg)

    if False: # Make True for debugging
        print(F"The number of CI states = {nci}") 
        for i in range(nci):
            print(F"CI block {i}")
            print(output[ci_beg[i]:ci_end[i]])


    #  Now, read the information about CI states
    E_CI = MATRIX(nci, nci)
    CI = MATRIX(nconfig, nci)
    for i in range(nci):
        for j in range(ci_beg[i], ci_end[i]):
            tmp = output[j].split()
            sz = len(tmp)
            if sz==7:
                ener = float( tmp[2] ) * units.ev2au  # convert to a.u.
                E_CI.set(i,i, ener)
            elif sz==4:
                if tmp[0]=="Config":
                    iconf = int(float(tmp[1])) - 1
                    coeff = float(tmp[2])
                    CI.set(iconf, i, coeff)
    

    return Es, MOs, E_CI, CI, configs




