#*********************************************************************************
#* Copyright (C) 2019-2020  Alexey V. Akimov, Brendan Smith
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 3 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
"""
.. module:: DFTB_methods
   :platform: Unix, Windows
   :synopsis: This module implements functions for dealing with the outputs from DFTB+ package

.. moduleauthor:: 
       Alexey V. Akimov 
  
"""


import os
import sys
import math
import re

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from . import units
from . import regexlib as rgl


def get_energy_forces(filename, nat):
    """Get forces from the input file 

    Args:
        filename ( string ): the name of the file to read, usually 
            this is a the file "detailed.out" produced with the input
            containing option "CalculateForces = Yes "

        nat ( int ): the number of atoms in the system

    Returns:
        tuple: ( E_ex, F, Flst ), where:

            * E_ex ( double ): excitation energy for the given state [ units: a.u. ]
            * F ( MATRIX(ndof, 1) ): the forces acting on all atoms [ units: a.u. of force ]
            * Flst ( list of [x, y, z] ): the forces acting on all atoms [ units: a.u. of force ]

    Warning: 
        it is likely not gonna work for files other than "detailed.out"
    
    """

    # Read the file
    f = open(filename, "r")
    A = f.readlines()
    f.close()
    

    # Returned properties initialized
    E_ex = 0.0
    F = MATRIX(3 * nat, 1)
    Flst = []
    
    # Lets look for what we need
    sz = len(A)
    for i in range(0,sz):    
        tmp = A[i].split()
        
        #=========== Look for forces =============
        if len(tmp)==2:
            if tmp[0]=="Total" and tmp[1]=="Forces":
                for j in range(i+1, i+1+nat):
                    ind = j - i - 1 
                    tmp1 = A[j].split()
                    x = float(tmp1[0])
                    y = float(tmp1[1])
                    z = float(tmp1[2])
                    F.set(3*ind+0, 0, x); F.set(3*ind+1, 0, y); F.set(3*ind+2, 0, z)
                    Flst.append([x, y, z])
                    
        #=========== Excitation energy =============
        if len(tmp)>3:
            if tmp[0]=="Excitation" and tmp[1]=="Energy:":
                E_ex = float(tmp[2])
        
    return E_ex, F, Flst



def get_dftb_matrices(filename, act_sp1=None, act_sp2=None):
    """Get the matrices printed out by the DFTB+

    Args: 
        filename ( string ): the name of the file to read, usually 
            these are any of the files: "hamsqr1.dat", "hamsqr2.dat", etc. or 
            "oversqrt.dat". Produced with the input containing option "WriteHS = Yes"
        act_sp1 ( list of ints or None): the row active space to extract from the original files
            Indices here start from 0. If set to None - the number of AOs will be
            determined automatically from the file. [default: None]
        act_sp2 ( list of ints or None): the cols active space to extract from the original files
            Indices here start from 0. If set to None - the number of AOs will be
            determined automatically from the file. [default: None]
    
    Returns: 
        list of CMATRIX(N, M): X: where N = len(act_sp1) and M = len(act_sp2) 
            These are the corresponding property matrices (converted to the complex type)
            for each k-point, such that X[ikpt] is a CMATRIX(N, M) containing the 
            overlaps/Hamiltonians for the k-point ```ikpt```.

    Warning:
        So far tested only for a single k-point!
    
    """

    
    # Copy the content of the file to a temporary location
    f = open(filename, "r")
    A = f.readlines()
    f.close()    
    norbs = int(float(A[1].split()[1]))
    nkpts = int(float(A[1].split()[2]))


    # Determine the dimensions of the output matrices
    if act_sp1==None:
        nstates1 = norbs
        act_sp1 = list(range(0, nstates1))
    else:
        nstates1 = len(act_sp1)

    if act_sp2==None:
        nstates2 = norbs
        act_sp2 = list(range(0, nstates2))
    else:
        nstates2 = len(act_sp2)

           
    # Output variables    
    X = []
    
    # For each k-point
    for ikpt in range(0,nkpts):

        # Write just the matrix data into a temporary files
        # This procedure is made fault-resistant, to avoid wrong working 
        # when some of the matrix elements printed out are nonsense.
        B = A[5+ikpt*norbs : 5+(1+ikpt)*norbs]
        
        f1 = open(filename+"_tmp", "w")  
        for i in range(0,norbs):
            
            tmp = B[i].split()
            line = ""            
            for j in range(0,norbs): 
                z = 0.0
                if tmp[j]=="NaN":
                    z = 0.0
                else:
                    try:
                        z = float(tmp[j])
                    except ValueError:
                        z = 0.0
                        pass
                    except TypeError:
                        z = 0.0
                        pass                    
                line = line + "%10.8f  " % (z)
            line = line + "\n"

            f1.write(line)
        f1.close()
                 
        # Read in the temporary file - get the entire matrix 
        x = MATRIX(norbs, norbs);  
        x.Load_Matrix_From_File(filename+"_tmp")        
        
        # Extract the sub-matrix of interest
        x_sub = MATRIX(nstates1, nstates2)
        pop_submatrix(x, x_sub, act_sp1, act_sp2)

        # Add the resulting matrices to the returned result
        X.append( CMATRIX(x_sub) )    
        
    return X


def xyz_traj2gen_sp(infile, outfile, md_iter, sys_type):
    """
 
    This file converts the xyz trajectory file in the format produced
    by DFTB+ to a gen file containing the `md_iter`-th step geometry


    Args:
        infile ( string ): the name of the input xyz trajectory file
        outfile ( string ): the name of a output gen single point file
        md_iter ( int ): index of the timeframe to extract
        sys_type ( string ): "C" or "P" - the system type - cluster or periodic
    
    Returns:
        none: but creates a file
        
    """

    f = open(infile, "r")
    A = f.readlines()
    f.close()
    sz = len(A)

    # Determine the key info
    nat = int( float(A[0].split()[0] ) )
    at_types = []

    for i in range(0,nat):
        at_typ = A[i+2].split()[0]   # atom type
        if at_typ not in at_types:
            at_types.append(at_typ)


    # Make up the output file
    line = "%5i  %s \n" % (nat, sys_type)

    for typ in at_types:
        line = line + " %s " % (typ)
    line = line + "\n"

    for i in range(0,nat):
        ln_indx = (nat+2)*md_iter + (i + 2)
        tmp = A[ln_indx].split()

        at_indx = i + 1
        at = tmp[0]
        typ_indx = at_types.index(at) + 1
        x = float(tmp[1])
        y = float(tmp[2])
        z = float(tmp[3])
        line = line + " %5i %3i  %15.8f  %15.8f  %15.8f\n" % (at_indx, typ_indx, x, y, z)


    # Write the file
    f1 = open(outfile, "w")
    f1.write(line)
    f1.close()

    
         
def xyz_traj2gen_ovlp(infile, outfile, md_iter1, md_iter2, sys_type):
    """
 
    This file converts the xyz trajectory file in the format produced
    by DFTB+ to a gen file containing the superimposed `md_iter1`-th 
    and `md_iter2`-th steps geometries

    Args:
        infile ( string ): the name of the input xyz trajectory file
        outfile ( string ): the name of a output gen single point file
        md_iter1 ( int ): index of the first timeframe to extract
        md_iter2 ( int ): index of the second timeframe to extract
        sys_type ( string ): "C" or "P" - the system type - cluster or periodic
    
    Returns:
        none: but creates a file
        
    """

    f = open(infile, "r")
    A = f.readlines()
    f.close()
    sz = len(A)

    # Determine the key info
    nat = int( float(A[0].split()[0] ) )
    at_types = []

    for i in range(0,nat):
        at_typ = A[i+2].split()[0]   # atom type
        if at_typ not in at_types:
            at_types.append(at_typ)


    # Make up the output file
    line = "%5i  %s \n" % (2*nat, sys_type)

    for typ in at_types:
        line = line + " %s " % (typ)
    line = line + "\n"


    for i in range(0,nat):
        ln_indx = (nat+2)*md_iter1 + (i + 2)
        tmp = A[ln_indx].split()

        at_indx = i + 1
        at = tmp[0]
        typ_indx = at_types.index(at) + 1
        x = float(tmp[1])
        y = float(tmp[2])
        z = float(tmp[3])
        line = line + " %5i %3i  %15.8f  %15.8f  %15.8f\n" % (at_indx, typ_indx, x, y, z)

    for i in range(0,nat):
        ln_indx = (nat+2)*md_iter2 + (i + 2)
        tmp = A[ln_indx].split()

        at_indx = nat + i + 1
        at = tmp[0]
        typ_indx = at_types.index(at) + 1
        x = float(tmp[1])
        y = float(tmp[2])
        z = float(tmp[3])
        line = line + " %5i %3i  %15.8f  %15.8f  %15.8f\n" % (at_indx, typ_indx, x, y, z)


    # Write the file
    f1 = open(outfile, "w")
    f1.write(line)
    f1.close()






def make_dftp_input(dftb_input_template_filename, istate):
    """
    This function creates an input file for DFTB+ package from a template file,
    it changes several placeholder lines to ensure the calculations are done
    for required electronic states
    
    Args: 
        dftb_input_template_filename ( strings ): the name of the generic input file (template) 
            for DFTB+ calculations
            
        istate ( int ): the index of the state, for which the calculations will be carried out
        
    Returns:
        None :  just creates the files
        
    """
    
    # Read in the template file
    f = open(dftb_input_template_filename, 'r')
    dftb_input_template = f.readlines()
    f.close() 

    
    nlines = len(dftb_input_template)
    
    # Create the actual output file - just replace few parameters in the template file
    dftb_input = open("dftb_in.hsd","w");     
    
    for i in range( nlines ):                
        
        dftb_input_template_line = dftb_input_template[i].split()
        
        
        if len(dftb_input_template_line) > 0:
                    
            if dftb_input_template_line[0] == "<<<":
                dftb_input.write( '    <<< "tmp.gen"\n ')

            elif dftb_input_template_line[0] == "StateOfInterest":
                dftb_input.write( F"    StateOfInterest    = {istate}\n" )

            elif istate == 0 and dftb_input_template_line[0] == "ExcitedStateForces":
                dftb_input.write( "    ExcitedStateForces = no\n" )

            else:
                dftb_input.write(dftb_input_template[i])
                
        else:
            dftb_input.write(dftb_input_template[i])

                    
    dftb_input.close()



def read_dftb_output(natoms, istate):
    """
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



def run_dftbp_adi(q, params_, full_id):
    """
   
    This function executes the DFTB+ quantum chemistry calculations and 
    returns the key properties needed for dynamical calculations.

    Args: 
        q ( MATRIX(ndof, ntraj) ): coordinates of the particle [ units: Bohr ]
        params ( dictionary ): model parameters

            * **params["labels"]** ( list of strings ): the labels of atomic symbolc - for all atoms,
                and in a order that is consistent with the coordinates (in triples) stored in `q`.
                The number of this labels is `natoms`, such that `ndof` = 3 * `natoms`. [ Required ]                
            * **params["nstates"]** ( int ): the total number of electronic states 
                in this model [ default: 1 - just the ground state]                     
            * **params["dftb_input_template_filename"]** ( string ):  the name of the input file template
                [ default: "dftb_input_template.hsd" ]                
            * **params["dftbp_exe"]** ( string ):  the full path to the DFTB+ executable 
                [ defaut: "dftb+" ]
            * **params["xyz2gen_exe"]** ( string ):  the full path to the xyz2gen executable 
                (part of the DFTB+ package) [ defaut: "xyz2gen" ]                
                Note: sometimes, especially if you are using conda-installed Python, you may need to
                edit the "xyz2gen" file in your DFTB+ installation to change the topmost line to 
                point to the correct python executable (e.g. #! /home/alexey/miniconda2/envs/py37/bin/python)            
                            
    Returns:       
        PyObject: obj, with the members:

            * obj.ham_adi ( CMATRIX(nstates,nstates) ): adiabatic Hamiltonian             
            * obj.d1ham_adi ( list of ndof CMATRIX(nstates, nstates) objects ): 
                derivatives of the adiabatic Hamiltonian w.r.t. the nuclear coordinate            
 
    """

    params = dict(params_)
    
    critical_params = [ "labels" ] 
    default_params = { "dftb_input_template_filename":"dftb_input_template.hsd",     
                       "nstates":1,
                       "dftb_exe":"dftb+",  "xyz2gen_exe":"xyz2gen"
                     }
    comn.check_input(params, default_params, critical_params)
        
    labels = params["labels"]      
    dftb_input_template_filename = params["dftb_input_template_filename"]
    nstates = params["nstates"]    
    dftb_exe = params["dftb_exe"]
    xyz2gen_exe = params["xyz2gen_exe"]
    
    natoms = len(labels)
    ndof = 3 * natoms
    
    
    obj = tmp()
    obj.ham_adi = CMATRIX(nstates, nstates)    
    obj.nac_adi = CMATRIX(nstates, nstates)    
    obj.hvib_adi = CMATRIX(nstates, nstates)            
    obj.d1ham_adi = CMATRIXList();            
    for idof in range(ndof):        
        obj.d1ham_adi.append( CMATRIX(nstates, nstates) )
  

    Id = Cpp2Py(full_id)
    indx = Id[-1]

    
    # Make an xyz file
    # since the DFTB+ expects the coordinates in Angstrom, but Libra 
    # goes with atomic units (Bohrs), hence expecting the `q` variable be 
    # in Bohrs, we need to convert the coordinates from Bohrs to Angstroms
    coords_str = scan.coords2xyz(labels, q, indx, 1.0/units.Angst)
    f = open("tmp.xyz", "w"); 
    f.write( F"{natoms}\nTemporary xyz file\n{coords_str}" )
    f.close()
        
    # Convert xyz to gen format: tmp.xyz -> tmp.gen
    # The temporary working file MUST be called "tmp.gen" since this is
    # what the DFTB+ input will list - see the `make_dftbp_input`
    os.system(F"{xyz2gen_exe} tmp.xyz")

        
    for istate in range(nstates):
        
        # Update the input file
        make_dftp_input(dftb_input_template_filename, istate)
        
        # We have written the dftb+ input file for a certain state in nstates. Now we must compute the 
        # state energies and forces. 
        os.system( dftb_exe )
        
        # At this point, we should have the "detailed.out" file created, so lets read it        
        E, grad = read_dftb_output(natoms, istate)

        # Now, populate the allocated matrices                
        obj.ham_adi.set(istate, istate, E * (1.0+0.0j) )
        obj.hvib_adi.set(istate, istate, E * (1.0+0.0j) )        
        for idof in range(ndof):        
            obj.d1ham_adi[idof].set(istate, istate, grad.get(idof, 0) * (1.0+0.0j) )                
    
                    
    return obj

