# *********************************************************************************
# * Copyright (C) 2019-2023  Alexey V. Akimov, Brendan Smith
# *
# * This file is distributed under the terms of the GNU General Public License
# * as published by the Free Software Foundation, either version 3 of
# * the License, or (at your option) any later version.
# * See the file LICENSE in the root directory of this distribution
# * or <http://www.gnu.org/licenses/>.
# *
# *********************************************************************************/
"""
.. module:: DFTB_methods
   :platform: Unix, Windows
   :synopsis: This module implements functions for dealing with the outputs from DFTB+ package

.. moduleauthor::
       Alexey V. Akimov, Thomas A. Niehaus

"""


import os
import sys
import math
import re
import struct
import numpy as np

if sys.platform == "cygwin":
    from cyglibra_core import *
elif sys.platform == "linux" or sys.platform == "linux2":
    from liblibra_core import *
import util.libutil as comn

from libra_py import units
from libra_py import scan
from libra_py import regexlib as rgl

import libra_py.packages.cp2k.methods as CP2K_methods
# import numpy as np


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
    for i in range(0, sz):
        tmp = A[i].split()

        # =========== Look for forces =============
        if len(tmp) == 2:
            if tmp[0] == "Total" and tmp[1] == "Forces":
                for j in range(i + 1, i + 1 + nat):
                    ind = j - i - 1
                    tmp1 = A[j].split()
                    x = float(tmp1[0])
                    y = float(tmp1[1])
                    z = float(tmp1[2])
                    F.set(3 * ind + 0, 0, x)
                    F.set(3 * ind + 1, 0, y)
                    F.set(3 * ind + 2, 0, z)
                    Flst.append([x, y, z])

        # =========== Excitation energy =============
        if len(tmp) > 3:
            if tmp[0] == "Excitation" and tmp[1] == "Energy:":
                E_ex = float(tmp[2])

    return E_ex, F, Flst


def get_dftb_matrices(filename, act_sp1=None, act_sp2=None, correction=0, tol=0.5):
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
        correction (int): 0 - don't do a correction, 1 - do it; this correction is useful for
            reading the time-overlap files from the super-molecule calculations. Because of the how
            integrals are evaluated in DFTB+, sometimes the overlaps for the same orbitals on the same
            position would not be the same as expected from the corresponding calculations using
            orbitals of only one atom. As a result, NaNs or incorrect values are observed in the off-diagonal
            blocks. When the correction is turned on, we'll first check for NaNs in the off-diagonal blocks. If
            NaNs are observed, the value is temporarily set to 1000. Then, the magnitudes of the elements in the
            off-diagonal blocks are compared to the the corresponding elements of the diagonal blocks. If the
            aboslute value of the difference is larger than the `tol` value, the elements in the off-diagonal
            blocks are reset to the values of either 0 or the corresponding value of the in-diagonal blocs. The
            latter is the case for the diagonal elements of the blocks, and zero is for the off-diagonal elements
            of the block. [default: 0]
        tol (float ): the threshold for fixing the elements of the off-diagonal blocks. The off-diagonal block elements
            are only reset to the corresponding values of the on-diagonal blocks (or to zeros) if they exceed them by this
            amount. The larger values of this parameter will favor keeping the off-diagonal block matrix elements to
            what they are in the direct calculations, which may be problematic for very close distances (e.g. identical geometries).
            On the contrary, the smaller values of this parameter will favor replacing the off-diagonal block matrix elements
            with the corresponding on-diagonal block matrix elements or zeroes. This may artificially decrease the NACs and may
            slow down the nonadiabatic transitions. E.g. in the extreme case of `tol = 0.0`, all NACs would be zero, so there will
            be no dynamics. [ default: 0.5 ]

    Returns:
        list of CMATRIX(N, M): X: where N = len(act_sp1) and M = len(act_sp2)
            These are the corresponding property matrices (converted to the complex type)
            for each k-point, such that X[ikpt] is a CMATRIX(N, M) containing the
            overlaps/Hamiltonians for the k-point ```ikpt```.

    Warning:
        So far tested only for a single k-point!

    """

    # Copy the content of the file to a temporary location
    print(F"Reading file {filename}")
    f = open(filename, "r")
    A = f.readlines()
    f.close()
    norbs = int(float(A[1].split()[1]))
    nkpts = int(float(A[1].split()[2]))

    # Determine the dimensions of the output matrices
    if act_sp1 is None:
        nstates1 = norbs
        act_sp1 = list(range(0, nstates1))
    else:
        nstates1 = len(act_sp1)

    if act_sp2 is None:
        nstates2 = norbs
        act_sp2 = list(range(0, nstates2))
    else:
        nstates2 = len(act_sp2)

    # Output variables
    X = []
    norbs2 = int(norbs / 2)

    print(F"norbs = {norbs}, norbs2 = {norbs2}, act_sp1 = {act_sp1}, act_sp2 = {act_sp2}")

    # For each k-point
    for ikpt in range(0, nkpts):

        # Write just the matrix data into a temporary files
        # This procedure is made fault-resistant, to avoid wrong working
        # when some of the matrix elements printed out are nonsense.
        B = A[5 + ikpt * norbs: 5 + (1 + ikpt) * norbs]

        f1 = open(filename + "_tmp", "w")
        x = MATRIX(norbs, norbs)

        for i in range(0, norbs):

            tmp = B[i].split()
            line = ""
            for j in range(0, norbs):

                z = 0.0
                if tmp[j] == "NaN" or tmp[j] == "NAN":
                    z = 1000.0
                    # z = float(tmp[int(j%norbs2)]) #-1000.0
                    sys.exit(0)
                else:
                    try:
                        z = float(tmp[j])
                        x.set(i,j, z)
                    except ValueError:
                        sys.exit(0)
                        #z = 1000.0
                        # z = float(tmp[int(j%norbs2)]) #-1000.0
                        #pass
                    except TypeError:
                        sys.exit(0)
                        #z = 1000.0
                        # z = float(tmp[int(j%norbs2)]) #-1000.0
                        #pass

                """
                if correction == 1:  # if we want to enforce correction of the matrix elements in blocks 01 and 10
                    #           norbs2    norbs
                    #       ______________
                    #       |__00__|__01__|    norbs2
                    #       |__10__|__11__|    norbs
                    #

                    if i < norbs2 and j >= norbs2:  # block 01
                        z_diag = float(tmp[j - norbs2])
                        if abs(z - z_diag) > tol:
                            if j - norbs2 == i:
                                z = z_diag
                            else:
                                z = 0.0
                        else:
                            pass  # keep z to what it is

                    if i >= norbs2 and j < norbs2:  # block 10
                        z_diag = float(tmp[j + norbs2])

                        if abs(z - z_diag) > tol:
                            if i - norbs2 == j:
                                z = z_diag
                            else:
                                z = 0.0
                        else:
                            pass
                """
                #line = line + "%10.8f  " % (z)
                line = F"{line} {z}" 
            line = line + "\n"

            f1.write(line)
        f1.close()

        # Read in the temporary file - get the entire matrix
        #x = MATRIX(norbs, norbs)
        #x.Load_Matrix_From_File(filename + "_tmp")
      

        # Extract the sub-matrix of interest
        x_sub = MATRIX(nstates1, nstates2)
        pop_submatrix(x, x_sub, act_sp1, act_sp2)

        # Add the resulting matrices to the returned result
        X.append(CMATRIX(x_sub))

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
    nat = int(float(A[0].split()[0]))
    at_types = []

    for i in range(0, nat):
        at_typ = A[i + 2].split()[0]   # atom type
        if at_typ not in at_types:
            at_types.append(at_typ)

    # Make up the output file
    line = "%5i  %s \n" % (nat, sys_type)

    for typ in at_types:
        line = line + " %s " % (typ)
    line = line + "\n"

    for i in range(0, nat):
        ln_indx = (nat + 2) * md_iter + (i + 2)
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
    nat = int(float(A[0].split()[0]))
    at_types = []

    for i in range(0, nat):
        at_typ = A[i + 2].split()[0]   # atom type
        if at_typ not in at_types:
            at_types.append(at_typ)

    # Make up the output file
    line = "%5i  %s \n" % (2 * nat, sys_type)

    for typ in at_types:
        line = line + " %s " % (typ)
    line = line + "\n"

    for i in range(0, nat):
        ln_indx = (nat + 2) * md_iter1 + (i + 2)
        tmp = A[ln_indx].split()

        at_indx = i + 1
        at = tmp[0]
        typ_indx = at_types.index(at) + 1
        x = float(tmp[1])
        y = float(tmp[2])
        z = float(tmp[3])
        line = line + " %5i %3i  %15.8f  %15.8f  %15.8f\n" % (at_indx, typ_indx, x, y, z)

    for i in range(0, nat):
        ln_indx = (nat + 2) * md_iter2 + (i + 2)
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


def make_dftb_input(dftb_input_template_filename, istate, timestep=-1):
    """
    This function creates an input file for DFTB+ package from a template file,
    it changes several placeholder lines to ensure the calculations are done
    for required electronic states

    Args:
        dftb_input_template_filename ( strings ): the name of the generic input file (template)
            for DFTB+ calculations

        istate ( int ): the index of the state, for which the calculations will be carried out

        timestep ( int ): the index of a given timestep along a MD trajectory. If not provided (-1),
                          then the defaul name of "tmp.gen" will be used as the coordinate file name.
                          If a time index is given, the name will be "coord-<timestep>.gen".

    Returns:
        None :  just creates the files

    """

    # Read in the template file
    f = open(dftb_input_template_filename, 'r')
    dftb_input_template = f.readlines()
    f.close()

    nlines = len(dftb_input_template)

    # Create the actual output file - just replace few parameters in the template file
    dftb_input = open("dftb_in.hsd", "w")

    for i in range(nlines):

        dftb_input_template_line = dftb_input_template[i].split()

        if len(dftb_input_template_line) > 0:

            if dftb_input_template_line[0] == "<<<" and timestep == -1:
                dftb_input.write('    <<< "tmp.gen"\n ')

            elif dftb_input_template_line[0] == "<<<" and timestep != -1:
                dftb_input.write('    <<< "coord-' + str(timestep) + '.gen"\n ')

            elif dftb_input_template_line[0] == "StateOfInterest":
                dftb_input.write(F"    StateOfInterest    = {istate}\n")

            elif istate == 0 and dftb_input_template_line[0] == "ExcitedStateForces":
                dftb_input.write("    ExcitedStateForces = no\n")

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
        # print("Found detailed.out")
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

    for i in range(nlines):

        output_line = output[i].split()

        if len(output_line) >= 2:

            if output_line[0] == "Total" and output_line[1] == "Forces":

                for j in range(natoms):

                    next_lines = output[i + j + 1].split()

                    grad.set(3 * j + 0, 0, -float(next_lines[1]))
                    grad.set(3 * j + 1, 0, -float(next_lines[2]))
                    grad.set(3 * j + 2, 0, -float(next_lines[3]))

            if output_line[0] == "Excitation" and output_line[1] == "Energy:":

                E += float(output_line[2])  # energy in a.u.
                # print(output_line[2])

            if output_line[0] == "Total" and output_line[1] == "energy:":

                E += float(output_line[2])  # energy in a.u.
                # print(output_line[2])

    # print(F"Final energy = {E}")
    return E, grad


class tmp:
    pass


def run_dftb_adi(q, params_, full_id):
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

    critical_params = ["labels"]
    default_params = {"dftb_input_template_filename": "dftb_input_template.hsd",
                      "nstates": 1,
                      "dftb_exe": "dftb+", "xyz2gen_exe": "xyz2gen"
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
    obj.basis_transform = CMATRIX(nstates, nstates)
    obj.time_overlap_adi = CMATRIX(nstates, nstates)
    obj.d1ham_adi = CMATRIXList()
    obj.dc1_adi = CMATRIXList()
    for idof in range(ndof):
        obj.d1ham_adi.append(CMATRIX(nstates, nstates))
        obj.dc1_adi.append(CMATRIX(nstates, nstates))

    Id = Cpp2Py(full_id)
    indx = Id[-1]

    # Make an xyz file
    # since the DFTB+ expects the coordinates in Angstrom, but Libra
    # goes with atomic units (Bohrs), hence expecting the `q` variable be
    # in Bohrs, we need to convert the coordinates from Bohrs to Angstroms
    coords_str = scan.coords2xyz(labels, q, indx, 1.0 / units.Angst)
    f = open("tmp.xyz", "w")
    f.write(F"{natoms}\nTemporary xyz file\n{coords_str}")
    f.close()

    # Convert xyz to gen format: tmp.xyz -> tmp.gen
    # The temporary working file MUST be called "tmp.gen" since this is
    # what the DFTB+ input will list - see the `make_dftbp_input`
    os.system(F"{xyz2gen_exe} tmp.xyz")

    for istate in range(nstates):

        # Update the input file
        make_dftb_input(dftb_input_template_filename, istate)

        # We have written the dftb+ input file for a certain state in nstates. Now we must compute the
        # state energies and forces.
        os.system(dftb_exe)

        # At this point, we should have the "detailed.out" file created, so lets read it
        E, grad = read_dftb_output(natoms, istate)

        # Now, populate the allocated matrices
        obj.ham_adi.set(istate, istate, E * (1.0 + 0.0j))
        obj.hvib_adi.set(istate, istate, E * (1.0 + 0.0j))
        obj.basis_transform.set(istate, istate, 1.0 + 0.0j)
        obj.time_overlap_adi.set(istate, istate, 1.0 + 0.0j)
        for idof in range(ndof):
            obj.d1ham_adi[idof].set(istate, istate, grad.get(idof, 0) * (1.0 + 0.0j))

    return obj


def dftb_distribute(istep, fstep, nsteps_this_job, trajectory_xyz_file, dftb_input, waveplot_input, curr_job_number):
    """
    Distributes dftb jobs for trivial parallelization

    Make sure that your dftb input file has absolute paths to the following input parameters:

        SlaterKosterFiles = Type2FileNames {
          Prefix = "/panasas/scratch/grp-alexeyak/brendan/dftbp_development/"
          Separator = "-"
          Suffix = ".skf"
        }

    Args:

        istep (integer): The initial time step in the trajectory xyz file.

        fstep (integer): The final time step in the trajectory xyz file.

        nsteps_this_job (integer): The number of steps for this job.

        trajectory_xyz_file (string): The full path to trajectory xyz file.

        dftb_input (string): the dftb_input_template.hsd file

        waveplot_input (string): the input file for the waveplot subpackage of dftbplus for generating cubefiles

        curr_job_number (integer): The current job number.

    Returns:

        None

    """

    # Now we need to distribute the jobs into job batches
    # First, make a working directory where the calculations will take place
    os.chdir("wd")

    nsteps = fstep - istep + 1
    njobs = int(nsteps / nsteps_this_job)

    # Initialize the curr_step to istep
    curr_step = istep

    # Make the job directory folders
    os.system("mkdir job" + str(curr_job_number) + "")
    os.chdir("job" + str(curr_job_number) + "")
    # Copy the trajectory file and input template there
    os.system("cp ../../" + trajectory_xyz_file + " .")
    os.system("cp ../../" + dftb_input + " .")
    os.system("cp ../../" + waveplot_input + " .")

    # Now, we need to edit the submit file
    # Now, in jobs folder njob, we should do only a certain number of steps
    for step in range(nsteps_this_job):

        # extract the curr_step xyz coordianates from the trajectory file and write it to another xyz file
        _, _ = CP2K_methods.read_trajectory_xyz_file(trajectory_xyz_file, curr_step)
        curr_step += 1

    # Go back to the main directory
    os.chdir("../../")


def get_dftb_ks_energies(_params):
    """
    This function reads the band.out file generated from a dftb+ computations. The band.out file
    stores the ks energies and their occupation.

    Args:
        params (dict): parameters controlling this calculation. Can contain:

        * **_params["logfile_name"]** (string): the file containing the output of the MD simulations [ default: "band.out"]
        * **_params["min_band"]** (int): index of the minimal orbital/band to include in the analysis, starting from 1 [ default: 1 ]
        * **_params["max_band"]** (int): index of the maximal orbital/band to include in the analysis, starting from 1 [ default: 2 ]
        * **_params["time"]** (int): Time step at which the properties will be read, for molecular dynamics it will read the energies
                 of time step 'time', but for static calculations the time is set to 0 [ default: 0]

    Returns:

        E (1D numpy array): The vector consisting of the KS energies from min_band to max_band.

    """

    params = dict(_params)

    # Critical parameters
    critical_params = []
    # Default parameters
    default_params = {"logfile_name": "band.out", "min_band": 1, "max_band": 2, "time": 0}
    # Check input
    comn.check_input(params, default_params, critical_params)

    dftb_outfile_name = params["logfile_name"]

    # Time step, for molecular dynamics it will read the energies
    # of time step 'time', but for static calculations the time is set to 0
    time = params["time"]

    # The minimum state number
    min_band = params["min_band"]  # ks_orbital_indicies[0]
    # The maximum state number
    max_band = params["max_band"]  # ks_orbital_indicies[-1]

    # First openning the file and reading its lines
    f = open(dftb_outfile_name, 'r')
    lines = f.readlines()
    f.close()

    # The lines containing the energies of the occupied states
    occ_energies = []
    # The lines containing the energies of the unoccupied states
    unocc_energies = []

    # For the dftb output file band.out, start from line 1
    for i in range(1, len(lines)):

        if i >= min_band and i <= max_band:

            b = lines[i].strip().split()
            if float(b[2]) > 0.0:
                occ_energies.append(float(b[1]) * units.ev2Ha)
            else:
                unocc_energies.append(float(b[1]) * units.ev2Ha)

    # Turn them into numpy arrays
    occ_energies = np.array(occ_energies)
    unocc_energies = np.array(unocc_energies)

    # Concatenate the occupied and unoccpied energies so we can choose from min_band to max_band
    ks_energies = np.concatenate((occ_energies, unocc_energies))

    # print("ks_energies", ks_energies)

    return ks_energies


def cube_generator_dftbplus(project_name, time_step, min_band, max_band, waveplot_exe, isUKS):
    """
    This function generates the cube files by first forming the 'fchk' file from 'chk' file.
    Then it will generate the cube files from min_band to max_band the same as CP2K naming.
    This will helps us to use the read_cube and integrate_cube functions easier.

    Args:

        project_name (string): The project_name.

        time_step (integer): The time step.

        min_band (integer): The minimum state number.

        max_band (integer): The maximum state number.

        waveplot_exe (str): Location of the executable for the waveplot program of the dftb+ software package

        isUKS (integer): The unrestricted spin calculation flag. 1 is for spin unrestricted calculations.
                         Other numbers are for spin restricted calculations

    Returns:

        None

    """

    # Make the cubefiles. For dftb+, this simply means executing the waveplot program
    # The min and max band are already defind in the waveplot template. So, we just use
    # them here for renaming the cubefiles
    os.system(waveplot_exe)

    # For now, only spin-restricted
    for state in range(min_band, max_band + 1):
        # Use cp2k names because the rest of the code expects this format
        state_name = CP2K_methods.state_num_cp2k(state)
        cube_name = '%s-%d-WFN_%s_1-1_0.cube' % (project_name, time_step, state_name)
        print('Renaming cube for state %d' % state)
        # Now, rename the cubefile from what waveplots calls it to the standard format
        os.system("mv *" + str(state) + "-real.cube " + cube_name)


def read_dftbplus_TRA_file(params):
    """
    This function reads TRA.dat files generated from TD-DFTB calculations using DFTB+ and returns
    the TD-DFTB excitation energies and CI-like coefficients

    Args:
        params ( dictionary ): parameters dictionary

            logfile_name ( string ): this is actualyl the TRA.dat file, but we keep it as logfile_name for ease, as many other similar
                                     type functions (for cp2k, gaussian) rely on this format internally.
            number_of_states ( int ): how many ci states to consider
            tolerance ( float ): the tolerance for accepting SDs are members of the CI wavefunctions
            isUKS ( Boolean ): flag for whether or not a spin-polarized Kohn-Sham basis is being used. TRUE means
                               that a spin-polarized Kohn-Sham basis is being used.
    Returns:
        excitation_energies ( list ): The excitation energies in the log file.
        ci_basis ( list ): The CI-basis configuration.
        ci_coefficients ( list ): The coefficients of the CI-states.
        spin_components (list): The spin component of the excited states.

    """

    # Critical parameters
    critical_params = ["logfile_name", "number_of_states"]
    # Default parameters
    default_params = {"tolerance": 0.001, "isUKS": 0}
    # Check input
    comn.check_input(params, default_params, critical_params)

    logfile_name = params["logfile_name"]
    number_of_states = int(params["number_of_states"])
    tolerance = float(params["tolerance"])
    isUKS = int(params["isUKS"])

    f = open(logfile_name, 'r')
    lines = f.readlines()
    f.close()

    # Make lists for excitation energies and the lines where the keyword "Energy" is
    excitation_energies = []
    energy_lines = []

    for i in range(0, len(lines)):
        tmp_line = lines[i].split()
        if 'Energy' in tmp_line:
            # print("Energy")
            # When found the line in which contains 'Energy'
            excitation_energies.append(float(tmp_line[2]))
            energy_lines.append(i)
        elif len(energy_lines) == number_of_states:
            break

    # print( energy_lines )
    # sys.exit(0)

    # Obtain CI-like coefficients
    # Spin-unpolarized only as of 11/6/2020
    # Start from 4 lines after finding the line contaning 'Energy'. This is how it is in DFTB v. 19.1
    nlines_to_skip = 4

    ci_basis = []
    ci_coefficients = []
    spin_components = []
    for i in energy_lines:

        tmp_spin = []
        tmp_state = []
        tmp_state_coefficients = []
        for j in range(i + nlines_to_skip, len(lines)):
            tmp_line = lines[j].split()
            if len(tmp_line) == 0:
                break
            else:
                ci_coefficient = float(tmp_line[3])
                if ci_coefficient**2 > tolerance:
                    tmp_spin.append("alp")
                    tmp_state.append([int(tmp_line[0]), int(tmp_line[2])])
                    tmp_state_coefficients.append(ci_coefficient)

        # Append the CI-basis and and their coefficients for
        # this state into the ci_basis and ci_coefficients lists
        ci_basis.append(tmp_state)
        ci_coefficients.append(tmp_state_coefficients)
        spin_components.append(tmp_spin)

    return excitation_energies[0:number_of_states], ci_basis, ci_coefficients, spin_components


def dftb_traj2xyz_traj(in_filename, out_filename):
    """
    This is an auxiliary function to convert the DFTB+-generated extended xyz trajectory file
    to a more coomon xyz trajectory file format

    Args:
        in_filename ( string ):  the name of the input file
        out_filename ( string ):  the name of the putput file

    Example:
        dftb_traj2xyz_traj("md.xyz", "md_reduced.xyz")

    """

    f = open(in_filename, "r")
    A = f.readlines()
    f.close()

    f = open(out_filename, "w")
    for line in A:
        tmp = line.split()
        res = line

        if len(tmp) == 8:
            res = F"  {tmp[0]}  {tmp[1]}  {tmp[2]}  {tmp[3]}\n"

        f.write(res)
    f.close()

def read_spx_mappings(filename):
    """
    Parses SPX.DAT and returns:
    1. rpa_map: [ini, fin, spin] -> rpa_idx
    2. transition_lookup: [rpa_idx] -> (ini, fin, spin)
    3. max_spin: 0 for unpolarized/Up-only, 1 if Down exists
    Indices follow Python usage and start with 0

    --- Usage ---
    rpa_map, trans_lookup, max_s = read_spx_mappings("SPX.DAT")
    print(f"Max Spin Index: {max_s}") # 0 for Restricted, 1 for Unrestricted
    """
    temp_data = []
    max_orb = 0
    max_rpa_idx = 0
    actual_max_spin = 0
    
    # Mapping for spin characters
    spin_map = {'U': 0, 'D': 1}

    with open(filename, 'r') as f:
        for line in f:
            if not line.strip() or '#' in line or '===' in line:
                continue
            
            parts = line.split()
            if len(parts) < 6:
                continue

            try:
                idx = int(parts[0]) - 1  # 0-indexed rpa_idx
                m = int(parts[3]) - 1    # 0-indexed ini
                n = int(parts[5]) - 1    # 0-indexed fin
                
                # Check for spin character
                s_val = 0
                if len(parts) > 6:
                    s_char = parts[6].upper()
                    s_val = spin_map.get(s_char, 0)
                    if s_val > actual_max_spin:
                        actual_max_spin = s_val
                
                temp_data.append((m, n, s_val, idx))
                
                # Track max values for array sizing
                max_orb = max(max_orb, m, n)
                max_rpa_idx = max(max_rpa_idx, idx)
                
            except (ValueError, IndexError):
                continue

    # 1. Forward Map: [ini, fin, spin] -> rpa_idx
    # Use actual_max_spin + 1 to determine depth (1 or 2)
    rpa_map = np.full((max_orb + 1, max_orb + 1, actual_max_spin + 1), -1, dtype=int)

    # 2. Reverse Map: [rpa_idx] -> [ini, fin, spin]
    transition_lookup = np.zeros((max_rpa_idx + 1, 3), dtype=int)

    # Fill both structures
    for m, n, s, idx in temp_data:
        rpa_map[m, n, s] = idx
        transition_lookup[idx] = [m, n, s]
        
    return rpa_map, transition_lookup, actual_max_spin


def parse_tagged_file(file_path):
    """
    --- Usage ---
    Load the file into a data dictionary
    results = parse_tagged_file("autotes.tag")

    Accessing a 0D (Scalar) value by tag name
    energy = results['mermin_energy']  # Returns an int or float

    Accessing a 1D or 2D (Array) value
    coords = results['end_coords']  # Returns a NumPy array

    """

    data_dict = {}
    
    with open(file_path, 'r') as f:
        content = f.read()

    # Split by lines and filter empty ones
    lines = [line.strip() for line in content.split('\n') if line.strip()]
    
    i = 0
    while i < len(lines):
        line = lines[i]
        
        # Identify header: contains tag and exactly 3 colons
        if line.count(':') >= 3:
            parts = [p.strip() for p in line.split(':')]
            tag = parts[0]
            dtype_str = parts[1]
            dim = int(parts[2])
            shape_str = parts[3]
            
            # 1. Determine Shape and Total Elements
            if dim == 0:
                shape = (1,)
                total_elements = 1
            else:
                # Convert "3,23" -> (3, 23)
                shape = tuple(int(s) for s in shape_str.split(',') if s)
                total_elements = np.prod(shape)
            
            # 2. Collect data values
            data_values = []
            i += 1
            while len(data_values) < total_elements and i < len(lines):
                if lines[i].count(':') >= 3: # Safety break if header is missing data
                    break
                # Replace Fortran 'D' with 'E' for float conversion
                clean_line = lines[i].replace('D', 'E')
                data_values.extend(clean_line.split())
                i += 1
            
            # 3. Convert and Return as Scalar or Array
            if dim == 0:
                # Return as a simple Python number
                val = float(data_values[0])
                data_dict[tag] = int(val) if dtype_str == 'integer' else val
            else:
                # Create NumPy array
                array = np.array(data_values, dtype=float)
                # Use Fortran order 'F' because Fortran fills columns first
                data_dict[tag] = array.reshape(shape, order='F')
        else:
            i += 1
            
    return data_dict


def read_mo_matrix(filename, ndim, spin_polarized=False):
    """
    Reads MO eigenvectors and always returns a 3D array of shape (spin, ndim, ndim).
    Spin 1: [4B Integer][ndim*ndim*8B Matrix]
    Spin 2: [ndim*ndim*8B Matrix]
    """
    mo_data = []

    with open(filename, 'rb') as f:
        # --- RECORD 1 (Spin Up / Total) ---
        
        # 1. Read the extra integer header (4 bytes)
        header_int = np.fromfile(f, dtype='i4', count=1)[0]
        
        # 2. Read the matrix
        matrix_up = np.fromfile(f, dtype='f8', count=ndim*ndim)
        mo_data.append(matrix_up.reshape((ndim, ndim), order='F'))
        
        # --- RECORD 2 (Spin Down) ---
        if spin_polarized:
            
            # 2. Read the matrix directly (no header integer here as specified)
            matrix_down = np.fromfile(f, dtype='f8', count=ndim*ndim)
            mo_data.append(matrix_down.reshape((ndim, ndim), order='F'))
            
    # Returns a 3D array: (1, ndim, ndim) or (2, ndim, ndim)
    return np.array(mo_data)

def read_xplusy_binary(filename):
    """
    Reads XplusY.BIN (Stream Unformatted)
    Matches: dp = real64 (8 bytes), ii = int32 (4 bytes), sign = char (1 byte)
    """
    results = []
    
    with open(filename, 'rb') as f:
        # Step 1: Read Header
        # nmat (i4), nExc (i4)
        header_data = f.read(8)
        if len(header_data) < 8: return None
        nmat, nexc = struct.unpack('ii', header_data)
        
        for _ in range(nexc):
            # Step 3: Root Header
            # ii (i4), sign (c1)
            ii = struct.unpack('i', f.read(4))[0]
            sign = f.read(1).decode('ascii')
            
            # Handling the Energy Slot
            # The Fortran code writes sqrt(eval) [8 bytes] OR '-' [1 byte]
            # This is tricky for stream binary. We check the next byte.
            peek = f.read(1)
            if peek == b'-':
                ene = None
                # No 8-byte float followed, move on
            else:
                # Put the peeked byte back and read as float64
                f.seek(-1, os.SEEK_CUR) 
                ene = struct.unpack('d', f.read(8))[0]
            
            # Step 4: Read xpy vector (nmat * 8 bytes)
            # Only read vector if energy was positive (based on your snippet logic)
            if ene is not None:
                vector = np.fromfile(f, dtype='f8', count=nmat)
            else:
                vector = np.array([])

            results.append({
                'root': ii,
                'sign': sign,
                'energy': ene,
                'vector': vector
            })
            
    return nmat, nexc, results

def read_xplusy_ascii(filename):
    """
    Reads XplusY.DAT (Formatted Text)
    """
    results = []
    with open(filename, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]
        
        # Header
        nmat, nexc = map(int, lines[0].split())
        
        idx = 1
        while idx < len(lines):
            # Root Header: ii, sign, sqrt(eval)
            parts = lines[idx].split()
            ii = int(parts[0])
            sign = parts[1]
            
            val_str = parts[2]
            if val_str == '-':
                ene = None
                vector = np.array([])
                idx += 1
            else:
                ene = float(val_str)
                # Vector follows on subsequent lines, 6 elements per line
                num_vec_lines = (nmat + 5) // 6
                vec_raw = " ".join(lines[idx+1 : idx+1+num_vec_lines])
                vector = np.fromstring(vec_raw, sep=' ')
                idx += 1 + num_vec_lines
            
            results.append({
                'root': ii,
                'sign': sign,
                'energy': ene,
                'vector': vector
            })
            
    return nmat, nexc, results

def read_overlap_matrix(filename, n_orb):
    """
    Reads the ASCII overlap matrix using np.loadtxt.
    n_orb: dimension of the square matrix.
    """
    # np.loadtxt handles the ES24.15 format automatically.
    # We skip 6 lines and read the entire flat block.
    flat_data = np.loadtxt(filename, skiprows=5)
    
    # Reshape: 'F' order is mandatory because Fortran flattens matrices column-first.
    # Shape becomes (n_orb, n_orb)
    matrix = flat_data.reshape((n_orb, n_orb), order='F')
    
    return matrix

def check_unity_deviation(A):
    """
    Calculates the deviation of matrix A from the Identity matrix.
    """
    n = A.shape[0]
    identity = np.eye(n)
    
    # Calculate the residual matrix
    residual = A - identity
    
    # Maximum absolute deviation
    max_dev = np.max(np.abs(residual))
    
    print(f"--- Matrix Deviation Analysis ({n}x{n}) ---")
    print(f"Max Absolute Deviation: {max_dev:.2e}")

def read_nacv(filename):
    """
    Parses NACV.DAT into a 4D array: [state_i, state_j, dim, atomindex]
    - state_i, state_j: Use raw values from file (1-based if file is 1-based)
    - dim, atomindex: 0-based

    # --- Example of how to access the data ---
    # data = read_nacv("NACV.DAT")
    # x_comp_atom_0 = data[1, 2, 0, 0] # State 1, State 2, X-dim, 1st Atom

    """
    couplings_dict = {}
    max_state = 0
    num_atoms = 0

    with open(filename, 'r') as f:
        lines = f.readlines()

    line_idx = 0
    while line_idx < len(lines):
        line = lines[line_idx].strip()
        if not line:
            line_idx += 1
            continue
        
        parts = line.split()
        # Header line: state_i state_j
        if len(parts) == 2:
            si, sj = int(parts[0]), int(parts[1])
            max_state = max(max_state, si, sj)
            
            coord_block = []
            line_idx += 1
            
            # Read vector components
            while line_idx < len(lines):
                sub_parts = lines[line_idx].split()
                if len(sub_parts) == 2: # Next state pair reached
                    break
                if len(sub_parts) == 3: # X, Y, Z
                    coord_block.append([float(x) for x in sub_parts])
                line_idx += 1
            
            couplings_dict[(si, sj)] = np.array(coord_block)
            num_atoms = len(coord_block)
        else:
            line_idx += 1

    # Array shape: [state_i, state_j, dim, atomindex]
    # We use max_state + 1 so that index 'n' corresponds to state 'n'
    shape = (max_state + 1, max_state + 1, 3, num_atoms)
    nacv = np.zeros(shape)

    for (si, sj), vectors in couplings_dict.items():
        # vectors.T converts (Atoms, 3) -> (3, Atoms)
        # This makes Dim 0-indexed (0=X, 1=Y, 2=Z) 
        # and AtomIndex 0-indexed (0 to N-1)
        nacv[si, sj, :, :] = vectors.T
        
        # Standard anti-symmetry for NACVs
        nacv[sj, si, :, :] = -vectors.T

    return nacv










