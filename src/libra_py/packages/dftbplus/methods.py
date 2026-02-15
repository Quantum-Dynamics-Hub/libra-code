# *********************************************************************************
# * Copyright (C) 2026  Alexey V. Akimov, Brendan Smith, Thomas A. Niehaus
# * Copyright (C) 2019-2025  Alexey V. Akimov, Brendan Smith
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
       Alexey V. Akimov, Brendan Smith, Thomas A. Niehaus

"""


import os, sys, math, re, struct, copy
import numpy as np

from liblibra_core import *
import util.libutil as comn

from libra_py import units
from libra_py import scan
from libra_py import regexlib as rgl
import libra_py.citools.ci as ci
from libra_py import data_conv
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


def make_dftb_input_obsolete(dftb_input_template_filename, istate, timestep=-1):
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
    Parse an SPX.DAT file and construct forward and reverse mappings
    between orbital transitions and RPA indices.

    This function reads an SPX.DAT file (as produced e.g. by TDDFT/RPA
    implementations) and builds:

    1. rpa_map[ini, fin, spin] -> rpa_idx
       A forward lookup table mapping a transition from orbital `ini`
       to orbital `fin` with a given spin to its corresponding RPA index.

    2. transition_lookup[rpa_idx] -> (ini, fin, spin)
       A reverse lookup table mapping an RPA index back to the
       corresponding transition.

    3. actual_max_spin
       Indicates whether spin-down transitions are present:
           - 0 : restricted / unpolarized (only spin-up or no spin label)
           - 1 : unrestricted (both U and D spins present)

    All indices returned by this function are **0-based**, following
    standard Python conventions, even if the SPX.DAT file uses 1-based
    indexing.

    Parameters
    ----------
    filename : str
        Path to the SPX.DAT file to be parsed.

    Returns
    -------
    rpa_map : ndarray of shape (n_orb, n_orb, n_spin)
        Integer array mapping (ini, fin, spin) -> rpa_idx.
        Entries without a corresponding transition are set to -1.

    transition_lookup : ndarray of shape (n_rpa, 3)
        Integer array where each row is (ini, fin, spin) corresponding
        to a given rpa_idx.

    actual_max_spin : int
        Maximum spin index encountered:
            0 for restricted / spin-unpolarized calculations
            1 if spin-down ('D') transitions exist

    Notes
    -----
    Expected SPX.DAT format (simplified):

        idx  ...  ini  ...  fin  [spin]

    where:
        - idx  : RPA transition index (1-based)
        - ini  : initial orbital index (1-based)
        - fin  : final orbital index (1-based)
        - spin : optional spin label ('U' or 'D')

    Lines that are empty, contain '#', or contain '===' are ignored.
    Malformed lines are skipped silently.

    Spin mapping used internally:
        'U' -> 0
        'D' -> 1

    If no spin label is present, spin is assumed to be 0.

    Examples
    --------
    **Restricted / spin-unpolarized calculation**

    >>> rpa_map, trans_lookup, max_spin = read_spx_mappings("SPX.DAT")
    >>> print(max_spin)
    0

    Look up the RPA index for a transition from orbital 2 to 5:
    >>> idx = rpa_map[2, 5, 0]
    >>> print(idx)
    17

    Reverse lookup (from RPA index to transition):
    >>> ini, fin, spin = trans_lookup[17]
    >>> print(ini, fin, spin)
    2 5 0

    **Unrestricted / spin-polarized calculation**

    >>> rpa_map, trans_lookup, max_spin = read_spx_mappings("SPX.DAT")
    >>> print(max_spin)
    1

    Spin-up transition (U):
    >>> idx_u = rpa_map[1, 4, 0]

    Spin-down transition (D):
    >>> idx_d = rpa_map[1, 4, 1]

    Reverse lookup:
    >>> trans_lookup[idx_d]
    array([1, 4, 1])

    **Iterating over all RPA transitions**

    >>> for idx, (ini, fin, spin) in enumerate(trans_lookup):
    ...     print(idx, ini, fin, spin)

    See Also
    --------
    TDDFT, RPA, Casida equations, transition-density analysis
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
    Parse a tagged ASCII data file and return its contents as a dictionary.

    The file is assumed to consist of *tagged data blocks*, where each block
    starts with a header line describing the data, followed by one or more
    lines containing the numerical values.

    File format
    -----------
    Each data block has the form:

        <tag> : <dtype> : <dim> : <shape>
        <data values ...>

    where:

    * ``tag``   : string identifier used as the dictionary key
    * ``dtype`` : data type specifier (e.g. ``integer`` or ``real``)
    * ``dim``   : dimensionality of the data
                  - 0 → scalar
                  - 1 → vector
                  - 2 → matrix (or higher-rank array if supported)
    * ``shape`` : comma-separated list of dimensions
                  (ignored for ``dim = 0``)

    Example
    -------
    Scalar value:
        mermin_energy : real : 0 :
        -123.456D+00

    2D array:
        end_coords : real : 2 : 3,23
        1.0D+00  2.0D+00  ...
        ...

    Parameters
    ----------
    file_path : str
        Path to the tagged ASCII file to be parsed.

    Returns
    -------
    data_dict : dict
        Dictionary mapping tag names to parsed data.

        * Scalars (``dim = 0``) are returned as Python ``int`` or ``float``.
        * Arrays (``dim ≥ 1``) are returned as ``numpy.ndarray`` objects,
          reshaped according to the specified shape.

    Usage
    -----
    >>> results = parse_tagged_file("autotes.tag")
    >>> energy = results["mermin_energy"]      # float or int
    >>> coords = results["end_coords"]         # NumPy array

    Notes
    -----
    * Empty lines are ignored.
    * Numerical values written using Fortran ``D`` exponents (e.g. ``1.0D+03``)
      are automatically converted to standard ``E`` notation.
    * Array data are reshaped using **Fortran (column-major) ordering**
      (``order='F'``), which matches how Fortran writes multi-dimensional arrays.
    * If fewer data values than expected are encountered, parsing stops at the
      next header line as a safety measure.
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
    Read molecular orbital (MO) eigenvector matrices from a binary file
    written in Fortran unformatted (stream) style.

    The function always returns a 3D NumPy array with shape
    ``(nspin, ndim, ndim)``, where ``nspin`` is 1 for non–spin-polarized
    calculations and 2 for spin-polarized ones.

    Binary layout
    -------------
    The file is assumed to contain MO coefficient matrices written in
    **column-major (Fortran) order** with the following layout:

    * Record 1 (spin-up or total):
        - 4-byte integer header (typically a Fortran record marker or
          auxiliary integer; its value is not used)
        - ``ndim × ndim`` double-precision floating-point numbers (float64)

    * Record 2 (spin-down, only if ``spin_polarized=True``):
        - ``ndim × ndim`` double-precision floating-point numbers (float64)
        - No leading integer header is assumed for this record

    Parameters
    ----------
    filename : str
        Path to the binary file containing the MO eigenvector matrices.
    ndim : int
        Dimension of the MO matrix (number of basis functions or orbitals).
        Each MO matrix has shape ``(ndim, ndim)``.
    spin_polarized : bool, optional
        If ``True``, the file is assumed to contain separate spin-up and
        spin-down MO matrices. If ``False`` (default), only a single
        matrix (spin-up or total) is read.

    Returns
    -------
    mo_matrices : numpy.ndarray, shape (nspin, ndim, ndim)
        Molecular orbital coefficient matrices.
        * ``nspin = 1`` for non–spin-polarized calculations
        * ``nspin = 2`` for spin-polarized calculations

    Notes
    -----
    * The reshape uses ``order='F'`` because the matrices are written
      column-by-column by Fortran. Using NumPy's default row-major order
      would produce an incorrect matrix.
    * No validation is performed on the integer header; it is read and
      discarded.
    * The function assumes the file uses native endianness and 8-byte
      floating-point values (``float64``).
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
    Read XplusY.BIN written in Fortran *stream unformatted* mode and
    extract excitation vectors and energies.

    This routine parses the XplusY.BIN file typically produced by
    TDDFT / RPA / Casida solvers, where the file contains a header
    followed by per-root data blocks with optional energies and
    X+Y vectors.

    The binary layout is assumed to be **Fortran stream unformatted**
    (i.e., no record markers).

    Data type conventions (matching Fortran code):
        - dp   : REAL*8   (float64, 8 bytes)
        - ii   : INTEGER*4 (int32, 4 bytes)
        - sign : CHARACTER*1 (1 byte)

    Parameters
    ----------
    filename : str
        Path to the `XplusY.BIN` file.

    Returns
    -------
    nmat : int
        Length of each X+Y vector (number of matrix elements).

    nexc : int
        Number of excitation roots stored in the file.

    results : list of dict
        One dictionary per excitation root with keys:

        - 'root'   : int
            Root index (as stored in the file; typically 1-based).
        - 'sign'   : str
            Sign character associated with the root (usually '+' or '-').
        - 'energy' : float or None
            Excitation energy. If the Fortran code wrote a '-' instead
            of a floating-point value, this is set to None.
        - 'vector' : ndarray, shape (nmat,)
            X+Y vector for this root. Empty if `energy is None`.

    Notes
    -----
    File structure (byte-level, sequential):

    1. Global header
       - nmat : int32 (4 bytes)
       - nexc : int32 (4 bytes)

    2. For each excitation root:
       - root index (ii) : int32 (4 bytes)
       - sign            : char*1 (1 byte)

       - energy field:
           * Either a single '-' character (1 byte), indicating
             an invalid or skipped root
           * OR sqrt(eigenvalue) written as float64 (8 bytes)

       - X+Y vector:
           * Present only if energy is written
           * nmat × float64 (8 × nmat bytes)

    Because this is a *stream* binary file, detecting whether an
    energy is present requires peeking at the next byte and
    conditionally rewinding the file pointer.

    Examples
    --------
    **Basic usage**

    >>> nmat, nexc, results = read_xplusy_binary("XplusY.BIN")
    >>> print(nmat, nexc)
    120 10

    **Access energy and vector of the first root**

    >>> root0 = results[0]
    >>> print(root0['root'], root0['sign'], root0['energy'])
    1 + 3.457812

    >>> xpy = root0['vector']
    >>> print(xpy.shape)
    (120,)

    **Handling missing energies**

    >>> for r in results:
    ...     if r['energy'] is None:
    ...         print(f"Root {r['root']} has no energy")

    **Iterating over all valid roots**

    >>> for r in results:
    ...     if r['energy'] is not None:
    ...         norm = np.linalg.norm(r['vector'])
    ...         print(r['root'], r['energy'], norm)

    **Typical TDDFT / RPA workflow**

    >>> nmat, nexc, roots = read_xplusy_binary("XplusY.BIN")
    >>> energies = np.array([r['energy'] for r in roots if r['energy'] is not None])
    >>> vectors  = np.array([r['vector'] for r in roots if r['energy'] is not None])

    See Also
    --------
    read_spx_mappings : Mapping between orbital transitions and RPA indices
    Casida equation, TDDFT, RPA, excitation vectors
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
    Read an ASCII overlap matrix written in Fortran format and return it as
    a NumPy array.

    The file is assumed to contain a square overlap matrix stored as a
    flat list of numbers in **column-major (Fortran) order**, preceded by
    a fixed-length header.

    Parameters
    ----------
    filename : str
        Path to the ASCII file containing the overlap matrix.
        The numerical values are assumed to be written in a format such as
        ES24.15 (or similar), which `numpy.loadtxt` can parse automatically.
    n_orb : int
        Number of orbitals. The overlap matrix dimension is
        (n_orb, n_orb).

    Returns
    -------
    matrix : numpy.ndarray, shape (n_orb, n_orb)
        Overlap matrix reconstructed in the correct two-dimensional form.

    Notes
    -----
    * The first 5 lines of the file are treated as a header and skipped.
    * The remaining data are read as a one-dimensional array.
    * The reshape uses `order='F'` because Fortran writes matrices
      column-by-column, whereas NumPy defaults to row-major (C) order.
      Using the wrong order would result in a transposed or scrambled matrix.
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
    Parse a NACV.DAT file and return nonadiabatic coupling vectors (NACVs)
    as a 4D NumPy array.

    The returned array has the shape:

        nacv[state_i, state_j, dim, atom_index]

    where:
        * state_i, state_j : electronic state indices as they appear in
          the file (typically 1-based if the file is 1-based)
        * dim              : Cartesian component (0 = X, 1 = Y, 2 = Z)
        * atom_index       : atom index (0-based)

    File format
    -----------
    The file is assumed to consist of repeated blocks of the form:

        state_i  state_j
        x_1  y_1  z_1
        x_2  y_2  z_2
        ...
        x_N  y_N  z_N

    where N is the number of atoms. Each block gives the NACV
    ⟨ψ_i | ∇_R | ψ_j⟩ for all atoms.

    Parameters
    ----------
    filename : str
        Path to the NACV.DAT file.

    Returns
    -------
    nacv : numpy.ndarray, shape (nstate+1, nstate+1, 3, natoms)
        Nonadiabatic coupling vectors.

        The array is sized as ``(max_state + 1, max_state + 1, 3, natoms)``
        so that index ``n`` corresponds directly to electronic state ``n``
        without shifting (i.e., state indices are *not* converted to 0-based).

    Usage example
    -------------
    >>> nacv = read_nacv("NACV.DAT")
    >>> x_comp_atom_0 = nacv[1, 2, 0, 0]  # State 1 → 2, X component, atom 1

    Notes
    -----
    * The number of atoms is inferred from the length of the coordinate block.
    * Cartesian components are stored as:
        - 0 → X
        - 1 → Y
        - 2 → Z
    * Atom indices are 0-based.
    * The function enforces the standard NACV antisymmetry:
        ⟨ψ_i | ∇ | ψ_j⟩ = −⟨ψ_j | ∇ | ψ_i⟩
    * All unassigned state pairs are initialized to zero.
    * No unit conversion is performed; values are returned in the units
      used in the input file (typically Bohr⁻¹ or a.u.).
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



def create_odin_inp(params):
    """
    Generate the input string for the ODIN overlap program.

    This function constructs the text input required by the ODIN code
    (https://github.com/thomas-niehaus/odin) for computing overlap and
    Hamiltonian-related quantities from Slater–Koster files and a GEN
    geometry file.

    The function:

        1. Reads the provided GEN file.
        2. Extracts the list of chemical elements from the second line.
        3. Appends the corresponding maximum angular momentum (ℓ_max)
           values for each element.
        4. Returns the complete ODIN input as a formatted string.

    Parameters
    ----------
    params : dict
        Dictionary containing ODIN input parameters.

        Optional keys
        -------------
        filename : str, optional
            Path to the GEN-format geometry file for which overlaps
            will be computed.
            Default: "x1.gen"

        slakos_prefix : str, optional
            Directory prefix containing the Slater–Koster (.skf) files.
            Default: "./"

        max_ang_mom : dict[str, int], optional
            Dictionary mapping element symbols to their maximum angular
            momentum quantum number ℓ_max.
            Example: {"H": 1, "C": 2}
            Default: {"H": 1, "C": 2}

    Returns
    -------
    str
        The formatted ODIN input file content as a single string.
        This string can be written directly to a file and used as
        standard input for the ODIN executable.

    Raises
    ------
    SystemExit
        If an element present in the GEN file is not found in the
        ``max_ang_mom`` dictionary.

    Notes
    -----
    - The GEN file must follow the standard DFTB+ GEN format.
    - The second line of the GEN file is assumed to contain the list
      of element symbols.
    - The function does not validate physical consistency of the
      angular momentum assignments.
    - The returned string is not automatically written to disk.
    """

    filename = params.get("filename", "x1.gen") # gen file
    slakos_prefix = params.get("slakos_prefix", "./") # prefix to the file where the Slater-Koster files are located
    max_l = params.get("max_ang_mom", {"H":1, "C":2} )  # maximal angular momentum numbers for elements

    res = F"'{filename}'\n'{slakos_prefix}'\n'-'\n'.skf'\n"
    f = open(filename)
    A = f.readlines()
    f.close()

    elts = A[1].split()
    print(elts)
    for elt in elts:
        if elt in max_l.keys():
            res = F"{res} {max_l[elt]} "
        else:
            print(F"Element: {elt} is not present in the input `max_l` parameter\n max_l = {max_l}")
            sys.exit(0)

    return res



def run_odin(params):
    """
    Execute the ODIN program to compute atomic orbital overlap matrices.

    This function prepares the required ODIN input file ("odin.inp"),
    runs the ODIN executable, and computes overlap-related quantities
    based on the provided GEN geometry file and Slater–Koster parameters.

    Specifically, the function:

        1. Generates the ODIN input string using ``create_odin_inp``.
        2. Writes the input to the file "odin.inp".
        3. Executes the ODIN program using standard input redirection.
        4. Produces overlap output files (e.g., "oversqr.dat").

    ODIN is a post-processing tool for DFTB that computes atomic orbital
    overlap matrices and related quantities from Slater–Koster files and
    molecular geometries.

    Parameters
    ----------
    params : dict
        Dictionary containing ODIN execution parameters.

        Optional keys
        -------------
        ODIN_EXE : str, optional
            Path or name of the ODIN executable.
            Must be accessible from the system PATH or specified explicitly.
            Default: "odin"

        filename : str, optional
            Path to the GEN-format geometry file used for overlap calculations.
            Default: "x1.gen"

        slakos_prefix : str, optional
            Directory prefix containing the Slater–Koster (.skf) files.
            Default: "./"

        max_ang_mom : dict[str, int], optional
            Dictionary mapping element symbols to their maximum angular
            momentum quantum number ℓ_max.
            Example: {"H": 1, "C": 2}
            Default: {"H": 1, "C": 2}

    Returns
    -------
    None
        This function does not return a value. It executes ODIN and produces
        output files in the current working directory.

    Side Effects
    ------------
    - Creates or overwrites the file "odin.inp".
    - Executes an external program via ``os.system``.
    - Generates ODIN output files, typically including:
        oversqr.dat : atomic orbital overlap matrix (square format)

    Notes
    -----
    - The ODIN executable must be installed and accessible.
    - The GEN file and Slater–Koster files must be present and consistent.
    - No error handling is performed for failed ODIN execution.
    - Output files are written to the current working directory.

    References
    ----------
    ODIN GitHub repository:
    https://github.com/thomas-niehaus/odin
    """

    EXE = params.get("ODIN_EXE", "odin")
    inp = create_odin_inp(params)
    f = open("odin.inp", "w")
    f.write(F"{inp}")
    f.close()
    os.system(F"{EXE} < odin.inp")





def write_dftb_gen(filename, atom_labels, coordinates, periodic=False):
    """
    Write a DFTB+ .gen file.

    Parameters
    ----------
    filename : str
        Output file name.
    atom_labels : list of str
        Atomic symbols, length N.
    coordinates : array-like, shape (3N, 1)
        Flattened Cartesian coordinates:
        [x1,y1,z1,x2,y2,z2,...]^T  (Angstrom)
    periodic : bool
        True -> periodic system (S)
        False -> cluster (C)

    Example
    -------
    > atom_labels = ["O", "H", "H"]
    > 
    > coordinates = np.array([
    > 0.000000, 0.000000, 0.000000,
    > 0.758602, 0.000000, 0.504284,
    > -0.758602, 0.000000, 0.504284
    > ]).reshape(-1, 1)
    >
    > write_dftb_gen("water.gen", atom_labels, coordinates)
    """

    coords = np.asarray(coordinates)

    if coords.ndim != 2 or coords.shape[1] != 1:
        raise ValueError("Coordinates must have shape (3N, 1)")

    n_atoms = len(atom_labels)

    if coords.shape[0] != 3 * n_atoms:
        raise ValueError("Coordinate length must be 3 * number of atoms")

    # Reshape to (N, 3)
    coords = coords.reshape(n_atoms, 3)

    # Unique species in order of appearance
    species = []
    for label in atom_labels:
        if label not in species:
            species.append(label)

    species_index = {s: i + 1 for i, s in enumerate(species)}
    system_type = "S" if periodic else "C"

    with open(filename, "w") as f:

        # Header
        f.write(f"{n_atoms} {system_type}\n")
        f.write(" ".join(species) + "\n\n")

        # Atom lines
        for i, (label, coord) in enumerate(zip(atom_labels, coords)):
            idx = i + 1
            sp_idx = species_index[label]
            x, y, z = coord
            f.write(f"{idx:5d} {sp_idx:3d} {x:15.8f} {y:15.8f} {z:15.8f}\n")

    print(f"DFTB+ GEN file written to {filename}")



def make_dftb_input(params):
    """
    Generate a DFTB+ input file ("dftb_in.hsd") for ground-state SCC-DFTB and
    excited-state TD-DFTB (Casida formalism) calculations.

    This function constructs a complete DFTB+ input file using parameters
    provided in the ``params`` dictionary. Missing parameters are automatically
    filled with reasonable defaults. The generated file is written to the
    current working directory as "dftb_in.hsd".

    The input includes the following sections:

        - Geometry (GenFormat)
        - Driver
        - Hamiltonian (SCC-DFTB)
        - ExcitedState (Casida TD-DFTB)
        - Options
        - Analysis

    Parameters
    ----------
    params : dict
        Dictionary containing DFTB+ input parameters. All keys are optional.
        The following keys are recognized:

        Geometry and general:
        ---------------------
        gen_file : str, optional
            Path to the GEN format geometry file.
            Default: "x1.gen"

        Driver : str, optional
            Driver block specifying geometry optimization, MD, or single-point
            calculation. Must be valid DFTB+ HSD syntax.
            Default: "{}" (single-point calculation)

        Hamiltonian parameters:
        -----------------------
        sk_prefix : str, optional
            Path prefix for Slater-Koster files.
            Default: "mio/FinalSK/"

        MaxAngularMomentum : str, optional
            Block specifying maximum angular momentum per element.
            Must be valid HSD syntax.
            Default:
                {
                    O = "p"
                    H = "s"
                }

        Filling : str, optional
            Electron filling specification block.
            Default:
                Fermi {
                    Temperature [K] = 40
                }

        Excited state (Casida TD-DFTB):
        -------------------------------
        Symmetry : str, optional
            Spin symmetry of excitations ("Singlet" or "Triplet").
            Default: "Singlet"

        NrOfExcitations : int, optional
            Number of excited states to compute.
            Default: 5

        StateOfInterest : int, optional
            Index of the excited state of interest.
            Default: 1

        WriteSPTransitions : str, optional
            Whether to write single-particle transitions ("Yes"/"No").
            Default: "Yes"

        WriteXplusY : str, optional
            Whether to write X+Y excitation vectors ("Yes"/"No").
            Default: "Yes"

        WriteXplusYAscii : str, optional
            Whether to write ASCII excitation vectors ("Yes"/"No").
            Default: "Yes"

        StateCouplings : str, optional
            State coupling specification in HSD format.
            Default: "{0 2}"

        Output and analysis:
        --------------------
        WriteAutotestTag : str, optional
            Enable autotest tag output ("Yes"/"No").
            Default: "Yes"

        WriteHS : str, optional
            Write Hamiltonian and overlap matrices ("Yes"/"No").
            Default: "No"

        WriteEigenvectors : str, optional
            Write eigenvectors ("Yes"/"No").
            Default: "Yes"

        EigenvectorsAsText : str, optional
            Write eigenvectors in ASCII format ("Yes"/"No").
            Default: "Yes"

        PrintForces : str, optional
            Print atomic forces ("Yes"/"No").
            Default: "Yes"

    Returns
    -------
    None
        Writes the file "dftb_in.hsd" to disk. No value is returned.

    Notes
    -----
    - The function assumes that all string parameters are valid DFTB+ HSD syntax.
    - This function does not validate physical correctness of the parameters.
    - Existing "dftb_in.hsd" file will be overwritten.
    - Designed for use in automated workflows such as molecular dynamics,
      excited-state dynamics, and high-throughput simulations.

    Examples
    --------
    Minimal usage:

    >>> params = {
    >>>     "gen_file": "geom.gen",
    >>>     "NrOfExcitations": 10,
    >>>     "StateOfInterest": 2
    >>> }
    >>> make_dftb_input(params)

    Custom angular momentum:

    >>> params = {
    >>>     "MaxAngularMomentum": '''
    >>>     {
    >>>         C = "p"
    >>>         O = "p"
    >>>         H = "s"
    >>>     }
    >>>     '''
    >>> }
    >>> make_dftb_input(params)

    """

    # Read in the parameters
    gen_file = params.get("gen_file", "x1.gen")
    sk_prefix = params.get("sk_prefix", "mio/FinalSK/")
    Driver = params.get("Driver", "{}")
    MaxAngularMomentum = params.get("MaxAngularMomentum",
                      """{ O = "p" 
                          H = "s" 
                     }
                     """)
    Symmetry = params.get("Symmetry", "Singlet")
    NrOfExcitations = params.get("NrOfExcitations", 5)
    StateOfInterest = params.get("StateOfInterest", 1)
    WriteSPTransitions = params.get("WriteSPTransitions", "Yes")
    WriteXplusY = params.get("WriteXplusY", "Yes")
    WriteXplusYAscii = params.get("WriteXplusYAscii", "Yes")
    StateCouplings = params.get("StateCouplings", "{0 2}")

    WriteAutotestTag = params.get("WriteAutotestTag", "Yes")
    WriteHS = params.get("WriteHS", "No")
    WriteEigenvectors = params.get("WriteEigenvectors", "Yes")
    EigenvectorsAsText = params.get("EigenvectorsAsText", "Yes")
    PrintForces = params.get("PrintForces", "Yes")

    Filling = params.get("Filling",
                         """Fermi {
                         Temperature [K] = 40    
                         }
                         """
                        )


    scf_template = F"""Geometry = GenFormat {{  
  <<< "{gen_file}"
}}
Driver = { Driver }

Hamiltonian = DFTB {{
    SCC = Yes
    SCCTolerance = 1.0E-10
    MaxAngularMomentum = { MaxAngularMomentum }
    SlaterKosterFiles = Type2FileNames {{
        Prefix = "{sk_prefix}"
        Separator = "-"
        Suffix = ".skf"
    }}
    Filling = {Filling}
}}

ExcitedState {{
    Casida {{
       Symmetry = {Symmetry}
       NrOfExcitations = {NrOfExcitations}
       StateOfInterest = {StateOfInterest}
       WriteSPTransitions = {WriteSPTransitions}
       WriteXplusY = {WriteXplusY}
       #WriteXplusYAscii = {WriteXplusYAscii}
       StateCouplings = {StateCouplings}
       Diagonaliser = Stratmann{{}}
    }}
}}
Options {{
  WriteAutotestTag = {WriteAutotestTag}
  WriteHS = {WriteHS}
}}
Analysis {{
  WriteEigenvectors = {WriteEigenvectors}
  EigenvectorsAsText = {EigenvectorsAsText}
  PrintForces = {PrintForces}
}}
"""

    f = open("dftb_in.hsd", "w")
    f.write(scf_template)
    f.close()



def run_dftb(coords, params):
    """
    Run a DFTB+ calculation for a given molecular geometry.

    This function prepares and executes a DFTB+ calculation in a dedicated
    working directory. It performs the following steps:

        1. Creates (if necessary) a working directory.
        2. Writes a GEN-format geometry file using the provided coordinates.
        3. Generates the DFTB+ input file ("dftb_in.hsd").
        4. Executes the DFTB+ binary.
        5. Returns to the original working directory.

    The function is designed for automated workflows such as molecular
    dynamics, nonadiabatic dynamics, or high-throughput simulations.

    Parameters
    ----------
    coords : ndarray
        Cartesian atomic coordinates. Expected shape is (3*N,) or (3*N, 1),
        where N is the number of atoms. Coordinates are assumed to be in
        atomic units (Bohr) and are internally converted to Angstroms
        before writing the GEN file.

    params : dict
        Dictionary containing run configuration parameters.

        Required keys:
        --------------
        labels : list of str
            Atomic symbols (length N), e.g., ["O", "H", "H"].

        Optional keys:
        --------------
        exe : str, optional
            Name or path of the DFTB+ executable.
            Default: "dftb+"

        dftb_run_params : dict, optional
            Dictionary of parameters passed to ``make_dftb_input`` to
            construct the "dftb_in.hsd" file.
            Default: {}

        working_directory : str, optional
            Directory where the calculation will be executed.
            If it does not exist, it will be created.
            Default: "wd"

        gen_file : str, optional
            Name of the GEN-format geometry file.
            Default: "x1.gen"

    Returns
    -------
    None
        The function executes DFTB+ and produces output files in the
        working directory. No values are returned.

    Notes
    -----
    - The function assumes that the DFTB+ executable is available in the
      system PATH unless an explicit path is provided via ``exe``.
    - Existing files in the working directory may be overwritten.
    - No error handling is performed for failed DFTB+ runs.
    - The function changes the current working directory during execution
      and restores it afterward (relative path "../").

    Side Effects
    ------------
    - Creates or modifies files in the specified working directory.
    - Executes an external system command via ``os.system``.

    Example
    -------
    >>> coords = np.array([
    >>>     0.0, 0.0, 0.0,
    >>>     0.0, 0.0, 1.0,
    >>>     1.0, 0.0, 0.0
    >>> ])
    >>> params = {
    >>>     "labels": ["O", "H", "H"],
    >>>     "working_directory": "run1",
    >>>     "dftb_run_params": {
    >>>         "NrOfExcitations": 5,
    >>>         "StateOfInterest": 1
    >>>     }
    >>> }
    >>> run_dftb(coords, params)
    """


    # Get the parameters
    labels = params["labels"]
    exe = params.get("exe", "dftb+")
    dftb_run_params = params.get("dftb_run_params", {}  )
    wd = params.get("working_directory", "wd")
    gen_file = params.get("gen_file", "x1.gen")

    natoms = len(labels)
    ndof = 3 * natoms

    # Create working directory, if doesn't exist
    if not os.path.exists(wd):
        os.mkdir(wd)

    # Go into that directory
    os.chdir(wd)

    # Make the gen file
    write_dftb_gen(gen_file, labels, coords/units.Angst, periodic=False)

    # Make the input file itself
    make_dftb_input(dftb_run_params)
    
    # Run DFTB+ calculations
    os.system(F"{exe}")

    # Go back to the original directory
    os.chdir("../")


def read_dftb_orbital_info(params_):
    """
    Read molecular orbital (MO) information, configurations, and CI amplitudes
    from DFTB+ excited-state output files.

    This function extracts the active-space MO coefficients and excited-state
    configuration interaction (CI) information from DFTB+ linear-response
    (Casida/RPA) output files.

    Parameters
    ----------
    params_ : dict
        Dictionary of input parameters. Recognized keys:

        source_directory : str, optional, default="calc"
            Directory containing the DFTB+ output files:
                - SPX.DAT      : orbital excitation mappings
                - eigenvec.bin : MO coefficient matrix
                - XplusY.DAT   : excitation energies and CI vectors

        orbital_space : list of int, optional, default=None
            List of molecular orbital indices (1-based indexing) defining
            the active orbital space. If None, all orbitals are used.

        nstates : int, optional, default=2
            Total number of electronic states to extract, including the ground
            state. Since DFTB+ Casida output contains only excited states,
            this function reads at most (nstates - 1) excited states.

        ci_threshold : float, optional, default=0.0
            Defines the minimial magnitude of CI amplitudes to be included in
            definition of TD-DFT excited states

    Returns
    -------
    info : dict
        Dictionary containing system and active-space metadata:

        nocc : int
            Number of occupied orbitals.

        nelec : int
            Number of electrons (assumes closed-shell: nelec = 2 * nocc).

        nao : int
            Number of atomic orbitals.

        nmo : int
            Total number of molecular orbitals.

        nci : int
            Number of excited states actually read.

        nact : int
            Number of active orbitals.

        actual_orbital_space : list of int
            List of active orbital indices (1-based).

        min_occ, max_occ : int
            Minimum and maximum occupied orbital indices involved in excitations.

        min_vir, max_vir : int
            Minimum and maximum virtual orbital indices involved in excitations.

    MOs : np.ndarray, shape (nao, nact)
        Active-space MO coefficient matrix:

            MOs[μ, i] = coefficient of atomic orbital μ in molecular orbital i

        Only active orbitals defined by `orbital_space` are included.

    data : list
        Container with CI information:

        data[0] : np.ndarray, shape (nci,)
            Excited-state energies in atomic units (Hartree).

        data[1] : np.ndarray, shape (nci, nconf, 2)
            Configuration list for each excited state.

            Each configuration is a pair:
                [i, j] = excitation from occupied orbital i → virtual orbital j

            Orbital indices use 1-based indexing.

        data[2] : np.ndarray, shape (nci, nconf)
            CI amplitudes corresponding to each configuration.

    Notes
    -----
    File dependencies (in source_directory)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SPX.DAT
        Contains excitation mappings between occupied and virtual orbitals.

    eigenvec.bin
        Binary file containing MO coefficients.

    XplusY.DAT
        ASCII file containing excitation energies and CI eigenvectors.

    Indexing conventions
    ~~~~~~~~~~~~~~~~~~~~
    - Orbital indices returned in `info` and `data` use 1-based indexing.
    - Internally, Python 0-based indexing is used for array access.

    Assumptions
    ~~~~~~~~~~~
    - Closed-shell system (nelec = 2 * nocc).
    - Spin-restricted DFTB calculation.
    - Single excitations only (Casida formalism).

    """

    # Copy parameters to avoid modifying input dictionary
    params = dict(params_)

    nstates = params.get("nstates", 2)
    orbital_space = params.get("orbital_space", None)
    wd = params.get("source_directory", "calc")
    ci_threshold = params.get("ci_threshold", 0.0)

    # ================= SPX: excitation mappings =================

    rpa_map, rpa_lookup, max_spin = read_spx_mappings(
        f"{wd}/SPX.DAT"
    )

    ndim = rpa_map.shape[1]
    nmo = ndim
    nao = nmo 

    # Define active orbital space (1-based indexing)
    if orbital_space is None:
        actual_orbital_space = list(range(1, nmo + 1))
    else:
        actual_orbital_space = list(orbital_space)

    nact = len(actual_orbital_space)

    # ================= Read MO coefficients =====================

    mo = read_mo_matrix(
        f"{wd}/eigenvec.bin",
        ndim,
        spin_polarized=False
    )
    # mo shape: (spin, AO, MO)

    active_space_indx = [i - 1 for i in actual_orbital_space]

    # Extract active-space MOs (closed shell → spin index 0)
    MOs = np.take(mo[0], active_space_indx, axis=1) # to avoid axis reordering due to
                                                    # advanced indexing

    # ================= Read CI amplitudes =======================

    nmat, nexc, results = read_xplusy_ascii(
        f"{wd}/XplusY.DAT"
    )

    # Number of excited states to include
    nci = min(nexc, nstates - 1)

    # Configuration list from SPX lookup (convert to 1-based indexing)
    all_confs = [[v[0] + 1, v[1] + 1] for v in rpa_lookup]
    nconf = len(all_confs)

    # Orbital ranges
    occ_indices = [v[0] + 1 for v in rpa_lookup]
    vir_indices = [v[1] + 1 for v in rpa_lookup]

    min_occ = min(occ_indices)
    max_occ = max(occ_indices)
    min_vir = min(vir_indices)
    max_vir = max(vir_indices)


    E_CI = []
    CI = []
    confs = []

    min_occ, max_occ = 10000000, 0
    min_vir, max_vir = 10000000, 0
    for i in range(nci):

        vec = np.array(results[i]["vector"], dtype=float)   # shape (nconf,)
        conf_arr = np.array(all_confs, dtype=int)           # shape (nconf, 2)

        # Find important configurations
        mask = np.abs(vec) > ci_threshold

        filtered_vec = vec[mask]
        filtered_confs = conf_arr[mask]

        for iconf in filtered_confs:
            if iconf[0]<min_occ:
                min_occ = iconf[0] 
            if iconf[0]>max_occ:
                max_occ = iconf[0]
            if iconf[1]<min_vir:
                min_vir = iconf[1]
            if iconf[1]>max_vir:
                max_vir = iconf[1]

        E_CI.append(results[i]["energy"])
        CI.append(filtered_vec)
        confs.append(filtered_confs)

    E_CI = np.array(E_CI, dtype=float)

    # Do NOT convert CI and confs to regular numpy arrays,
    # because each state may have different number of configurations
    data = [E_CI, confs, CI]

    # ================= Build info dictionary ====================

    nocc = max_occ
    nelec = 2 * nocc

    info = {
        "nocc": nocc,
        "nelec": nelec,
        "nao": nao,
        "nmo": nmo,
        "nci": nci,
        "nact": nact,
        "actual_orbital_space": actual_orbital_space,
        "min_occ": min_occ,
        "max_occ": max_occ,
        "min_vir": min_vir,
        "max_vir": max_vir,
    }

    return info, MOs, data


class tmp:
    pass

def dftb_compute_adi(q, params, full_id):
    """
    Compute adiabatic Hamiltonian, vibronic Hamiltonian, and time-overlap matrices
    using DFTB+ electronic structure and ODIN orbital overlap calculations.

    This function performs a single electronic structure step for one trajectory,
    computes molecular orbitals (MOs), CI states, and their time-overlaps between
    consecutive time steps, and constructs the adiabatic vibronic Hamiltonian.

    This is designed for trajectory-based nonadiabatic dynamics methods such as:
        - FSSH
        - Ehrenfest dynamics
        - Mapping-based methods
        - Exact factorization / quantum trajectory methods

    Parameters
    ----------
    q : MATRIX
        Nuclear coordinates for all trajectories.

        Shape:
            (3*N_atoms, N_trajectories)

        Units:
            Bohr

        Column `itraj` corresponds to the coordinates of trajectory `itraj`.

    params : list of dict
        List of trajectory-specific parameter dictionaries.

        Each params[itraj] must contain:

        Required keys
        -------------
        orbital_space : list of int
            Active-space molecular orbital indices (1-based indexing).

        atom_labels : list of str
            Atomic symbols, e.g. ["O", "H", "H"]

        Optional keys
        -------------

        dt : float, default = 41.0
            Nuclear time step in atomic units (1 fs ≈ 41 a.u.)

        working_directory : str, default = "wd"
            Directory used for DFTB+ and ODIN calculations.

        dftb_exe : str, default = "dftb+"
            Path to DFTB+ executable.

        odin_exe : str, default = "odin"
            Path to ODIN executable.

        dftb_run_params : dict
            Parameters controlling DFTB+ calculation, including:

                NrOfExcitations : int
                    Number of excited states requested from Casida calculation

                gen_file : str
                    Geometry file name

                sk_prefix : str
                    Slater–Koster files directory

        nelec_act_space : int, optional
            Number of electrons in active space.

        ci_threshold : float, default = 0.01
            CI amplitude threshold for truncating configurations.

        is_first_time : bool, default = True
            Indicates whether this is the first step of trajectory.

        Internal state keys (automatically updated)
        -------------------------------------------
        MO_prev : ndarray (nao, nact)
        data_prev : list
        coordinates_prev : ndarray (3*N_atoms, 1)

    full_id : int
        Encoded trajectory identifier.

    Returns
    -------
    obj : object
        Object containing adiabatic electronic properties.

        Attributes
        ----------

        ham_adi : CMATRIX (nstates, nstates)
            Adiabatic Hamiltonian matrix (Hartree)

        hvib_adi : CMATRIX (nstates, nstates)
            Vibronic Hamiltonian including derivative couplings

        time_overlap_adi : CMATRIX (nstates, nstates)
            Time-overlap matrix between electronic states:

                S_ij(t, t+dt) = <Ψ_i(t) | Ψ_j(t+dt)>

        basis_transform : CMATRIX (nstates, nstates)
            Basis transformation matrix (currently identity)

        overlap_adi : CMATRIX (nstates, nstates)
            Placeholder for state overlap matrix

    Notes
    -----

    Electronic structure workflow
    -----------------------------

    1. Run DFTB+ to compute:
        - Molecular orbitals
        - Excitation energies
        - CI coefficients

    2. Run ODIN to compute orbital overlaps between consecutive geometries

    3. Construct:
        - MO overlap matrix
        - CI state overlap matrix
        - Adiabatic Hamiltonian
        - Vibronic Hamiltonian

    Vibronic Hamiltonian definition
    -------------------------------

    Hvib_ij = E_i δ_ij - i d_ij

    where

        d_ij = (S_ij - S_ji) / (2 dt)

    is the time-derivative coupling.

    Units
    -----

    Energies:
        Hartree

    Time:
        atomic units

    Coordinates:
        Bohr

    Overlaps:
        dimensionless

    State indexing
    --------------

    state 0 = ground state
    state 1..nstates-1 = excited states

    """

    # ================= Decode trajectory index =================

    Id = Cpp2Py(full_id)
    itraj = Id[-1]

    # ================= Validate params structure =================

    if not isinstance(params, (list, tuple)):
        raise TypeError(
            "params must be a list (or tuple) of per-trajectory dictionaries"
        )

    if itraj >= len(params):
        raise IndexError(
            f"Trajectory index {itraj} out of range (len={len(params)})"
        )

    if not isinstance(params[itraj], dict):
        raise TypeError(
            f"params[{itraj}] must be a dictionary"
        )

    # ================= Extract nuclear coordinates =================

    coords = q.col(itraj)
    coordinates = data_conv.MATRIX2nparray(coords, float)
    # shape: (3*N_atoms, 1)
    # units: Bohr

    # ================= Read parameters =================

    dt = float(params[itraj].get("dt", 41.0))
    dftb_run_params = params[itraj].get("dftb_run_params", {})
    atom_labels = params[itraj]["atom_labels"]
    dftb_exe = params[itraj].get("dftb_exe", "dftb+")
    wd = params[itraj].get("working_directory", "wd")
    gen_file = dftb_run_params.get("gen_file", "x1.gen")
    sk_prefix = dftb_run_params.get("sk_prefix", "../mio/FinalSK/")
    odin_max_ang_mom = params[itraj].get(
        "odin_max_ang_mom",
        {"H": 1, "O": 2}
    )
    odin_exe = params[itraj].get("odin_exe", "odin")
    orbital_space = params[itraj]["orbital_space"]
    nstates = dftb_run_params.get("NrOfExcitations", 1) + 1
    is_first_time = params[itraj].get("is_first_time", True)
    nelec_act_space = params[itraj].get("nelec_act_space", None)
    ci_threshold = params[itraj].get("ci_threshold", 0.01)

    # ================= Run DFTB+ =================

    dftb_params = dftb_run_params
    make_dftb_input(dftb_params)
    prms1 = {
        "labels": atom_labels,
        "exe": dftb_exe,
        "dftb_run_params": dftb_params,
        "working_directory": wd,
        "gen_file": gen_file,
    }
    run_dftb(coordinates, prms1)

    # ================= Run ODIN overlap calculation =================

    if is_first_time:
        coordinates_prev = copy.deepcopy(coordinates)
    else:
        coordinates_prev = copy.deepcopy(
            params[itraj]["coordinates_prev"]
        )

    double_labels = atom_labels * 2

    double_coords = np.concatenate(
        (coordinates_prev, coordinates),
        axis=0
    )

    os.chdir(wd)
    write_dftb_gen( "doubled_geometry.gen",double_labels, double_coords / units.Angst )
    odin_params = {
        "ODIN_EXE": odin_exe,
        "filename": "doubled_geometry.gen",
        "slakos_prefix": sk_prefix,
        "max_ang_mom": odin_max_ang_mom,
    }

    create_odin_inp(odin_params)
    run_odin(odin_params)
    os.chdir("../")

    # ================= Read electronic structure =================

    read_params = {
        "nstates": nstates,
        "orbital_space": orbital_space,
        "source_directory": wd,
        "ci_threshold": ci_threshold,
    }

    info, MO_curr, data_curr = read_dftb_orbital_info(read_params)

    # ================= Read AO overlap matrix =================

    ndim = info["nmo"]
    S = read_overlap_matrix( f"{wd}/oversqr.dat", 2 * ndim )
    s_ao_curr = S[ndim:, ndim:]
    st_ao = S[:ndim, ndim:]

    # ================= Construct active space =================

    if nelec_act_space is None:
        active_space = info["actual_orbital_space"]
    else:
        min_indx = info["nocc"] - nelec_act_space // 2 + 1
        if min_indx > info["min_occ"]:
            raise ValueError("Active space too small for requested electrons")
        active_space = list( range(min_indx, info["nmo"] + 1) )

    # ================= Retrieve previous step data =================

    if is_first_time:
        MO_prev = copy.deepcopy(MO_curr)
        data_prev = copy.deepcopy(data_curr)
    else:
        MO_prev = copy.deepcopy(params[itraj]["MO_prev"])
        data_prev = copy.deepcopy(params[itraj]["data_prev"])

    # ================= Allocate output object =================

    obj = tmp()
    obj.ham_adi = CMATRIX(nstates, nstates)
    obj.nac_adi = CMATRIX(nstates, nstates)
    obj.hvib_adi = CMATRIX(nstates, nstates)
    obj.basis_transform = CMATRIX(nstates, nstates)
    obj.time_overlap_adi = CMATRIX(nstates, nstates)
    obj.overlap_adi = CMATRIX(nstates, nstates)

    # ================= Compute MO overlap =================

    st_mo_orb = MO_prev.T @ st_ao @ MO_curr
    st_mo = np.kron(np.eye(2), st_mo_orb)

    # ================= Compute CI overlaps =================

    ovlp_params = {
        "homo_indx": info["nocc"],
        "nocc": info["nocc"] - 1,
        "nvirt": info["nmo"] - info["nocc"],
        "nelec": info["nelec"],
        "nstates": nstates,
        "active_space": active_space,
    }

    st_ci = ci.overlap(st_mo, data_prev, data_curr,ovlp_params )

    # ================= Populate Hamiltonian =================

    for i in range(nstates):
        energy = 0.0
        if i > 0:
            energy = 0.5 * ( data_prev[0][i - 1] + data_curr[0][i - 1] )
        obj.ham_adi.set(i, i, energy)
        obj.hvib_adi.set(i, i, energy)
        obj.basis_transform.set(i, i, 1.0)

        for j in range(nstates):
            obj.time_overlap_adi.set( i, j, float(st_ci[i, j])  )

    # ================= Compute derivative couplings =================

    for i in range(nstates):
        for j in range(i + 1, nstates):
            dij = ( obj.time_overlap_adi.get(i, j) - obj.time_overlap_adi.get(j, i) ) / (2.0 * dt)
            obj.hvib_adi.set(i, j, -1j * dij)
            obj.hvib_adi.set(j, i, +1j * dij)

    # ================= Store state for next step =================
    params[itraj]["MO_prev"] = copy.deepcopy(MO_curr)
    params[itraj]["data_prev"] = copy.deepcopy(data_curr)
    params[itraj]["coordinates_prev"] = copy.deepcopy(coordinates)
    params[itraj]["is_first_time"] = False

    return obj


