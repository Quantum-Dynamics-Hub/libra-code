# *********************************************************************************
# * Copyright (C) 2024  Alexey V. Akimov
# *
# * This file is distributed under the terms of the GNU General Public License
# * as published by the Free Software Foundation, either version 3 of
# * the License, or (at your option) any later version.
# * See the file LICENSE in the root directory of this distribution
# * or <http://www.gnu.org/licenses/>.
# *
# *********************************************************************************/
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

if sys.platform == "cygwin":
    from cyglibra_core import *
elif sys.platform == "linux" or sys.platform == "linux2":
    from liblibra_core import *
import util.libutil as comn

from libra_py import units
from libra_py import scan
from libra_py import regexlib as rgl
from libra_py import data_conv

import libra_py.packages.cp2k.methods as CP2K_methods
import libra_py.workflows.nbra.mapping2 as mapping2
import libra_py.workflows.nbra.mapping3 as mapping3
import libra_py.workflows.nbra.step3 as step3

import libra_py.citools.slatdet as sd
import libra_py.citools.interfaces as interfaces


def make_mopac_input(mopac_input_filename, mopac_run_params, labels, coords):
    """
    This function creates an input file for MOPAC package using the
    parameters passed in the `mopac_input_params` dictionary

    Args:
        * mopac_input_filename ( string ): the name of the input file to create

        * mopac_run_params ( string ): the string containing the specification for the MOPAC run.
        E.g. one can use: "INDO C.I.=(6,3) CHARGE=0 RELSCF=0.000001 ALLVEC  WRTCONF=0.00  WRTCI=2"

        * labels (list of stings): element symbols for atoms in the system (N items), e.g.
         ["C", "H", "H", "H", "H"] for methane

        * coords ( MATRIX(3N, 1) ): Cartesian coordinates of all atoms ordered in triples x, y, z [ units: Bohr ]

    Returns:
        None :  just creates the files

    """

    # Create the actual output file
    mopac_input = open(mopac_input_filename, "w")

    mopac_input.write(F"{mopac_run_params}\n\n\n")

    nat = len(labels)  # how many atoms
    for i in range(nat):
        x = coords.get(3 * i + 0, 0) / units.Angst
        y = coords.get(3 * i + 0, 1) / units.Angst
        z = coords.get(3 * i + 0, 2) / units.Angst
        mopac_input.write(F"{labels[i]}   {x}  1  {y} 1   {z}  1\n")
    mopac_input.write("\n")
    mopac_input.close()


class tmp:
    pass


def run_mopac(coords, params_):
    """

    This function executes the MOPAC quantum chemistry calculations

    Args:
        coords ( MATRIX(ndof, 1) ): coordinates of the particle [ units: Bohr ]
        params ( dictionary ): model parameters

            * **params_["labels"]** ( list of strings ): the labels of atomic symbolc - for all atoms,
                and in a order that is consistent with the coordinates (in triples) stored in `q`.
                The number of this labels is `natoms`, such that `ndof` = 3 * `natoms`. [ Required ]
            * **params_["mopac_exe"]** ( string ):  the full path to `the mopac` executable [ defaut: "mopac" ]
            * **params_["mopac_run_params"]** ( string ): the control string to define the MOPAC job
                [default: "INDO C.I.=(6,3) CHARGE=0 RELSCF=0.000001 ALLVEC  WRTCONF=0.00  WRTCI=2"]
            * **params_["mopac_working_directory"]** ( string ) [ default: "mopac_wd"]
            * **params_["mopac_jobid"]** ( string ) [ default: "job_0000" ]
            * **params_["mopac_input_prefix"]** ( string ) [ default: "input_" ]
            * **params_["mopac_output_prefix"]** ( string ) [ default: "output_" ]

    Returns:
        None

    """

    params = dict(params_)

    critical_params = ["labels"]
    default_params = {"mopac_exe": "mopac",
                      "mopac_run_params": "INDO C.I.=(6,3) CHARGE=0 RELSCF=0.000001 ALLVEC  WRTCONF=0.00  WRTCI=2",
                      "mopac_working_directory": "mopac_wd",
                      "mopac_jobid": "job_0000",
                      "mopac_input_prefix": "input_", "mopac_output_prefix": "output_"
                      }
    comn.check_input(params, default_params, critical_params)

    labels = params["labels"]
    mopac_exe = params["mopac_exe"]
    mopac_run_params = params["mopac_run_params"]
    mopac_wd = params["mopac_working_directory"]
    mopac_jobid = params["mopac_jobid"]
    mopac_input_prefix = params["mopac_input_prefix"]
    mopac_output_prefix = params["mopac_output_prefix"]

    natoms = len(labels)
    ndof = 3 * natoms

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
    os.system(F"{mopac_exe} {mopac_input_filename} > {mopac_output_filename}")

    # Go back to the original directory
    os.chdir("../")


def make_ref(nelec, active_space=None):
    """
    Makes the reference determinant based on the number of electrons

    Args:
        nelec (int) : the total number of electrons in the system
        active_space (list of ints): the indices of allowed orbitals, starting from 1 [default: None]

    Returns:
        list of ints: representation of the reference (ground-state) determinant

        E.g. [1, -1, 2, -2] is the determinant of 4 electrons with 2 lowest orbitals doubly-filled
        Here, the indexing starts with 1, and negative values correspond to the beta-spin electron
        while positive values - to the alpha-spin electron.

        If the active_space is [2], then the above reference determinant will become just [2, -2]
    """

    res, orb = [], 1
    if active_space is None:
        for i in range(nelec):
            if i % 2 == 0:
                res.append(orb)
            else:
                res.append(-orb)
                orb = orb + 1

    else:
        for i in range(nelec):
            if i % 2 == 0:
                if orb in active_space:
                    res.append(orb)
            else:
                if orb in active_space:
                    res.append(-orb)
                orb = orb + 1
    return res


def make_alpha_excitation(ref_determinant, config):
    """
    This function creates an excitation of the alpha electron in the reference determinant
    exciting the electron from and to orbitals determined by the `config` argument

    Args:

        ref_determinant (list of ints): representation of the reference determinants

        config (list of 2 ints): [src, targt]  - source and target orbitals for single excitation, the
            indexing starts from 1, not from 0.

    Returns:
        list of ints: the representation of the excited configuration

    Example:
        make_alpha_excitation([1, -1, 2, -2], [2, 3])  corresponds to 2->3 excitation and should return

        [1, -1, 3, -2]
    """
    src = config[0]
    trgt = config[1]

    res = list(ref_determinant)
    nelec = len(res)
    for i in range(nelec):
        if res[i] == src:
            res[i] = trgt

    return res


def read_mopac_orbital_info(params_):
    """
    This function reads the MOs, configurations, and CI from the output files

    Args:

        params_ ( dict ): the dictionary containing key parameters

            * **params_["filename"]** ( string ) : the name of the file to read
            * **params_["orbital_space"]** (list of ints): the orbital numbers to be included, the indexing starts with 1, not 0 [default: None]
            * **params_["nstates"]** ( int ): the number of CI states + ground state to use, `nstates = 2` means 1 ground and 1 excited states

    Returns:
        (Es, MOs, E_CI, CI, configs):

            * Es - MATRIX(nact, nact): the matrix of the MO energies for active MOs
            * MOs - MATRIX(nao, nact): the matrix of MO-LCAO coefficients for active MOs
            * E_CI - MATRIX(nci, nci): the matrix of CI energies
            * CI - MATRIX(nconf, nci): the matrix of CI coefficients in the basis of spin-adapted configurations
            * sd_basis - (list of lists of `nact` ints): representation of the Slater determinants of `nelec`-electronic
                 states but with the restriction of using only the orbitals belonging to the active space. This is the reduced basis
                 The indexing starts from 1.
            * sd_basis_raw - (list of lists of `nelec` ints): representation of the Slater determinants of `nelec`-electronic
                 states. The indexing starts from 1. This is the original (raw) basis without the restriction.
            * actual_orbital_space - (list of ints): the actual active space considered. The indexing starts with 1, not 0.

    Notes:
            * nact - the number of active MOs to be included
            * nao - the number of AOs, it is the same as nmo
            * nmo - the number of MOs
            * nconf  - the number of configurations (spin-adapted Slater determinants)
            * nci - the number of CI states
            * nelec - the number of electrons


    """
    params = dict(params_)

    critical_params = []
    default_params = {"filename": "output", "orbital_space": None, "nstates":2 }
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

    # Determine the number of electrons:
    nocc = 0
    for i in range(nlines):
        line = output[i]
        if line.find("RHF CALCULATION, NO. OF DOUBLY OCCUPIED LEVELS") != -1:
            nocc = int(float(line.split()[8]))
    nelec = 2 * nocc

    # First, let's find where the MOs are and count how many of them we have
    ibeg, iend, nmo = 0, nlines - 1, 0
    for i in range(nlines):
        line = output[i]

        if line.find("MOLECULAR ORBITALS") != -1:
            ibeg = i
        if line.find("Reference determinate nber") != -1:
            iend = i

        if line.find("ROOT NO.") != -1:
            nmo = int(float(line.split()[-1]))

    if False:  # Make True for debugging
        print("nmo = ", nmo)
        print(output[ibeg:iend])

    nao = nmo
    actual_orbital_space = list(range(1, nmo + 1))
    if params["orbital_space"] is None:
        pass  # defualt - use all orbitals
    else:
        actual_orbital_space = list(params["orbital_space"])

    nact = len(actual_orbital_space)

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
    for j in range(nblocks - 1):
        i = break_lines[j]
        tmp = output[i + 2].split()
        for e in tmp:
            ener = float(e)
            E.append(ener)

    # Es = MATRIX(nmo, nmo)
    Es = MATRIX(nact, nact)
    for i in range(nact):
        j = actual_orbital_space[i] - 1
        Es.set(i, i, E[j])

    if False:  # Make True for debugging
        print(E)

    # Now read MOs:
    mo_indx = 0
    # MOs = MATRIX(nao, nmo)
    MOs = MATRIX(nao, nact)

    for j in range(nblocks - 1):
        i = break_lines[j]
        tmp = output[i + 2].split()
        ncols = len(tmp)
        ao_indx = 0
        for i in range(break_lines[j] + 5, break_lines[j + 1]):
            tmp = output[i].split()
            sz = len(tmp)
            if (sz == ncols + 3):
                for a in range(ncols):
                    coeff = float(tmp[3 + a])
                    # MOs.set(ao_indx, mo_indx + a, coeff)
                    if mo_indx + a + 1 in actual_orbital_space:
                        indx = actual_orbital_space.index(mo_indx + a + 1)
                        MOs.set(ao_indx, indx, coeff)
                ao_indx += 1
        mo_indx += ncols

    if False:  # Make True for debugging
        MOs.show_matrix("MO.txt")

    # Find the configurations
    nconfig = 0
    configs = []
    for i in range(nlines):
        line = output[i]
        if line.find("The lowest") != -1 and line.find("spin-adapted configurations of multiplicity=  1"):
            tmp = line.split()
            if (len(tmp) > 3):
                nconfig = int(float(tmp[2]))
                # ========= Found the configurations info, now read the configurations ========
                for iconfig in range(nconfig):
                    tmp = output[i + 4 + iconfig].split()
                    if len(tmp) == 12:
                        configs.append([1, 1])  # this is the ground state determinant
                    if len(tmp) == 13:
                        i_orb = int(float(tmp[10].split(")->(")[0]))
                        j_orb = int(float(tmp[11]))
                        configs.append([i_orb, j_orb])  # orbital indexing from 1 - as in MOPAC

    # Raw SD basis - using the indixes (starting with 1) used in the original N-electron system
    sd_basis_raw = []
    ref = make_ref(nelec, actual_orbital_space)  # make a restricted reference determinant
    for cnf in configs:
        sd = make_alpha_excitation(ref, cnf)
        sd_basis_raw.append(sd)

    # Make the SD basis restructed - so reindex everything according to active_space
    sd_basis = []
    for cnf_raw in sd_basis_raw:
        cnf = []
        for a in cnf_raw:
            sign = 1
            if a < 0.0:
                sign = -1
            b = actual_orbital_space.index(abs(a)) + 1
            cnf.append(sign * b)
        sd_basis.append(cnf)

    if False:  # Make True for debugging
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
    if nci < params["nstates"]:
        print("Not enough CI states requested in the INDO input line\nExiting...\n")
        sys.exit(0)
    else:
        nci = params["nstates"]


    if False:  # Make True for debugging
        print(F"The number of CI states = {nci}")
        for i in range(nci):
            print(F"CI block {i}")
            print(output[ci_beg[i]:ci_end[i]])

    # Now, read the information about CI states
    E_CI = MATRIX(nci, nci)
    CI = MATRIX(nconfig, nci)
    for i in range(nci):
        for j in range(ci_beg[i], ci_end[i]):
            tmp = output[j].split()
            sz = len(tmp)
            if sz == 7:
                ener = float(tmp[2]) * units.ev2au  # convert to a.u.
                E_CI.set(i, i, ener)
            elif sz == 4:
                if tmp[0] == "Config":
                    iconf = int(float(tmp[1])) - 1
                    coeff = float(tmp[2])
                    CI.set(iconf, i, coeff)


    return Es, MOs, E_CI, CI, sd_basis, sd_basis_raw, actual_orbital_space


def mopac_compute_adi(q, params, full_id):
    """

    This function creates an input for MOPAC, runs such calculations, extracts the current information on
    state energies and CI vectors, computes the required properties (such as time-overlaps, or NACs, etc.)
    and stores the current information in the "previous variables", so that we could compute the dependent properties
    on the next time-step

    Args:
        q ( MATRIX(ndof, ntraj) ): coordinates of the particle [ units: Bohr ]
        params ( list of dictionaries ): model parameters, for each trajectory; this parameters variable will be used to
            also store the previous calculations, but that has to be done separately for each trajectory, i = full_id[-1]. That's why we
            are making it into a list of dictionaries

            * **params[i]["orbital_space"]** (list of ints): the orbitals to read and use in MO overlaps, the indexing starts with 1, not 0 [default: None]
            * **params[i]["timestep"]** (int): the index of the timestep for trajectory i [ Required ]
            * **params[i]["is_first_time"]** (int): the flag indicating if this is the new calculation for this trajectory or not
              if it is True (1), the current values will be used as if they were previous; if False (0) - the previously stored values
              will be used [ default: True ]
            * **params[i]["labels"]** ( list of strings ): the labels of atomic symbolc - for all atoms,
                and in a order that is consistent with the coordinates (in triples) stored in `q`.
                The number of this labels is `natoms`, such that `ndof` = 3 * `natoms`. [ Required ]
            * **params[i]["mopac_exe"]** ( string ):  the full path to `the mopac` executable [ defaut: "mopac" ]
            * **params_[i]["mopac_run_params"]** ( string ): the control string to define the MOPAC job
                [default: "INDO C.I.=(6,3) CHARGE=0 RELSCF=0.000001 ALLVEC  WRTCONF=0.00  WRTCI=2"]
            * **params[i]["mopac_working_directory"]** ( string ) [ default: "mopac_wd"]
            * **params[i]["mopac_jobid"]** ( string ) [ default: F"timestep_{timestep}_traj_{i}" ]
            * **params[i]["mopac_input_prefix"]** ( string ) [ default: "input_" ]
            * **params[i]["mopac_output_prefix"]** ( string ) [ default: "output_" ]
            * **params[i]["dt"]** ( float ) - the time interval between the snapshots [ units: a.u.; default: 41 a.u. = 1 fs ]
            * **params[i]["do_Lowdin"]** ( bool or int): 0 - don't do - use the raw inputs [ defualt ]; 1 - do it - correct for rounding errors
            * **params[i]["CAS"]**  (list of of ints, int): the list represents the indices of spatial orbitals used in CAS definition in INDO-CI
                calculations, absoluted values; the "int" represents the number of electrons in the active space. For example, if the HOMO is 8
                and and CAS was (6, 3) - 6 orbitals with 3 doubly-filled orbitals, then one would use: [ [6,7,8,  9, 10, 11], 6]
            * **params[i]["is_single_excitation"]** ( bool or int): 0 - consider full configuration space, 1 - use single excitation specific methods
    Returns:
        PyObject: obj, with the members:

            * obj.ham_adi ( CMATRIX(nstates,nstates) ): adiabatic Hamiltonian, in the CI basis
            * obj.hvib_adi ( CMATRIX(nstates,nstates) ): adiabatic vibronic Hamiltonian, in the CI basis
            * obj.basis_transform ( CMATRIX(nstates,nstates) ): assumed the identity - don't use it yet!
            * obj.time_overlap_adi. ( CMATRIX(nstates,nstates) ): time-overlap in the CI basis

        Also, the following key-value pairs right in the input parameters dictionaries will be created/updated:

           * **params[i]["E_prev"]** (MATRIX(nact, nact)): the diagonal matrix of MO energies in the MO basis belonging to the active space
           * **params[i]["MO_prev"]** (MATRIX(nao, nact)): the MOs belonging to the active space
           * **params[i]["E_CI_prev"]** (MATRIX(nstates, nstates)): the diagonal matrix of the CI state energies
           * **params[i]["CI_prev"]** (MATRIX(nconf, nstates)): the CI eigenvectors in the basis of configuration functions
           * **params[i]["configs_prev"]** (list of lists): the list of reduced-notation configurations
           * **params[i]["configs_raw_prev"]** (list of lists): the list of full-notation configurations

    """

    # params = dict(params_)
    sqt2 = math.sqrt(2.0)

    Id = Cpp2Py(full_id)
    itraj = Id[-1]
    coords = q.col(itraj)

    critical_params = ["labels", "timestep"]
    default_params = {"mopac_exe": "mopac", "is_first_time": True, "orbital_space": None, "nstates":2,
                      "mopac_run_params": "INDO C.I.=(6,3) CHARGE=0 RELSCF=0.000001 ALLVEC  WRTCONF=0.00  WRTCI=2",
                      "mopac_working_directory": "mopac_wd",
                      "mopac_input_prefix": "input_", "mopac_output_prefix": "output_",
                      "dt": 1.0 * units.fs2au, "do_Lowdin": 0, 
                      "CAS":[ [1,2], 2], "mult_S":0, "mult_Ms":0
                      }
    comn.check_input(params[itraj], default_params, critical_params)

    timestep = params[itraj]["timestep"]
    labels = params[itraj]["labels"]
    is_first_time = params[itraj]["is_first_time"]
    mopac_exe = params[itraj]["mopac_exe"]
    mopac_run_params = params[itraj]["mopac_run_params"]
    mopac_wd = params[itraj]["mopac_working_directory"]
    mopac_jobid = params[itraj]["mopac_jobid"] = F"timestep_{timestep}_traj_{itraj}"
    mopac_input_prefix = params[itraj]["mopac_input_prefix"]
    mopac_output_prefix = params[itraj]["mopac_output_prefix"]
    orbital_space = params[itraj]["orbital_space"]
    nstates = params[itraj]["nstates"]
    dt = params[itraj]["dt"]
    do_Lowdin = params[itraj]["do_Lowdin"]
    CAS = params[itraj]["CAS"]
    mult_S = params[itraj]["mult_S"]
    mult_Ms = params[itraj]["mult_Ms"]
    is_single_excitation = params[itraj]["is_single_excitation"]
    natoms = len(labels)
    ndof = 3 * natoms

    # Run the calculations
    # print("================ RUN MOPAC =================\n")
    run_mopac(coords, params[itraj])

    # Read the current output
    filename = F"{mopac_wd}/{mopac_input_prefix}{mopac_jobid}.out"

    # print("================ READ MOPAC =================\n")
    # print("cwd = ", os.getcwd(), " reading the file ", filename)
    E_curr, MO_curr, E_CI_curr, CI_curr, _, configs_raw_curr, actual_orbital_space = read_mopac_orbital_info({
      "filename": filename, 
      "orbital_space":orbital_space,
      "nstates":nstates
    })
    if is_single_excitation:
        configs_curr, CSF_coeff_curr = interfaces.configs_and_T_matrix_singlet(configs_raw_curr, 
                                                                   active_space=CAS[0], 
                                                                   orbital_space = actual_orbital_space,
                                                                   nelec=CAS[1], S=mult_S, Ms=mult_Ms)
    else
        configs_curr, CSF_coeff_curr = interfaces.configs_and_T_matrix(configs_raw_curr, 
                                                                   active_space=CAS[0], 
                                                                   orbital_space = actual_orbital_space,
                                                                   nelec=CAS[1], S=mult_S, Ms=mult_Ms)
    # Get the properties at the previous time-step
    E_prev, MO_prev, E_CI_prev, CI_prev, configs_prev, configs_raw_prev = None, None, None, None, None, None
    CSF_coeff_prev = None

    # print("================ THE REST =================\n")
    if is_first_time:
        # On the first step, assume the current properties are as the previous
        E_prev, MO_prev = MATRIX(E_curr), MATRIX(MO_curr)
        E_CI_prev, CI_prev = MATRIX(E_CI_curr), MATRIX(CI_curr)
        CSF_coeff_prev = CMATRIX(CSF_coeff_curr)
        configs_prev = list(configs_curr)
    else:
        # Otherwise, retrieve the previously-stored data
        E_prev = params[itraj]["E_prev"]
        MO_prev = params[itraj]["MO_prev"]
        E_CI_prev = params[itraj]["E_CI_prev"]
        CI_prev = params[itraj]["CI_prev"]
        CSF_coeff_prev = params[itraj]["CSF_coeff_prev"]
        configs_prev = params[itraj]["configs_prev"]

    #nstates = min(CI_curr.num_of_cols, Nstates)
    # print("nstates = ", nstates)

    # Do the calculations - time-overlaps, energies, and Hvib
    obj = tmp()
    obj.ham_adi = CMATRIX(nstates, nstates)
    obj.nac_adi = CMATRIX(nstates, nstates)
    obj.hvib_adi = CMATRIX(nstates, nstates)
    obj.basis_transform = CMATRIX(nstates, nstates)
    obj.time_overlap_adi = CMATRIX(nstates, nstates)
    obj.overlap_adi = CMATRIX(nstates, nstates)
    # Don't do the derivatives yet
    # obj.d1ham_adi = CMATRIXList();
    # obj.dc1_adi = CMATRIXList();
    # for idof in range(ndof):
    #    obj.d1ham_adi.append( CMATRIX(nstates, nstates) )
    #    obj.dc1_adi.append( CMATRIX(nstates, nstates) )

    #======================= MO ===============================
    # MO overlaps - needed for Lowdin orthonormalization
    S_prev, S_curr, U_prev, U_curr = None, None, None, None
    
    if do_Lowdin:
        S_prev = MO_prev.T() * MO_prev
        im = MATRIX(S_prev);  im *= 0.0
        S_prev = CMATRIX(S_prev, im)

        S_curr = MO_curr.T() * MO_curr
        im = MATRIX(S_curr);  im *= 0.0
        S_curr = CMATRIX(S_curr, im)

        U_prev = step3.get_Lowdin_general(S_prev).real()
        U_curr = step3.get_Lowdin_general(S_curr).real()
    
    # Time-overlap in the MO basis
    mo_st = MO_prev.T() * MO_curr
    if do_Lowdin:
        mo_st = U_prev.T() * mo_st * U_curr  # Lowdin correction on MOs

    # Overlaps in the SD basis - for the Lowdin:
    ci_ovlp_curr = None
    if do_Lowdin:
        ident_curr = U_curr.T() * S_curr.real() * U_curr
        # The original (older) approach
        #ovlp_sd_curr = mapping2.ovlp_mat_arb(configs_curr, configs_curr, ident_curr, False).real()

        # The new way:
        # I don't like this way - we need to fix it later
        ident_curr = data_conv.MATRIX2nparray(ident_curr).real  # MATRIX -> real np array
        ovlp_sd_curr = mapping3.ovlp_mat_arb(configs_curr, configs_curr, ident_curr, actual_orbital_space)  # do the calculations
        ovlp_sd_curr = data_conv.nparray2MATRIX( ovlp_sd_curr )  # real np array -> MATRIX

        # 10/31/2025 - no need to scale the SD overlaps
        #ovlp_sd_curr.scale(-1, 0, sqt2)
        #ovlp_sd_curr.scale(0, -1, sqt2)
        #ovlp_sd_curr.scale(0, 0, 0.5)
 
        # CSF:
        T_curr = MATRIX(CSF_coeff_curr.real())    
        ovlp_csf_curr = T_curr.T() * ovlp_sd_curr * T_curr


        ident_prev = U_prev.T() * S_prev.real() * U_prev
        #ovlp_sd_prev = mapping2.ovlp_mat_arb(configs_prev, configs_prev, ident_prev, False).real()

        ident_prev = data_conv.MATRIX2nparray(ident_prev).real  # MATRIX -> real np array
        ovlp_sd_prev = mapping3.ovlp_mat_arb(configs_prev, configs_prev, ident_prev, actual_orbital_space)  # do the calculations
        ovlp_sd_prev = data_conv.nparray2MATRIX( ovlp_sd_prev )  # complex np array -> CMATRIX

        #ovlp_sd_prev.scale(-1, 0, sqt2)
        #ovlp_sd_prev.scale(0, -1, sqt2)
        #ovlp_sd_prev.scale(0, 0, 0.5)

        # CSF:
        T_prev = MATRIX(CSF_coeff_prev.real())
        ovlp_csf_prev = T_prev.T() * ovlp_sd_prev * T_prev


        ovlp_ci_curr = CI_curr.T() * ovlp_csf_curr * CI_curr
        im = MATRIX(ovlp_ci_curr);  im *= 0.0
        ovlp_ci_curr = CMATRIX(ovlp_ci_curr, im)

        ovlp_ci_prev = CI_prev.T() * ovlp_csf_prev * CI_prev
        im = MATRIX(ovlp_ci_prev);  im *= 0.0
        ovlp_ci_prev = CMATRIX(ovlp_ci_prev, im)

        U_prev = step3.get_Lowdin_general(ovlp_ci_prev).real()
        U_curr = step3.get_Lowdin_general(ovlp_ci_curr).real()

        # Now corrected CI states overlap
        ci_ovlp_curr = U_curr.T() * ovlp_ci_curr.real() * U_curr

    else:
        # Overlap in MO basis:
        S_curr = MO_curr.T() * MO_curr

        # Overlap in SD basis:
        S_curr_np = data_conv.MATRIX2nparray(S_curr).real  # MATRIX -> real np array
        ovlp_sd_curr_np = sd.slater_overlap_matrix(configs_curr, configs_curr, S_curr_np)  # do the calculations
        ovlp_sd_curr = data_conv.nparray2MATRIX( ovlp_sd_curr_np )  # complex np array -> MATRIX

        # CSF:
        T_curr = MATRIX(CSF_coeff_curr.real())
        ovlp_csf_curr = T_curr.T() * ovlp_sd_curr * T_curr
        ci_ovlp_curr = CI_curr.T() * ovlp_csf_curr * CI_curr


    #========== Time-overlap in the SD basis ===============
    # First, let's convert the MO overlaps to numpy format
    mo_st_np = data_conv.MATRIX2nparray(mo_st).real  # MATRIX -> real np array
    # Next, compute the matrix of SD overlaps
    time_ovlp_sd_np = sd.slater_overlap_matrix(configs_prev, configs_curr, mo_st_np)
    # Convert them back
    time_ovlp_sd = data_conv.nparray2MATRIX( time_ovlp_sd_np )  # complex np array -> MATRIX


    #=========== Time-overlap in the CSF basis =============
    T_prev = MATRIX(CSF_coeff_prev.real())
    T_curr = MATRIX(CSF_coeff_curr.real())
    time_ovlp_csf = T_prev.T() * time_ovlp_sd * T_curr

    #============ Time-overlap in the CI basis ============
    time_ovlp_ci = CI_prev.T() * time_ovlp_csf * CI_curr
    if do_Lowdin:
        time_ovlp_ci = U_prev.T() * time_ovlp_ci * U_curr  # Lowdin correction on CIs

    # Now, populate the allocated matrices
    for istate in range(nstates):
        energ = 0.5 * (E_CI_curr.get(istate, istate) + E_CI_prev.get(istate, istate))
        obj.ham_adi.set(istate, istate, energ * (1.0 + 0.0j))
        obj.hvib_adi.set(istate, istate, energ * (1.0 + 0.0j))
        obj.basis_transform.set(istate, istate, 1.0 + 0.0j)  # assume identity

        for jstate in range(nstates):
            obj.time_overlap_adi.set(istate, jstate, time_ovlp_ci.get(istate, jstate) * (1.0 + 0.0j))
            if ci_ovlp_curr is not None:
                obj.overlap_adi.set(istate, jstate, ci_ovlp_curr.get(istate, jstate) * (1.0 + 0.0j))
        # Don't do this yet
        # for idof in range(ndof):
        #    obj.d1ham_adi[idof].set(istate, istate, grad.get(idof, 0) * (1.0+0.0j) )

    # Update the Hvib:
    for istate in range(nstates):
        for jstate in range(istate + 1, nstates):
            dij = (obj.time_overlap_adi.get(istate, jstate) - obj.time_overlap_adi.get(jstate, istate)) / (2.0 * dt)
            obj.hvib_adi.set(istate, jstate, dij * (0.0 - 1.0j))
            obj.hvib_adi.set(jstate, istate, dij * (0.0 + 1.0j))

    # Now, make the current the previous and reset the flag `is_first_time` to False
    # Note - we directly modify the input parameters
    params[itraj]["E_prev"] = E_curr
    params[itraj]["MO_prev"] = MO_curr
    params[itraj]["E_CI_prev"] = E_CI_curr
    params[itraj]["CI_prev"] = CI_curr
    params[itraj]["configs_prev"] = configs_curr
    params[itraj]["CSF_coeff_prev"] = CSF_coeff_curr
    params[itraj]["is_first_time"] = False

    return obj
