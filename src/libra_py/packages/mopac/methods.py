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
import copy
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
import libra_py.citools.ci as ci

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
        (info, MOs, data):

            * info (dict): the key numbers (see the notes below)
            * MOs - (np.ndarray; shape = [nao, nact]): the matrix of MO-LCAO coefficients for active MOs
            * data ( list of 3 elements):
              - data[0] - list of `nstates-1` excited (not including the ground state) state energies (a.u.)
              - data[1] - list of `nstates-1` lists, each containing lists of pairs [i, j] denoting the single excitations 
                 of the i->j kind entering the expression of the corresponding multiconfigurational state
              - data[2] - list of the amplitudes of single-particle excitations entering all the excited states, isomorphic
                 to `data[1]`

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
    min_occ, max_occ = 0, 0
    min_vir, max_vir = 0, 0
    for i in range(nlines):
        line = output[i]
        if line.find("RHF CALCULATION, NO. OF DOUBLY OCCUPIED LEVELS") != -1:
            nocc = int(float(line.split()[8]))
        if line.find("SINGLE excitations FROM orbs") != -1:
            min_occ = int(float(line.split()[4]))
            max_occ = int(float(line.split()[6]))
            min_vir = int(float(line.split()[9]))
            max_vir = int(float(line.split()[11]))
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

    #================ Now read MOs: ====================
    mo_indx = 0
    MOs = np.zeros( (nao, nact), dtype=np.float64)

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
                    if mo_indx + a + 1 in actual_orbital_space:
                        indx = actual_orbital_space.index(mo_indx + a + 1)
                        MOs[ao_indx, indx] = coeff
                ao_indx += 1
        mo_indx += ncols


    # Find the configurations
    nconfig = 0
    configs_dict = {}
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
                        configs_dict[1] = [1, 1]
                    if len(tmp) == 13:
                        i_orb = int(float(tmp[10].split(")->(")[0]))
                        j_orb = int(float(tmp[11]))
                        conf_indx = int(float(tmp[0]))
                        configs_dict[conf_indx] = [i_orb, j_orb]


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
    confs = []
    E_CI = []
    CI = []
    for i in range(1, nci):
        ci_i = []
        conf_i = []
        for j in range(ci_beg[i], ci_end[i]):
            tmp = output[j].split()
            sz = len(tmp)
            if sz == 7:
                ener = float(tmp[2]) * units.ev2au  # convert to a.u.
                E_CI.append(ener)
            elif sz == 4:
                if tmp[0] == "Config":
                    iconf = int(float(tmp[1])) # - 1
                    coeff = float(tmp[2])
                    if iconf > 1:
                        ci_i.append(coeff)
                        conf_i.append( configs_dict[ iconf ] )
                    
        CI.append(ci_i)
        confs.append(conf_i)

    data = [E_CI, confs, CI ]

    info = { "nocc":nocc, 
             "nelec":nelec,
             "nao":nao, "nmo":nmo,
             "nci":nci,
             "min_occ":min_occ, "max_occ":max_occ, 
             "min_vir":min_vir, "max_vir":max_vir,
             "nact":nact, "actual_orbital_space": list(actual_orbital_space)
           }


    return info, MOs, data


def mopac_compute_adi(q, params, full_id):
    """
    Run MOPAC calculations for a given trajectory, extract electronic structure
    information, compute time-dependent electronic properties (energies,
    overlaps, NACs), and store the current results for use in the next time step.

    This function is **trajectory-aware**: all parameters and cached data are
    stored separately for each trajectory and accessed using the trajectory
    index `itraj = full_id[-1]`.

    Parameters
    ----------
    q : MATRIX(ndof, ntraj)
        Nuclear coordinates for all trajectories [Bohr].
        Each column corresponds to a single trajectory.

    params : list of dict
        Per-trajectory parameter dictionaries.
        `params[i]` contains **both input settings and cached data**
        for trajectory `i`.

        This function **modifies** `params[i]` in place by storing
        quantities from the current time step for use at the next step.

        Required keys in `params[i]`:
            * labels : list[str]
                Atomic symbols, length = natoms
            * timestep : int
                Time step index for this trajectory

        Optional keys in `params[i]` (defaults applied if missing):
            * mopac_exe : str
                Path to MOPAC executable [default: "mopac"]
            * mopac_run_params : str
                MOPAC control string
            * mopac_working_directory : str
                Directory for MOPAC I/O [default: "mopac_wd"]
            * mopac_input_prefix : str
            * mopac_output_prefix : str
            * dt : float
                Time step [a.u.]
            * orbital_space : list[int] or None
            * nstates : int
                Number of electronic states (including ground)
            * CAS : list
                CAS definition for INDO-CI
            * mult_S, mult_Ms : int
                Spin multiplicity parameters
            * is_singlet_excitation : bool or int
            * is_first_time : bool
                If True, current data are treated as previous (first step)

        Cached quantities created/updated in `params[i]`:
            * MO_prev : ndarray
                Active-space molecular orbitals from previous step
            * data_prev : tuple
                CI energies, configurations, and amplitudes from previous step
            * is_first_time : bool
                Set to False after first invocation

    full_id : list[int]
        Full trajectory identifier; the trajectory index is taken as
        `itraj = full_id[-1]`.

    Returns
    -------
    obj : PyObject
        Container with computed adiabatic quantities:
            * obj.ham_adi : CMATRIX(nstates, nstates)
                Adiabatic electronic Hamiltonian
            * obj.hvib_adi : CMATRIX(nstates, nstates)
                Adiabatic vibronic Hamiltonian
            * obj.time_overlap_adi : CMATRIX(nstates, nstates)
                CI time-overlap matrix
            * obj.basis_transform : CMATRIX(nstates, nstates)
                Basis transformation (currently identity)

    """

    Id = Cpp2Py(full_id)
    itraj = Id[-1]

    # Sanity check on params structure
    if not isinstance(params, (list, tuple)):
        raise TypeError(
            "params must be a list (or tuple) of per-trajectory dictionaries; "
            f"got {type(params)}"
        )

    if itraj >= len(params):
        raise IndexError(
            f"Trajectory index itraj={itraj} out of range for params (len={len(params)})"
        )

    if not isinstance(params[itraj], dict):
        raise TypeError(
            f"params[{itraj}] must be a dictionary; got {type(params[itraj])}"
        )

    coords = q.col(itraj)

    critical_params = ["labels", "timestep"]
    default_params = {"mopac_exe": "mopac", "is_first_time": True, "orbital_space": None, "active_space": None,
                      "nstates":2,
                      "mopac_run_params": "INDO C.I.=(6,3) CHARGE=0 RELSCF=0.000001 ALLVEC  WRTCONF=0.00  WRTCI=2",
                      "mopac_working_directory": "mopac_wd",
                      "mopac_input_prefix": "input_", "mopac_output_prefix": "output_",
                      "dt": 1.0 * units.fs2au, "do_Lowdin": 0, 
                      "CAS":[ [1,2], 2], "mult_S":0, "mult_Ms":0, "is_singlet_excitation":0
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
    active_space = params[itraj]["active_space"]
    nstates = params[itraj]["nstates"]
    dt = params[itraj]["dt"]
    do_Lowdin = params[itraj]["do_Lowdin"]
    CAS = params[itraj]["CAS"]
    mult_S = params[itraj]["mult_S"]
    mult_Ms = params[itraj]["mult_Ms"]
    is_singlet_excitation = params[itraj]["is_singlet_excitation"]
    natoms = len(labels)
    ndof = 3 * natoms

    # Run the calculations
    # print("================ RUN MOPAC =================\n")
    run_mopac(coords, params[itraj])

    # Read the MOPAC output
    # print("================ READ MOPAC =================\n")
    filename = F"{mopac_wd}/{mopac_input_prefix}{mopac_jobid}.out"
    read_params = {"nstates":nstates, "filename":filename, "orbital_space":None}
    info, MO_curr, data_curr = read_mopac_orbital_info(read_params)

    # Get the properties at the previous time-steps
    MO_prev, data_prev = None, None
    if is_first_time:
        # On the first step, assume the current properties are as the previous
        MO_prev = copy.deepcopy(MO_curr)
        data_prev = copy.deepcopy(data_curr)
    else:
        # Otherwise, retrieve the previously-stored data
        MO_prev = copy.deepcopy(params[itraj]["MO_prev"])
        data_prev = copy.deepcopy(params[itraj]["data_prev"])

    # Do the calculations - time-overlaps, energies, and Hvib
    obj = tmp()
    obj.ham_adi = CMATRIX(nstates, nstates)
    obj.nac_adi = CMATRIX(nstates, nstates)
    obj.hvib_adi = CMATRIX(nstates, nstates)
    obj.basis_transform = CMATRIX(nstates, nstates)
    obj.time_overlap_adi = CMATRIX(nstates, nstates)
    obj.overlap_adi = CMATRIX(nstates, nstates)

    #======================= MO ===============================
    # MO overlaps
    st_mo_orb = MO_prev.T @ MO_curr

    # Make it doubled - block-matrix
    st_mo = np.kron(np.eye(2), st_mo_orb)

    #================= Compute CI time-overlaps =============
    ovlp_params = {"homo_indx":info["nocc"], 
                   "nocc":info["nocc"] - 1, 
                   "nvirt":info["nmo"] - info["nocc"], 
                   "nelec":info["nelec"], "nstates":nstates,
                   "active_space":active_space
                   }
    st_ci = ci.overlap(st_mo, data_prev, data_curr, ovlp_params)

    #=============== Now, populate the allocated matrices ======================
    for istate in range(nstates):
        energ = 0.0
        if istate > 0:
            energ = float(0.5 * (data_prev[0][istate-1] + data_curr[0][istate-1]))

        obj.ham_adi.set(istate, istate, energ * (1.0 + 0.0j))
        obj.hvib_adi.set(istate, istate, energ * (1.0 + 0.0j))
        obj.basis_transform.set(istate, istate, 1.0 + 0.0j)  # assume identity

        for jstate in range(nstates):
            obj.time_overlap_adi.set(istate, jstate, float(st_ci[istate, jstate]) * (1.0 + 0.0j))

    # Update the Hvib:
    for istate in range(nstates):
        for jstate in range(istate + 1, nstates):
            dij = (obj.time_overlap_adi.get(istate, jstate) - obj.time_overlap_adi.get(jstate, istate)) / (2.0 * dt)
            obj.hvib_adi.set(istate, jstate, dij * (0.0 - 1.0j))
            obj.hvib_adi.set(jstate, istate, dij * (0.0 + 1.0j))

    # Now, make the current the previous and reset the flag `is_first_time` to False
    # Note - we directly modify the input parameters
    params[itraj]["MO_prev"] = copy.deepcopy(MO_curr)
    params[itraj]["data_prev"] = copy.deepcopy(data_curr)
    params[itraj]["is_first_time"] = False

    return obj
