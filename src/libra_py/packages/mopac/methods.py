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


import os, sys, math, copy, re, subprocess
import numpy as np

from liblibra_core import *
import util.libutil as comn

from libra_py import units
from libra_py import scan
from libra_py import regexlib as rgl
from libra_py import data_conv

import libra_py.packages.cp2k.methods as CP2K_methods
import libra_py.workflows.nbra.step3 as step3

import libra_py.citools.slatdet as sd
import libra_py.citools.interfaces as interfaces
import libra_py.citools.ci as ci
import libra_py.orthogonalizations as ortho

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


def make_open_shell_reference(n_doubly, multiplicity):
    """Construct the highest-Ms restricted open-shell reference determinant."""
    return interfaces.make_open_shell_reference(n_doubly, multiplicity)


def spin_quantum_numbers(multiplicity, spin_projection=None):
    """Validate a multiplicity and return its total spin and selected Ms."""
    return interfaces.spin_quantum_numbers(multiplicity, spin_projection)


def add_mopac_spin_keyword(run_params, multiplicity):
    """Add the requested MOPAC spin keyword and reject conflicting keywords."""
    spin_names = {
        1: "SINGLET", 2: "DOUBLET", 3: "TRIPLET", 4: "QUARTET",
        5: "QUINTET", 6: "SEXTET", 7: "SEPTET", 8: "OCTET", 9: "NONET",
    }
    multiplicity = int(multiplicity)
    tokens = run_params.upper().split()
    name_to_multiplicity = {name: value for value, name in spin_names.items()}
    for token in tokens:
        selected = None
        if token in name_to_multiplicity:
            selected = name_to_multiplicity[token]
        elif token.startswith("MS="):
            selected = int(round(2 * float(token.split("=", 1)[1]) + 1))
        if selected is not None:
            if selected != multiplicity:
                raise ValueError(
                    f"mopac_run_params selects multiplicity {selected}, but "
                    f"multiplicity={multiplicity} was requested"
                )
            return run_params
    keyword = spin_names.get(multiplicity, f"MS={0.5 * (multiplicity - 1):g}")
    return f"{run_params} {keyword}"


def run_mopac(coords, params_):
    """
    Execute a MOPAC calculation in a thread-safe manner within a specified
    working directory.

    This function prepares a MOPAC input file, executes the MOPAC quantum
    chemistry package, and stores all generated files in a dedicated working
    directory. It is designed for use in parallel workflows (e.g., multiple
    trajectories in nonadiabatic dynamics simulations), where each calculation
    runs independently without changing the global working directory.

    The function performs the following steps:

        1. Ensures the working directory exists.
        2. Creates a MOPAC input file containing the molecular geometry and
           requested calculation settings.
        3. Executes the MOPAC binary within the working directory.
        4. Redirects the program output to a designated output file.
        5. Captures errors and raises an exception if the calculation fails.

    Parameters
    ----------
    coords : MATRIX(ndof, 1)
        Cartesian atomic coordinates of all atoms, in atomic units (Bohr).
        Here, ndof = 3 * natoms, where natoms is the number of atoms.
        Coordinates are ordered as

            [x1, y1, z1, x2, y2, z2, ..., xN, yN, zN]^T

        and stored as a Libra ``MATRIX`` object.

    params_ : dict
        Dictionary containing MOPAC execution parameters.

        Required keys
        -------------
        labels : list of str
            Atomic symbols (length N), e.g., ["O", "H", "H"].

        Optional keys
        -------------
        mopac_exe : str, optional
            Name or full path of the MOPAC executable.

            Default: "mopac"

        mopac_run_params : str, optional
            MOPAC keyword string defining the calculation type and settings.

            Example:
            "INDO C.I.=(6,3) CHARGE=0 RELSCF=0.000001 ALLVEC
             WRTCONF=0.00 WRTCI=2"

            Default:
            "INDO C.I.=(6,3) CHARGE=0 RELSCF=0.000001 ALLVEC
             WRTCONF=0.00 WRTCI=2"

        mopac_working_directory : str, optional
            Directory in which the calculation will be executed.
            All input and output files are created in this directory.

            Default: "mopac_wd"

        mopac_jobid : str, optional
            Identifier appended to input and output file names.

            Default: "job_0000"

        mopac_input_prefix : str, optional
            Prefix used when generating the MOPAC input file.

            Default: "input_"

        mopac_output_prefix : str, optional
            Prefix used when generating the MOPAC output file.

            Default: "output_"

    Returns
    -------
    None
        The function executes MOPAC and produces output files in the
        working directory. No value is returned.

    Raises
    ------
    FileNotFoundError
        If the specified MOPAC executable cannot be found.

    subprocess.CalledProcessError
        If the MOPAC calculation terminates with a non-zero exit code.
        Standard output and error streams are attached to the exception
        object and can be inspected for debugging.

    Side Effects
    ------------
    - Creates or modifies files in the working directory, including:
        - MOPAC input file (*.mop)
        - Output log file
        - Auxiliary files generated by MOPAC
          (*.out, *.arc, *.aux, *.den, etc., depending on settings)

    - Executes the MOPAC code.

    Notes
    -----
    - The function avoids the use of ``os.chdir`` and instead relies on
      the ``cwd`` argument of ``subprocess.run``. This makes the function
      safe for concurrent execution in multithreaded or multiprocess
      workflows.

    - Each calculation should use a unique working directory to avoid
      file collisions.

    - Coordinates are assumed to be provided in atomic units (Bohr).

    Examples
    --------
    >>> params = {
    >>>     "labels": ["O", "H", "H"],
    >>>     "mopac_working_directory": "wd_itraj0",
    >>>     "mopac_jobid": "traj0",
    >>>     "mopac_run_params":
    >>>         "INDO C.I.=(6,3) CHARGE=0 RELSCF=0.000001 "
    >>>         "ALLVEC WRTCONF=0.00 WRTCI=2"
    >>> }
    >>>
    >>> run_mopac(coords, params)

    """
    
    params = dict(params_)

    labels = params["atom_labels"]
    exe = params.get("exe", "mopac")
    mopac_run_params = params.get("mopac_run_params", "INDO C.I.=(6,3) CHARGE=0 RELSCF=0.000001 ALLVEC WRTCONF=0.00 WRTCI=2" )
    mopac_wd = params.get("working_directory", "mopac_wd")
    mopac_jobid = params.get("mopac_jobid", "job_0000")
    mopac_input_prefix = params.get("mopac_input_prefix", "input_")
    mopac_output_prefix = params.get("mopac_output_prefix", "output_")

    
    # Create working directory
    os.makedirs(mopac_wd, exist_ok=True)

    # File names
    mopac_input_filename = f"{mopac_input_prefix}{mopac_jobid}"
    mopac_output_filename = f"{mopac_output_prefix}{mopac_jobid}"

    # Full paths
    mopac_input_path = os.path.join(mopac_wd, mopac_input_filename)
    mopac_output_path = os.path.join(mopac_wd, mopac_output_filename)

    # Create input
    make_mopac_input(
        mopac_input_path,
        mopac_run_params,
        labels,
        coords
    )

    # Run MOPAC
    with open(mopac_output_path, "w") as fout:
        subprocess.run(
            [exe, mopac_input_filename],
            cwd=mopac_wd,
            check=True,
            stdout=fout,
            stderr=subprocess.STDOUT,
        )

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
    default_params = {
        "filename": "output", "orbital_space": None, "nstates": 2,
        "multiplicity": 1, "spin_projection": None,
    }
    comn.check_input(params, default_params, critical_params)

    out_file = params["filename"]
    multiplicity = int(params["multiplicity"])
    spin, spin_projection = spin_quantum_numbers(
        multiplicity, params["spin_projection"]
    )

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
    reference_det = make_open_shell_reference(nocc, multiplicity)
    nelec = len(reference_det)
    homo_indx = max(abs(orb) for orb in reference_det)

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
        match = re.search(r"spin-adapted configurations of multiplicity=\s*(\d+)", line)
        if line.find("The lowest") != -1 and match and int(match.group(1)) == multiplicity:
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
                        # Excitation into a singly occupied alpha orbital can
                        # only occur in the beta channel. Other one-
                        # determinant configurations use the alpha channel;
                        # spin adaptation supplies their partners as needed.
                        if j_orb in reference_det and -j_orb not in reference_det:
                            configs_dict[conf_indx] = [-i_orb, -j_orb]
                        else:
                            configs_dict[conf_indx] = [i_orb, j_orb]

    if not configs_dict:
        raise RuntimeError(
            f"Could not find MOPAC spin-adapted configurations for multiplicity {multiplicity}"
        )


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
             "homo_indx": homo_indx,
             "multiplicity": multiplicity,
             "spin": spin,
             "spin_projection": spin_projection,
             "reference_det": reference_det,
             "nao":nao, "nmo":nmo,
             "nci":nci,
             "min_occ":min_occ, "max_occ":max_occ, 
             "min_vir":min_vir, "max_vir":max_vir,
             "nact":nact, "actual_orbital_space": list(actual_orbital_space)
           }


    return info, MOs, data



def mopac_compute_adi(q, params, full_id):
    """
    Compute adiabatic-state energies, overlaps, time-overlaps, and vibronic
    Hamiltonian matrix elements using MOPAC electronic structure calculations.

    This function serves as a Libra-compatible electronic structure interface
    for nonadiabatic molecular dynamics simulations. For a given trajectory,
    it executes a MOPAC calculation, reads molecular orbital and configuration
    interaction (CI) information, computes adiabatic-state energies and
    time-overlaps between consecutive time steps, and constructs the vibronic
    Hamiltonian in the adiabatic representation.

    The function maintains trajectory-specific electronic structure data
    (molecular orbitals, CI vectors, and overlap matrices) between successive
    calls through the ``params`` dictionary, enabling the evaluation of
    nonadiabatic couplings via finite differences of wavefunction overlaps.

    Parameters
    ----------
    q : MATRIX(ndof, ntraj)
        Nuclear coordinates for all trajectories, stored as a Libra
        ``MATRIX`` object. Here, ``ndof = 3 * natoms`` and ``ntraj`` is
        the number of trajectories. Coordinates are assumed to be in
        atomic units (Bohr).

    params : dict
        Dictionary containing simulation and MOPAC parameters.

        Required keys
        -------------
        atom_labels : list of str
            Atomic symbols corresponding to the molecular geometry.

        Optional keys
        -------------
        timestep : int
            Current simulation time step.

            Default: 0

        energy_zero : float
            Energy shift applied to the computed electronic energies.

            Default: 0.0

        orbital_space : list of int or None
            User-defined orbital space to use in the electronic structure
            calculations.

            Default: None

        nstates : int
            Number of electronic states to include in the adiabatic basis.

            Default: 2

        dt : float
            Nuclear time step in atomic units.

            Default: 1.0 * units.fs2au

        working_directory_prefix : str
            Prefix used when creating trajectory-specific working
            directories.

            Default: "wd"

        mopac_input_prefix : str
            Prefix used for generated MOPAC input files.

            Default: "input_"

        mopac_output_prefix : str
            Prefix used for generated MOPAC output files.

            Default: "output_"

        mopac_run_params : str
            MOPAC keyword string defining the electronic structure
            calculation.

            Default:
            "INDO C.I.=(6,3) CHARGE=0 RELSCF=0.000001
             ALLVEC WRTCONF=0.00 WRTCI=2"

        do_Lowdin : bool
            Whether to apply Löwdin orthogonalization to the CI overlap
            matrices.

            Default: True

        nelec_act_space : int or None
            Number of active-space electrons used when constructing the
            determinant basis for overlap calculations. If ``None``,
            the active space reported by MOPAC is used.

            Default: None

        multiplicity : int
            Spin multiplicity ``2*S+1`` of the MOPAC CI calculation.
            Default: 1.

        spin_projection : int or float or None
            Desired spin projection ``Ms``. If omitted, the highest-weight
            component ``Ms=S`` is used.

        MO_prev : dict
            Trajectory-indexed storage of molecular orbital coefficients
            from the previous time step.

        data_prev : dict
            Trajectory-indexed storage of CI expansion data from the
            previous time step.

        s_ci_inv_prev : dict
            Trajectory-indexed storage of inverse square roots of CI
            overlap matrices used for Löwdin orthogonalization.

        is_first_time : dict
            Trajectory-indexed flags indicating whether the current call
            corresponds to the first simulation step.

        act_state : dict
            Trajectory-indexed active electronic state indices.

    full_id : intList
        Libra trajectory identifier. The last element is interpreted as
        the trajectory index.

    Returns
    -------
    obj : tmp
        Object containing electronic structure quantities in the adiabatic
        representation. The returned object contains the following members:

        ham_adi : CMATRIX(nstates, nstates)
            Adiabatic Hamiltonian matrix.

        hvib_adi : CMATRIX(nstates, nstates)
            Vibronic Hamiltonian matrix.

        nac_adi : CMATRIX(nstates, nstates)
            Nonadiabatic coupling matrix. Currently allocated but not
            explicitly populated.

        basis_transform : CMATRIX(nstates, nstates)
            Basis transformation matrix. Currently assumed to be the
            identity matrix.

        overlap_adi : CMATRIX(nstates, nstates)
            Adiabatic-state overlap matrix at the current time step.

        time_overlap_adi : CMATRIX(nstates, nstates)
            Time-overlap matrix between electronic states at consecutive
            time steps.

    Raises
    ------
    ValueError
        If the requested active space contains fewer electrons than
        required by the occupied orbitals reported by MOPAC.

    subprocess.CalledProcessError
        If the underlying MOPAC calculation fails.

    Side Effects
    ------------
    - Executes a MOPAC electronic structure calculation.
    - Creates or updates trajectory-specific working directories.
    - Modifies the input ``params`` dictionary by updating:

        * ``MO_prev``
        * ``data_prev``
        * ``s_ci_inv_prev``
        * ``is_first_time``

    Notes
    -----
    - The electronic energies are taken directly from the current MOPAC
      calculation and assigned to the diagonal elements of the adiabatic
      Hamiltonian.

    - Nonadiabatic couplings are computed from antisymmetrized
      time-overlaps using

        dij = (Sij(t,t+dt) - Sji(t,t+dt)) / (2*dt)

      and incorporated into the off-diagonal elements of the vibronic
      Hamiltonian.

    - For the first time step of a trajectory, the current electronic
      structure information is reused as the previous-step data so that
      overlap calculations remain well defined.

    - Trajectory-specific working directories make the function suitable
      for concurrent execution in ensemble and surface-hopping
      simulations.

    Examples
    --------
    >>> obj = mopac_compute_adi(q, params, full_id)
    >>> E0 = obj.ham_adi.get(0, 0).real
    >>> S01 = obj.time_overlap_adi.get(0, 1)
    >>> Hvib = obj.hvib_adi
    """

    # ================= Decode trajectory index =================
    Id = Cpp2Py(full_id)
    itraj = Id[-1]
    
    # ================= Extract coordinates =================
    coords = q.col(itraj)

    ndof = coords.num_of_rows
    nat = ndof // 3
    
    # ================= Safe param access =================
    params.setdefault("MO_prev", {})
    params.setdefault("data_prev", {})
    params.setdefault("s_ci_inv_prev", {})
    params.setdefault("is_first_time", {})
    params.setdefault("act_state", {})

    # ================= Read parameters =================
    # General: trajectory-agnostic
    atom_labels = params["atom_labels"]
    timestep = params.get("timestep", 0)
    energy_zero = params.get("energy_zero", 0.0 )
    orbital_space = params.get("orbital_space", None)
    nstates = params.get("nstates", 2)
    dt = params.get("dt", 1.0 * units.fs2au)
    wd_prefix = params.get("working_directory_prefix", "wd")
    mopac_input_prefix = params.get("mopac_input_prefix", "input_")
    mopac_output_prefix = params.get("mopac_output_prefix", "output_")
    mopac_run_params = params.get("mopac_run_params",
                                  "INDO C.I.=(6,3) CHARGE=0 RELSCF=0.000001 ALLVEC  WRTCONF=0.00  WRTCI=2")
    do_Lowdin = params.get("do_Lowdin", True)
    nelec_act_space = params.get("nelec_act_space", None)
    multiplicity = int(params.get("multiplicity", 1))
    spin, spin_projection = spin_quantum_numbers(
        multiplicity, params.get("spin_projection", None)
    )
    

    # Trajectory-specific
    is_first_time = params["is_first_time"].get(itraj, True)
    act_state = params["act_state"].get(itraj, 0)
    
    wd = f"{wd_prefix}_itraj{itraj}"
    
    # ================= Run MOPAC =================
    mopac_params = add_mopac_spin_keyword(
        copy.deepcopy(mopac_run_params), multiplicity
    )
    #mopac_params["StateOfInterest"] = act_state

    mopac_jobid = F"_timestep_{timestep}_traj_{itraj}"
    prms1 = {
        "atom_labels": atom_labels,
        "exe": params.get("exe", "mopac"),
        "mopac_run_params": mopac_params,
        "working_directory": wd,
        "mopac_jobid" : mopac_jobid,
        "mopac_input_prefix" : mopac_input_prefix,
        "mopac_output_prefix" : mopac_output_prefix
    }
    
    run_mopac(coords, prms1)
        
    # Read the MOPAC output
    # This is counterintuitive, but the actual output file name is derived from
    # that of the input
    read_params = {"nstates":nstates,
                   "filename":F"{wd}/{mopac_input_prefix}{mopac_jobid}.out", 
                   "orbital_space":orbital_space,
                   "multiplicity": multiplicity,
                   "spin_projection": spin_projection}
    info, MO_curr, data_curr = read_mopac_orbital_info(read_params)

    #================= Construct active space ==================
    active_space = None
    if nelec_act_space is None:
        active_space = info["actual_orbital_space"]
    else: 
        occupied = sorted(abs(orb) for orb in info["reference_det"])
        if nelec_act_space < 1 or nelec_act_space > len(occupied):
            raise ValueError(
                f"nelec_act_space must be between 1 and {len(occupied)}"
            )
        min_indx = occupied[-nelec_act_space]
        # Keep the requested number of active reference electrons and only
        # the virtual orbitals that MOPAC included in its CI expansion.  The
        # MO overlap itself may still be evaluated in the full MO space.
        active_space = list(range(min_indx, info["max_vir"] + 1))

    # Get the properties at the previous time-steps
    MO_prev, data_prev = None, None
    if is_first_time:
        # On the first step, assume the current properties are as the previous
        MO_prev = copy.deepcopy(MO_curr)
        data_prev = copy.deepcopy(data_curr)
        coordinates_prev = copy.deepcopy(coords)
    else:
        # Otherwise, retrieve the previously-stored data
        MO_prev = params["MO_prev"].get(itraj, MO_curr).copy()
        data_prev = params["data_prev"].get(itraj, data_curr)
        
    # Do the calculations - time-overlaps, energies, and Hvib
    obj = tmp()
    obj.ham_adi = CMATRIX(nstates, nstates)
    obj.nac_adi = CMATRIX(nstates, nstates)
    obj.hvib_adi = CMATRIX(nstates, nstates)
    obj.basis_transform = CMATRIX(nstates, nstates)
    obj.time_overlap_adi = CMATRIX(nstates, nstates)
    obj.overlap_adi = CMATRIX(nstates, nstates)

    #======================= MO ===============================
    s_mo_orb = MO_curr.T @ MO_curr
    s_mo = np.kron(np.eye(2), s_mo_orb) # Make it doubled - block-matrix

    st_mo_orb = MO_prev.T @ MO_curr
    st_mo = np.kron(np.eye(2), st_mo_orb) # Make it doubled - block-matrix


    #================= Compute CI time-overlaps ============= 
    lowest_orbital = info["actual_orbital_space"][0]
    highest_orbital = info["actual_orbital_space"][-1]
    if info["actual_orbital_space"] != list(range(lowest_orbital, highest_orbital + 1)):
        raise ValueError("orbital_space must be a contiguous range")

    ovlp_params = {"homo_indx":info["homo_indx"],
                   "nocc":info["homo_indx"] - lowest_orbital,
                   "nvirt":highest_orbital - info["homo_indx"],
                   "nelec":info["nelec"], "nstates":nstates,
                   "active_space":active_space,
                   "spin": info["spin"],
                   "spin_projection": info["spin_projection"],
                   }
    if info["multiplicity"] != 1:
        ovlp_params["reference_det"] = info["reference_det"]
    st_ci = ci.overlap(st_mo, data_prev, data_curr, ovlp_params)
    s_ci = ci.overlap(s_mo, data_curr, data_curr, ovlp_params)

    s_ci_inv_curr, s_ci_inv_prev = None, None
    if do_Lowdin==True:
        # Lowding orthogonalization to fight the rounding errors
        s_ci_inv_curr = ortho.lowdin_inverse_sqrt(s_ci)
        
        if is_first_time:
            s_ci_inv_prev = copy.deepcopy(s_ci_inv_curr)
        else:
            s_ci_inv_prev = params["s_ci_inv_prev"].get(itraj, s_ci_inv_curr)

        s_ci = s_ci_inv_curr @ s_ci @ s_ci_inv_curr
        st_ci = s_ci_inv_prev @ st_ci @ s_ci_inv_curr


    #=============== Now, populate the allocated matrices ======================
    for istate in range(nstates):
        energ = 0.0
        if istate > 0:
            #energ = float(0.5 * (data_prev[0][istate-1] + data_curr[0][istate-1]))
            energ = float(data_curr[0][istate-1])

        obj.ham_adi.set(istate, istate, energ * (1.0 + 0.0j))
        obj.hvib_adi.set(istate, istate, energ * (1.0 + 0.0j))
        obj.basis_transform.set(istate, istate, 1.0 + 0.0j)  # assume identity

        for jstate in range(nstates):
            obj.time_overlap_adi.set(
                istate, jstate, float(np.real(st_ci[istate, jstate])) * (1.0 + 0.0j)
            )
            obj.overlap_adi.set(
                istate, jstate, float(np.real(s_ci[istate, jstate])) * (1.0 + 0.0j)
            )

    # Update the Hvib:
    for istate in range(nstates):
        for jstate in range(istate + 1, nstates):
            dij = (obj.time_overlap_adi.get(istate, jstate) - obj.time_overlap_adi.get(jstate, istate)) / (2.0 * dt)
            obj.hvib_adi.set(istate, jstate, dij * (0.0 - 1.0j))
            obj.hvib_adi.set(jstate, istate, dij * (0.0 + 1.0j))

    # Now, make the current the previous and reset the flag `is_first_time` to False
    # Note - we directly modify the input parameters
    
    # ================= Store state =================
    params["MO_prev"][itraj] = MO_curr.copy()
    params["data_prev"][itraj] = copy.deepcopy(data_curr)
    params["s_ci_inv_prev"][itraj] = copy.deepcopy(s_ci_inv_curr)
    params["is_first_time"][itraj] = False

    return obj
