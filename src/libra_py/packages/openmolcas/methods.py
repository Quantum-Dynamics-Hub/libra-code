# This module implements functions for dealing with the outputs from openmolcas package

import os, sys, math, copy, re, subprocess
import numpy as np
import h5py
import warnings
import json
from liblibra_core import *
import util.libutil as comn

from types import SimpleNamespace

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



"""
How the Data Flows — Step by Step

1. Libra reads the trajectory
cp2k.read_trajectory_xyz_file() extracts nuclear coordinates q from each frame of your .xyz file.
This is where Libra's CP2K module acts as a geometry reader.

2. Libra writes + runs OpenMolcas
make_molcas_input() generates a complete OpenMolcas input file using your molcas_run_params.
run_molcas() spawns pymolcas as a subprocess to run the SA-CASSCF calculation.
OpenMolcas produces:
job.out -> CASSCF energies and CI vectors
job.RasOrb -> MO coefficients

3. Libra parses the results
read_molcas_orbital_info() reads both output files and returns:
info -> orbital space metadata (inactive, active, virtual boundaries)
MO_curr -> MO coefficient matrix (shape nbas × nmo) where nbas is Total number of AO basis functions
data_curr -> CI vectors and energies for each state

4. Libra builds determinant cache
build_det_cache() selects the most important Slater determinants (above ci_coeff_thresh), reducing computational cost.
This cache is reused when computing overlaps between consecutive timesteps.

5. Libra computes overlaps
ci_overlap_general() computes the overlap matrix between CI wavefunctions using:
Slater determinant overlaps (derived from MO overlaps)
CI coefficient products
Two overlap matrices are produced:
time_overlap_adi -> ⟨Ψᵢ(t) | Ψⱼ(t+Δt)⟩ (time-overlap for nonadiabatic dynamics)
overlap_adi -> ⟨Ψᵢ(t) | Ψⱼ(t)⟩ (instantaneous overlap)

6. Results saved to disk
ham_adi — adiabatic Hamiltonian (diagonal = CASSCF energies)
hvib_adi — nonadiabatic coupling vector approximation
st_adi — time-overlap between consecutive timesteps
s_adi — state overlap at the same timestep
"""

def _write_gateway_coord(f, labels, coords, nat):
    """
    Internal helper: write the &GATEWAY Coord block into an open file handle.
    Shared between str-mode and dict-mode branches to avoid code duplication.
    """
    f.write("&GATEWAY\n")
    f.write("Coord\n")
    f.write(f"{nat}\n")
    f.write("Bohr\n")                    
    for i in range(nat):
        x = coords.get(3 * i + 0, 0)
        y = coords.get(3 * i + 1, 0)
        z = coords.get(3 * i + 2, 0)
        f.write(f"{labels[i]:4s}  {x:14.8f}  {y:14.8f}  {z:14.8f}\n")


def _normalize_property_requests(nstates, gradient_states=None, nac_pairs=None):
    """Normalize zero-based Libra state selections for OpenMolcas properties."""
    if gradient_states is None:
        gradients = []
    elif isinstance(gradient_states, str) and gradient_states.lower() == "all":
        gradients = list(range(nstates))
    elif isinstance(gradient_states, (int, np.integer)):
        gradients = [int(gradient_states)]
    else:
        gradients = [int(state) for state in gradient_states]

    if nac_pairs is None:
        pairs = []
    elif isinstance(nac_pairs, str) and nac_pairs.lower() == "all":
        pairs = [(i, j) for i in range(nstates) for j in range(i + 1, nstates)]
    elif (isinstance(nac_pairs, (list, tuple)) and len(nac_pairs) == 2
          and all(isinstance(state, (int, np.integer)) for state in nac_pairs)):
        pairs = [tuple(map(int, nac_pairs))]
    else:
        pairs = [tuple(map(int, pair)) for pair in nac_pairs]

    gradients = list(dict.fromkeys(gradients))
    pairs = list(dict.fromkeys(tuple(sorted(pair)) for pair in pairs))
    if any(state < 0 or state >= nstates for state in gradients):
        raise ValueError(f"gradient_states must be in [0, {nstates})")
    if any(len(pair) != 2 or pair[0] == pair[1] or pair[0] < 0 or pair[1] >= nstates
           for pair in pairs):
        raise ValueError(f"nac_pairs must contain distinct state pairs in [0, {nstates})")
    return gradients, pairs


def _resolve_scf_method(method, spin):
    """Resolve the preliminary SCF method for an OpenMolcas RASSCF run."""
    method = str(method).lower()
    if method == "auto":
        return "uhf" if spin > 1 else "rhf"
    if method not in {"rhf", "uhf", "rohf"}:
        raise ValueError("scf_method must be 'auto', 'rhf', 'uhf', or 'rohf'")
    if method == "rhf" and spin > 1:
        raise ValueError("Open-shell multiplicities require scf_method='uhf' or 'rohf'")
    return method


def _active_electron_count(nactel):
    """Extract the number of active electrons from an OpenMolcas NACTEL value."""
    if isinstance(nactel, (int, np.integer)):
        return int(nactel)
    if isinstance(nactel, str):
        fields = nactel.replace(",", " ").split()
    else:
        fields = list(nactel)
    if not fields:
        raise ValueError("nactel must specify the number of active electrons")
    return int(fields[0])
        
        
def make_molcas_input(molcas_input_filename, molcas_run_params, labels, coords):
    """
    Write an OpenMolcas input file (.in) from parameters and coordinates.

    Parameters
    ----------
    molcas_input_filename : str
        Path to the output .in file.
    molcas_run_params : dict or str
        If dict  — must contain keys: basis, charge, spin, title,
                   nactel, inactive, ras2, ciroot, nac_states.
        If str   — treated as a raw keyword string appended after the
                   Coord block (legacy compatibility mode).
    labels : list of str
        Atom labels, e.g. ["N", "N", "C", "H", ...].
    coords : dict-like with .get(key, default)
        Atomic coordinates in Angstroms, indexed as:
            coords[3*i + 0] = x
            coords[3*i + 1] = y
            coords[3*i + 2] = z
    """
    nat = len(labels)

    with open(molcas_input_filename, "w") as f:

        _write_gateway_coord(f, labels, coords, nat)

        if isinstance(molcas_run_params, str):
            f.write(f"{molcas_run_params}\n\n")
            return

        p = molcas_run_params

        # --- &GATEWAY (finish) ---
        f.write(f"Basis={p.get('basis', '6-31G*')}\n")
        f.write("Group=NoSym\n\n")

        # --- &SEWARD ---
        f.write("&SEWARD\n\n")

        spin = int(p.get("spin", 1))
        nactel = p.get("nactel", "6 0 0")
        scf_method = _resolve_scf_method(p.get("scf_method", "auto"), spin)
        if _active_electron_count(nactel) % 2 != (spin - 1) % 2:
            raise ValueError(
                f"nactel={nactel!r} and spin multiplicity {spin} have "
                "incompatible electron parity"
            )

        # OpenMolcas SCF supports RHF and UHF, but not ROHF. For ROHF-like
        # spin-adapted starting orbitals, proceed directly to RASSCF.
        if scf_method != "rohf":
            f.write("&SCF\n")
            if scf_method == "uhf":
                f.write("UHF\n")
            f.write(f"Charge={p.get('charge', 0)}\n")
            f.write(f"Spin={spin}\n\n")

        # --- &RASSCF ---
        f.write("&RASSCF\n")
        f.write(f"Title={p.get('title', 'Molcas Job')}\n")
        f.write("Symmetry=1\n")
        f.write(f"Spin={spin}\n")
        f.write(f"Nactel={nactel}\n")
        f.write(f"Inactive={p.get('inactive', 15)}\n")
        f.write(f"Ras2={p.get('ras2', 6)}\n")
        f.write(f"Ciroot={p.get('ciroot', '2 2 1')}\n")
        f.write(f"Maxiter={p.get('maxiter', 200)}\n")   # increased from 100

        # FIX: tighten threshold — loose THRE differentiates badly
        f.write(f"Thre={p.get('thre', '1.0e-10 1.0e-6 1.0e-6')}\n")

        # FIX: respect prwf from params, don't hardcode
        # PRSD/VecDet must see the full CSF list. Determinant selection is
        # performed later with ci_threshold after repeated SDs are combined.
        f.write(f"PRWF={p.get('prwf', '0.0')}\n")
        f.write("PRSD\n")

        # only write LSHIFT if explicitly requested and small
        lshift = p.get('lshift', None)
        if lshift is not None:
            if float(lshift) > 0.3:
                warnings.warn(
                    f"[make_molcas_input] LSHIFT={lshift} is large (>0.3). "
                    "This can mask non-convergence. Consider removing it "
                    "once FILEORB restart is working.",
                    UserWarning
                )
            f.write(f"LSHIFT={lshift}\n")

        # Feed previous converged orbitals as starting guess.
        fileorb = p.get('fileorb', None)
        if fileorb is not None:
            if os.path.isfile(fileorb):
                f.write(f"FILEORB={fileorb}\n")
                print(f"[make_molcas_input] Using restart orbitals: {fileorb}")
            else:
                warnings.warn(
                    f"[make_molcas_input] FILEORB path does not exist: {fileorb}\n"
                    "Proceeding without orbital restart — "
                    "convergence may be unreliable.",
                    UserWarning
                )
        elif scf_method == "uhf":
            orbital_set = str(p.get("uhf_orbital_set", "beta")).lower()
            if orbital_set not in {"alpha", "beta"}:
                raise ValueError("uhf_orbital_set must be 'alpha' or 'beta'")
            f.write("FILEORB=$Project.UhfOrb\n")
            f.write(f"AlphaOrBeta={'1' if orbital_set == 'alpha' else '-1'}\n")
        elif scf_method == "rohf":
            f.write("FILEORB=$Project.GssOrb\n")

        f.write("\n")

        nstates = p.get("nstates", _infer_nstates_from_ciroot(p.get("ciroot")))
        if nstates is None:
            nstates = 1
        gradient_states, nac_pairs = _normalize_property_requests(
            nstates, p.get("gradient_states"), p.get("nac_pairs")
        )

        # ``nac_states`` is the legacy one-based single-pair option.
        if not nac_pairs and p.get("nac_states") is not None:
            i, j = p["nac_states"]
            nac_pairs = [(int(i) - 1, int(j) - 1)]

        for state in gradient_states:
            f.write("&ALASKA\n")
            f.write(f"ROOT={state + 1}\n")
            f.write("SHOW\n\n")

        for i, j in nac_pairs:
            f.write("&ALASKA\n")
            f.write(f"NAC={i + 1} {j + 1}\n")
            if p.get("nac_nocsf", False):
                f.write("NOCSF\n")
            f.write("SHOW\n\n")

def _print_output_tail(path, n=80):
    """Print last n lines of a file, then scan for error keywords for debugging"""
    if not os.path.isfile(path):
        print(f"  [Output file missing: {path}]")
        return

    with open(path, "r") as f:
        lines = f.readlines()

    print(f"── Last {n} lines of {path} (total: {len(lines)} lines) ──")
    for line in lines[-n:]:
        print(line, end="")

    
    error_kws = ["convergence", "error", "fatal", "abort", "not found", "stop", "--- stop"]
    hits = [
        (i + 1, line.rstrip())
        for i, line in enumerate(lines)
        if any(kw in line.lower() for kw in error_kws)
    ]

    if hits:
        print(f"\n── Error keyword scan ({len(hits)} hits) ──")
        for lineno, text in hits:
            print(f"  Line {lineno:5d}: {text}")

            
def run_molcas(coords, params_):
    """
    Run an OpenMolcas single-point + NAC calculation.
    Returns (output_path, rasorb_path).
    """
    params = dict(params_)

    labels        = params["atom_labels"]
    exe           = params.get("exe", "pymolcas")
    molcas_wd     = params.get("working_directory", "molcas_wd")
    molcas_jobid  = params.get("molcas_jobid", "job_0000")
    input_prefix  = params.get("input_prefix", "input_")

    mrp = params.get("molcas_run_params")
    if mrp is not None:
        molcas_run_params = dict(mrp)
        molcas_run_params.setdefault("title", f"Molcas_{molcas_jobid}")
        molcas_run_params.setdefault("nac_states", params.get("nac_states"))
        if "fileorb" not in molcas_run_params:
            molcas_run_params["fileorb"] = params.get("fileorb", None)
    else:
        molcas_run_params = {
            "basis":      params["basis"],
            "charge":     params["charge"],
            "spin":       params["spin"],
            "title":      f"Molcas_{molcas_jobid}",
            "nactel":     params["nactel"],
            "inactive":   params["inactive"],
            "ras2":       params["ras2"],
            "ciroot":     params["ciroot"],
            "maxiter":    params.get("maxiter", 200),
            "thre":       params.get("thre", "1.0e-10 1.0e-6 1.0e-6"),
            "prwf":       params.get("prwf", "1.0d-10"),
            "lshift":     params.get("lshift", None),
            "scf_method": params.get("scf_method", "auto"),
            "uhf_orbital_set": params.get("uhf_orbital_set", "beta"),
            "nac_states": params.get("nac_states"),
            "fileorb":    params.get("fileorb", None),
        }

    # critical fix
    if molcas_run_params.get("fileorb") is not None:
        molcas_run_params["fileorb"] = os.path.abspath(molcas_run_params["fileorb"])

    print("[DEBUG run_molcas] final fileorb =", molcas_run_params.get("fileorb"))
    print("[DEBUG run_molcas] final thre    =", molcas_run_params.get("thre"))
    print("[DEBUG run_molcas] final lshift  =", molcas_run_params.get("lshift"))

    project_name = f"{input_prefix}{molcas_jobid}"

    molcas_input_filename  = f"{project_name}.in"
    molcas_output_filename = f"{project_name}.out"

    molcas_input_path  = os.path.join(molcas_wd, molcas_input_filename)
    molcas_output_path = os.path.join(molcas_wd, molcas_output_filename)

    rasorb_flat = os.path.join(molcas_wd, f"{project_name}.RasOrb")
    rasorb_sub  = os.path.join(molcas_wd, project_name, f"{project_name}.RasOrb")

    os.makedirs(molcas_wd, exist_ok=True)

    molcas_env = os.environ.copy()
    molcas_env["MOLCAS_PROJECT"]    = project_name
    molcas_env["MOLCAS_SCRATCH"]    = os.path.abspath(molcas_wd)
    molcas_env["MOLCAS_SUBMIT_DIR"] = os.path.abspath(molcas_wd)

    make_molcas_input(molcas_input_path, molcas_run_params, labels, coords)

    print(f"[run_molcas] Running    : {exe} {molcas_input_filename}")
    print(f"[run_molcas] Working dir: {molcas_wd}")
    print(f"[run_molcas] Output file: {molcas_output_path}")
    print(f"[run_molcas] Restart orb: {molcas_run_params.get('fileorb')}")

    with open(molcas_output_path, "w") as fout:
        result = subprocess.run(
            [exe, molcas_input_filename],
            cwd=molcas_wd,
            env=molcas_env,
            check=False,
            stdout=fout,
            stderr=subprocess.STDOUT,
        )

    if result.returncode != 0:
        print(f"\n[run_molcas] OpenMolcas FAILED (exit code {result.returncode})")
        print(f"[run_molcas] Input : {molcas_input_path}")
        print(f"[run_molcas] Output: {molcas_output_path}\n")
        _print_output_tail(molcas_output_path, n=80)
        raise RuntimeError(
            f"OpenMolcas failed (exit code {result.returncode}).\n"
            f"  Input : {molcas_input_path}\n"
            f"  Output: {molcas_output_path}"
        )

    if os.path.isfile(rasorb_flat):
        rasorb_path = rasorb_flat
    elif os.path.isfile(rasorb_sub):
        rasorb_path = rasorb_sub
    else:
        raise FileNotFoundError(
            f"OpenMolcas completed but RasOrb not found.\n"
            f"  Tried: {rasorb_flat}\n"
            f"  Tried: {rasorb_sub}\n"
        )

    return molcas_output_path, rasorb_path

def read_rasscf_energies(h5_path):
    """
    Read RASSCF electronic state energies from an OpenMolcas HDF5 file.

    This function extracts the energies of the RASSCF roots from an
    OpenMolcas-generated HDF5 file. It first searches for the standard
    dataset locations and, if they are not found, performs a recursive
    search for datasets whose names contain the keywords "energy" or "root".

    Parameters
    ----------
    h5_path : str
        Path to the OpenMolcas HDF5 (.h5) file.

    Returns
    -------
    list of float
        A list containing the energies of all available RASSCF roots.
        Returns an empty list if no energy datasets are found.

    Notes
    -----
    The function checks the following datasets in order:
        - /RASSCF/ROOT_ENERGIES
        - /RASSCF/CM_ENERGY

    If none of these datasets exist, it searches the entire HDF5 file
    for datasets with names containing "energy" or "root".

    Example
    -------
    >>> energies = read_rasscf_energies("molcas.rasscf.h5")
    >>> print(energies)
    [-382.745231, -382.612184, -382.503917]
    """
    energies = []
    with h5py.File(h5_path, "r") as f:
        candidate_paths = [
            "/RASSCF/ROOT_ENERGIES",
            "/RASSCF/CM_ENERGY",
            "/RASSCF/ROOT_ENERGIES",
        ]
        for path in candidate_paths:
            if path in f:
                data = f[path][()]
                energies = data.flatten().tolist()
                break
        if not energies:
            def _search(name, obj):
                if isinstance(obj, h5py.Dataset):
                    name_lower = name.lower()
                    if "energy" in name_lower or "root" in name_lower:
                        data = obj[()].flatten()
                        energies.extend(data.tolist())
            f.visititems(_search)
    return energies


def read_alaska_vectors(out_file):
    """Read Cartesian gradients and derivative couplings printed by ALASKA.

    The returned list contains ``(kind, values)`` pairs in output order, where
    ``kind`` is ``"gradient"`` or ``"nac"`` and ``values`` has shape
    ``(natoms, 3)``. OpenMolcas prints gradients in Hartree/Bohr and total
    derivative couplings in Bohr^-1.
    """
    if not os.path.isfile(out_file):
        raise FileNotFoundError(f"OpenMolcas output file not found: {out_file}")

    marker = re.compile(r"(Molecular gradients|Total derivative coupling)", re.I)
    number = re.compile(r"^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?$")
    blocks = []
    current_kind = None
    rows = []

    def finish_block():
        nonlocal current_kind, rows
        if current_kind is not None and rows:
            blocks.append((current_kind, np.asarray(rows, dtype=np.float64)))
        current_kind = None
        rows = []

    with open(out_file, "r") as output:
        for line in output:
            match = marker.search(line)
            if match:
                finish_block()
                current_kind = ("gradient" if "molecular" in match.group(1).lower()
                                else "nac")
                continue
            if current_kind is None:
                continue

            numeric = [
                float(token.replace("D", "E").replace("d", "e"))
                for token in line.split() if number.fullmatch(token)
            ]
            if len(numeric) >= 3:
                rows.append(numeric[-3:])
            elif rows and set(line.strip()) == {"-"}:
                finish_block()
    finish_block()
    return blocks


def read_alaska_properties(out_file, natoms, gradient_states=None, nac_pairs=None,
                           nstates=None):
    """Associate ALASKA tables with requested zero-based states and pairs."""
    if nstates is None:
        requested = list(gradient_states or [])
        requested += [state for pair in (nac_pairs or []) for state in pair]
        nstates = max(requested, default=-1) + 1
    gradient_states, nac_pairs = _normalize_property_requests(
        nstates, gradient_states, nac_pairs
    )
    blocks = read_alaska_vectors(out_file)
    gradient_blocks = [values for kind, values in blocks if kind == "gradient"]
    nac_blocks = [values for kind, values in blocks if kind == "nac"]

    if len(gradient_blocks) != len(gradient_states):
        raise ValueError(
            f"Found {len(gradient_blocks)} ALASKA gradient tables, expected "
            f"{len(gradient_states)}"
        )
    if len(nac_blocks) != len(nac_pairs):
        raise ValueError(
            f"Found {len(nac_blocks)} ALASKA NAC tables, expected {len(nac_pairs)}"
        )
    for values in gradient_blocks + nac_blocks:
        if values.shape != (natoms, 3):
            raise ValueError(
                f"ALASKA Cartesian table has shape {values.shape}, expected ({natoms}, 3)"
            )

    return dict(zip(gradient_states, gradient_blocks)), dict(zip(nac_pairs, nac_blocks))


def check_energy_continuity(energies_curr, energies_prev, itraj, timestep,
                             thresh_eV=0.3, hartree_to_eV=27.2114):
    """
    Flag any state where the energy jumps by more than thresh_eV between steps.
    Every flagged step is a bad-convergence geometry.
    Returns True if all states pass, False if any flag is raised.
    """
    if energies_prev is None:
        return True   # nothing to compare on the first step

    all_ok = True
    for i, (e_curr, e_prev) in enumerate(zip(energies_curr, energies_prev)):
        delta_eV = abs(e_curr - e_prev) * hartree_to_eV
        if delta_eV > thresh_eV:
            print(
                f"[ENERGY SANITY]  traj={itraj} step={timestep} "
                f"State {i+1}: |ΔE| = {delta_eV:.4f} eV  "
                f"({e_prev:.8f} → {e_curr:.8f} Ha) "
                f"— likely bad-convergence geometry. "
                f"Check RasOrb restart and active-space continuity."
            )
            all_ok = False

    if all_ok:
        print(f"[ENERGY SANITY]  traj={itraj} step={timestep} "
              f"All states continuous (thresh={thresh_eV} eV)")
    return all_ok


def _combine_determinants(determinants, coefficients, ci_threshold=0.0):
    """Combine repeated determinants and discard weights below ``ci_threshold``."""
    combined = {}
    order = []
    for determinant, coefficient in zip(determinants, coefficients):
        if determinant not in combined:
            combined[determinant] = 0.0
            order.append(determinant)
        combined[determinant] += coefficient

    kept_determinants = []
    kept_coefficients = []
    for determinant in order:
        coefficient = combined[determinant]
        if abs(coefficient) ** 2 >= ci_threshold:
            kept_determinants.append(determinant)
            kept_coefficients.append(coefficient)
    return kept_determinants, kept_coefficients


def read_vecdet(filepath, ci_threshold=0.0):
    """Read an OpenMolcas PRSD ``VecDet`` determinant expansion.

    The returned occupation tuples use ``2`` for a doubly occupied spatial
    orbital, ``1`` for alpha, ``-1`` for beta, and ``0`` for empty. For example,
    the file row ``0.70710678 ab0`` becomes ``((1, -1, 0), 0.70710678)``.
    Repeated determinants are combined, then determinants with
    ``|C|**2 < ci_threshold`` are discarded.

    A VecDet file describes the M_S component printed by OpenMolcas, normally
    the highest-weight component M_S=S. Use :func:`expand_ms_projections` to
    construct the other components of a spin multiplet.
    """
    if not os.path.isfile(filepath):
        raise FileNotFoundError(f"OpenMolcas VecDet file not found: {filepath}")

    occupation_map = {"2": 2, "0": 0, "a": 1, "b": -1}
    determinants = []
    coefficients = []
    with open(filepath, "r") as vecdet:
        next(vecdet, None)  # inactive orbitals per symmetry
        for line_number, line in enumerate(vecdet, start=2):
            fields = line.split()
            if len(fields) != 2:
                if fields:
                    raise ValueError(f"Malformed VecDet row {line_number} in {filepath}")
                continue
            coefficient = float(fields[0].replace("D", "E").replace("d", "e"))
            occupation = fields[1].lower()
            if not occupation or any(char not in occupation_map for char in occupation):
                raise ValueError(
                    f"Invalid determinant occupation {occupation!r} at "
                    f"{filepath}:{line_number}"
                )
            determinants.append(tuple(occupation_map[char] for char in occupation))
            coefficients.append(coefficient)

    return _combine_determinants(determinants, coefficients, ci_threshold)


def _read_hdf5_ci_vectors(h5_path):
    """Return the full-precision CSF coefficient array when available."""
    if not h5_path or not os.path.isfile(h5_path):
        return None
    with h5py.File(h5_path, "r") as h5_file:
        if "CI_VECTORS" not in h5_file:
            return None
        return np.asarray(h5_file["CI_VECTORS"], dtype=np.float64)


def _spin_orbital_occupation(determinant):
    """Return occupied indices in OpenMolcas's interleaved a,b convention."""
    occupied = []
    for orbital, occupation in enumerate(determinant):
        if occupation in (1, 2):
            occupied.append(2 * orbital)
        if occupation in (-1, 2):
            occupied.append(2 * orbital + 1)
    return occupied


def _lower_spin_state(determinants, coefficients, spin, ms):
    """Apply S_- to a normalized |S,M_S> determinant expansion."""
    lowered_determinants = []
    lowered_coefficients = []
    normalization = math.sqrt(spin * (spin + 1.0) - ms * (ms - 1.0))

    for determinant, coefficient in zip(determinants, coefficients):
        occupied = _spin_orbital_occupation(determinant)
        for orbital, occupation in enumerate(determinant):
            if occupation != 1:
                continue

            alpha = 2 * orbital
            beta = alpha + 1
            annihilation_position = occupied.index(alpha)
            intermediate = [index for index in occupied if index != alpha]
            creation_position = sum(index < beta for index in intermediate)
            sign = -1.0 if (annihilation_position + creation_position) % 2 else 1.0

            lowered = list(determinant)
            lowered[orbital] = -1
            lowered_determinants.append(tuple(lowered))
            lowered_coefficients.append(sign * coefficient / normalization)

    return _combine_determinants(lowered_determinants, lowered_coefficients, 0.0)


def expand_ms_projections(determinants, coefficients, spin_multiplicity):
    """Expand each highest-weight spin state into all distinct M_S projections.

    Parameters
    ----------
    determinants, coefficients : list[list]
        Outer lists are spin-free CASSCF roots. Each input root must represent
        a pure-spin, normalized highest-weight state ``|S, M_S=S>``. Inner
        determinant tuples use ``2, 1, -1, 0`` for double, alpha, beta, and
        empty spatial-orbital occupations, respectively.
    spin_multiplicity : int
        ``2*S + 1``. Any positive multiplicity is supported: 1 is a singlet,
        2 a doublet, 3 a triplet, 4 a quartet, and so forth.

    Returns
    -------
    expanded_determinants, expanded_coefficients : list[list]
        One determinant expansion per explicit M_S state. Ordering is
        root-major; within each root it is ``M_S=S, S-1, ..., -S``.
    labels : list[tuple[int, int]]
        ``(root_index, 2*M_S)`` for every returned state. Twice M_S is stored
        as an integer so half-integer projections have an exact label.

    Notes
    -----
    Lower projections are generated recursively with the normalized spin
    lowering operator ``S_-``. Thus a triplet root has labels
    ``[(root, 2), (root, 0), (root, -2)]`` and a doublet has
    ``[(root, 1), (root, -1)]``. With two triplet roots, the complete ordering
    is ``[(0, 2), (0, 0), (0, -2), (1, 2), (1, 0), (1, -2)]``.

    Examples
    --------
    A two-orbital triplet highest-weight determinant expands as follows::

        dets = [[(1, 1)]]
        coeffs = [[1.0]]
        dets_ms, coeffs_ms, labels = expand_ms_projections(dets, coeffs, 3)
        # labels == [(0, 2), (0, 0), (0, -2)]
        # dets_ms[0] == [(1, 1)]                 # M_S = +1
        # dets_ms[1] == [(-1, 1), (1, -1)]      # M_S =  0
        # dets_ms[2] == [(-1, -1)]               # M_S = -1

    The signs and determinant order in an expansion follow the interleaved
    spin-orbital convention ``(1a, 1b, 2a, 2b, ...)`` used by OpenMolcas.
    """
    multiplicity = int(spin_multiplicity)
    if multiplicity < 1:
        raise ValueError("spin_multiplicity must be positive")
    spin = 0.5 * (multiplicity - 1)

    expanded_determinants = []
    expanded_coefficients = []
    labels = []
    for root, (root_determinants, root_coefficients) in enumerate(
            zip(determinants, coefficients)):
        projection_determinants = list(root_determinants)
        projection_coefficients = list(root_coefficients)
        ms = spin
        for projection in range(multiplicity):
            expanded_determinants.append(projection_determinants)
            expanded_coefficients.append(projection_coefficients)
            labels.append((root, multiplicity - 1 - 2 * projection))
            if projection + 1 < multiplicity:
                projection_determinants, projection_coefficients = _lower_spin_state(
                    projection_determinants, projection_coefficients, spin, ms
                )
                ms -= 1.0
    return expanded_determinants, expanded_coefficients, labels


def read_ci_vectors(out_file, expected_states=None, h5_path=None, ci_threshold=0.0,
                    vecdet_prefix=None):
    """Read OpenMolcas wave functions expanded in Slater determinants.

    ``VecDet.<root>`` files are preferred because PRSD writes determinant
    coefficients directly. Otherwise, CSF coefficients are read from HDF5 when
    possible and multiplied by the CSF-to-determinant expansions printed in the
    output. Repeated determinants are combined before ``ci_threshold`` is
    applied to their final squared coefficient magnitudes.

    Returns two root-indexed nested lists. For example, two roots may produce
    ``confs = [[(2, 0)], [(1, -1), (-1, 1)]]`` and
    ``CI = [[1.0], [0.7071, -0.7071]]``. These are spin-free roots, not an
    enumeration of all M_S components. OpenMolcas normally prints the
    highest-weight M_S=S component for a requested multiplicity.
    """
    if ci_threshold < 0.0:
        raise ValueError("ci_threshold must be non-negative")

    if vecdet_prefix is None:
        vecdet_prefix = os.path.splitext(out_file)[0] + ".VecDet"

    n_vecdet = 0
    while os.path.isfile(f"{vecdet_prefix}.{n_vecdet + 1}"):
        n_vecdet += 1
    if n_vecdet:
        nstates = max(n_vecdet, expected_states or 0)
        all_determinants = []
        all_coefficients = []
        for root in range(1, nstates + 1):
            path = f"{vecdet_prefix}.{root}"
            if os.path.isfile(path):
                determinants, coefficients = read_vecdet(path, ci_threshold)
            else:
                determinants, coefficients = [], []
            all_determinants.append(determinants)
            all_coefficients.append(coefficients)
        return all_determinants, all_coefficients

    csf_coefficients = _read_hdf5_ci_vectors(h5_path)
    determinant_map = {"2": 2, "0": 0, "a": 1, "b": -1}

    root_pat = re.compile(
        r"printout\s+of\s+CI-coefficients.*for\s+root\s+(\d+)",
        re.IGNORECASE
    )

    expansion_pat = re.compile(
        r"^([+-])\s+sqrt\(\s*(\d+)\s*/\s*(\d+)\s*\)\s+\|([20ab]+)\|$",
        re.IGNORECASE,
    )
    root_confs = {}
    root_CIs = {}

    current_root   = None
    current_confs  = []
    current_CIs    = []
    current_csf_coefficient = None

    CI_END = [
        "Natural orbitals",
        "RASSCF results",
        "CASSCF results",
        "--- Stop Module",
        "++  Convergence",
    ]

    with open(out_file, "r") as f:
        for raw in f:
            stripped = raw.strip()
            m = root_pat.search(stripped)
            if m:
                if current_root is not None:
                    root_confs[current_root], root_CIs[current_root] = _combine_determinants(
                        current_confs, current_CIs, ci_threshold
                    )
                current_root  = int(m.group(1)) - 1   
                current_confs = []
                current_CIs   = []
                current_csf_coefficient = None
                continue

            if current_root is None:
                continue

          
            if any(p in stripped for p in CI_END):
                root_confs[current_root], root_CIs[current_root] = _combine_determinants(
                    current_confs, current_CIs, ci_threshold
                )
                current_root  = None
                current_confs = []
                current_CIs   = []
                current_csf_coefficient = None
                continue

            if not stripped:
                continue
            expansion = expansion_pat.match(stripped)
            if expansion and current_csf_coefficient is not None:
                sign = 1.0 if expansion.group(1) == "+" else -1.0
                factor = sign * math.sqrt(int(expansion.group(2)) / int(expansion.group(3)))
                occupation = expansion.group(4).lower()
                determinant = tuple(determinant_map[char] for char in occupation)
                current_confs.append(determinant)
                current_CIs.append(current_csf_coefficient * factor)
                continue
            if any(h in stripped for h in ["energy=", "conf/sym", "Coeff", "Weight"]):
                continue

            parts = stripped.split()
            if len(parts) < 3:
                continue

            try:
                int(parts[0])   
            except ValueError:
                continue

            try:
                csf_index = int(parts[0]) - 1
                printed_coefficient = float(
                    parts[-2].replace("D", "E").replace("d", "e")
                )
                if (csf_coefficients is not None
                        and current_root < csf_coefficients.shape[0]
                        and csf_index < csf_coefficients.shape[1]):
                    current_csf_coefficient = csf_coefficients[current_root, csf_index]
                else:
                    current_csf_coefficient = printed_coefficient

            except (ValueError, IndexError, KeyError):
                continue

    if current_root is not None:
        root_confs[current_root], root_CIs[current_root] = _combine_determinants(
            current_confs, current_CIs, ci_threshold
        )

    if not root_confs:
        n_found = 0
    else:
        n_found = max(root_confs.keys()) + 1

    if expected_states is not None:
        n_found = max(n_found, expected_states)

    all_confs = [root_confs.get(i, []) for i in range(n_found)]
    all_CIs   = [root_CIs.get(i, [])   for i in range(n_found)]

    n_total = sum(len(c) for c in all_confs)
    if root_confs and n_total == 0:
        raise ValueError(
            "No CSF-to-determinant expansions were found in the OpenMolcas "
            "output. Add PRSD to &RASSCF or provide the generated VecDet files."
        )
    print(f"[DEBUG] Found {len(all_confs)} CI vectors with {n_total} total determinants")

    for i, state_confs in enumerate(all_confs):
        if state_confs:
            print(f"[DEBUG] State {i+1}: first determinant = {state_confs[0]}, "
                  f"norb = {len(state_confs[0])}")
        else:
            print(f"[WARNING]  State {i+1}: 0 configurations — "
                  f"check PRWF threshold and Ciroot in Molcas input")

    return all_confs, all_CIs


def read_rasorb(filepath):
    """
    Reads Molecular Orbital (MO) coefficients from an OpenMolcas .RasOrb / .ScfOrb file.
    This parses the standard INPORB format. Assumes C1 symmetry (Group=NoSym).
    
    Returns
    -------
    mos : np.ndarray
        MO coefficient matrix of shape (nbas, nmo).
    """
    if not os.path.isfile(filepath):
        raise FileNotFoundError(f"RasOrb file not found: {filepath}")

    with open(filepath, 'r') as f:
        lines = f.readlines()

    orb_start = -1
    for i, line in enumerate(lines):
        if line.startswith("#ORB"):
            orb_start = i
            break

    if orb_start == -1:
        raise ValueError(f"Could not find '#ORB' section in {filepath}")

    coeffs = []
    current_orb = []

    for line in lines[orb_start+1:]:
        line = line.strip()
        
        if line.startswith("#"):
            break

        if line.startswith("* ORBITAL"):
           
            if current_orb:
                coeffs.append(current_orb)
                current_orb = []
        elif line.startswith("*"):
          
            continue
        elif line:
           
            parts = line.upper().replace('D', 'E').split()
            current_orb.extend([float(x) for x in parts])

    # Append the last orbital
    if current_orb:
        coeffs.append(current_orb)

    # Convert to numpy array. 
    # 'coeffs' is shape (nmo, nbas). Transpose to get (nbas, nmo)
    mos = np.array(coeffs).T
    
    return mos

def read_molcas_orbital_info(params):
    """
    Parse an OpenMolcas output file and RasOrb file to extract
    orbital space information, MO coefficients, energies, and CI vectors.

    Parameters
    ----------
    params : dict
        Must contain 'filename' or 'output_file', and optionally 'rasorb_file'.

    Returns
    -------
    info : dict
        Orbital space metadata.
    mos : np.ndarray
        MO coefficient matrix, shape (nbas, nmo).
    calc_data : dict
        Contains keys 'energies', 'confs', 'CI'.
    """
    params = dict(params)
    out_file = params.get("filename", params.get("output_file", "job.out"))
    rasorb_file = params.get("rasorb_file", "job.RasOrb")
    h5_path = params.get("h5_path", rasorb_file.replace(".RasOrb", ".rasscf.h5"))
    ci_threshold = float(params.get("ci_threshold", 0.0))

    if not os.path.isfile(out_file):
        raise FileNotFoundError(
            f"Cannot find the OpenMolcas output file: '{out_file}'\n"
            f"Hint: Current working directory is: {os.getcwd()}"
        )

    if not os.path.isfile(rasorb_file):
        raise FileNotFoundError(
            f"Cannot find the RasOrb file: '{rasorb_file}'\n"
            f"Hint: Current working directory is: {os.getcwd()}"
        )

    n_inactive = 0
    n_active_orb = 0
    n_act_elec = 0
    nbas = 0

    with open(out_file, "r") as f:
        for line in f:
            if "Number of inactive orbitals" in line:
                n_inactive = int(line.split()[-1])
            elif "Number of active orbitals" in line:
                n_active_orb = int(line.split()[-1])
            elif "Number of electrons in active shells" in line:
                n_act_elec = int(line.split()[-1])
            elif "Number of basis functions" in line:
                nbas = int(line.split()[-1])

            if all(x > 0 for x in [n_inactive, n_active_orb, n_act_elec, nbas]):
                break

    if nbas == 0:
        raise RuntimeError(
            f"Could not determine 'Number of basis functions' from '{out_file}'. "
            "Is this a valid OpenMolcas output file?"
        )

    nelec = 2 * n_inactive + n_act_elec

    min_occ = 1
    max_occ = n_inactive
    min_act = n_inactive + 1
    max_act = n_inactive + n_active_orb
    min_vir = max_act + 1
    max_vir = nbas

    mos = read_rasorb(rasorb_file)
    energies = read_rasscf_energies(h5_path)
    confs, ci = read_ci_vectors(
        out_file,
        expected_states=params.get("nstates"),
        h5_path=h5_path,
        ci_threshold=ci_threshold,
    )

    info = {
        "nao": nbas,
        "nmo": nbas,
        "nelec": nelec,
        "nocc": n_inactive,
        "nact": n_active_orb,
        "nact_elec": n_act_elec,

        "min_occ": min_occ,
        "max_occ": max_occ,
        "min_active": min_act,
        "max_active": max_act,
        "min_vir": min_vir,
        "max_vir": max_vir,

        "actual_orbital_space": list(range(min_act, max_act + 1)),

        "boundaries": {
            "inactive_range": (min_occ, max_occ),
            "active_range": (min_act, max_act),
            "virtual_range": (min_vir, max_vir),
        },
    }

    calc_data = {
        "energies": energies,
        "confs": confs,
        "CI": ci,
    }

    print(json.dumps(info, indent=4))

    return info, mos, calc_data


def read_ao_overlap(job_file, nbas):
    """
    This may not be the perfect method to read the overlap matrix for now it is fine.
    Read AO overlap matrix from OpenMolcas.

    Tries in order:
      1. scf.h5  (primary — AO_OVERLAP_MATRIX lives here)
      2. rasscf.h5 (fallback HDF5)
      3. Text output (if 'Overlap' keyword was used in &SEWARD)
      4. Identity matrix fallback (with warning)
    """
    h5_candidates = [
        job_file.replace(".out", ".scf.h5"),    
        job_file.replace(".out", ".rasscf.h5"), 
        job_file.replace(".out", ".h5"),
    ]

    for h5_path in h5_candidates:
        if not os.path.isfile(h5_path):
            continue
        try:
            import h5py
            with h5py.File(h5_path, "r") as fh:
                print(f"[DEBUG] HDF5 keys in {h5_path}: {list(fh.keys())}")

                for key in ["AO_OVERLAP_MATRIX", "OVERLAP", "overlap", "S"]:
                    if key not in fh:
                        continue

                    raw = np.array(fh[key][:], dtype=np.float64).ravel()

                    expected_packed = nbas * (nbas + 1) // 2

                    if raw.size == expected_packed:
                        
                        S = np.zeros((nbas, nbas), dtype=np.float64)
                        idx = 0
                        for row in range(nbas):
                            for col in range(row + 1):
                                S[row, col] = raw[idx]
                                S[col, row] = raw[idx]
                                idx += 1

                    elif raw.size == nbas * nbas:
                       
                        S = raw.reshape(nbas, nbas)

                    else:
                        print(
                            f"[DEBUG] Unexpected array size {raw.size} "
                            f"(expected {expected_packed} packed or "
                            f"{nbas*nbas} full) — skipping key '{key}'"
                        )
                        continue

                    diag = np.diag(S)
                    if not np.allclose(diag, 1.0, atol=1e-3):
                        print(
                            f"[DEBUG] Key '{key}' diagonal not ~1.0 "
                            f"(max={diag.max():.4f}, min={diag.min():.4f}) — skipping"
                        )
                        continue

                    print(
                        f"[DEBUG] Loaded S_ao from HDF5 key '{key}' "
                        f"in {os.path.basename(h5_path)}, shape {S.shape}"
                    )
                    return S

        except ImportError:
            print("[DEBUG] h5py not available — skipping HDF5 read")
        except Exception as e:
            print(f"[DEBUG] HDF5 read failed for {h5_path}: {e}")

    # ── Option 2: Parse text output 
    S = _parse_overlap_from_output(job_file, nbas)
    if S is not None:
        return S

    # ── Option 3: Identity fallback 
    warnings.warn(
        "[read_ao_overlap] Could not find AO overlap matrix. "
        "Falling back to identity — MO overlaps will be WRONG for non-orthonormal AO basis.\n"
        "Fix: add 'Overlap' keyword to &SEWARD in your OpenMolcas input.",
        UserWarning,
        stacklevel=2,
    )
    return np.eye(nbas, dtype=np.float64)

def _parse_overlap_from_output(out_file, nbas):
    """
    Parse AO overlap matrix from OpenMolcas text output.
    Requires 'Overlap' keyword in &SEWARD.
    """
    S = np.zeros((nbas, nbas), dtype=np.float64)
    filled = np.zeros((nbas, nbas), dtype=bool)

    in_ovlp   = False
    current_cols = []

    col_header_re = re.compile(r"^\s*(\d+)(\s+\d+)*\s*$")

    with open(out_file, "r") as f:
        for line in f:
            stripped = line.strip()

            if "OVERLAP MATRIX" in stripped.upper():
                in_ovlp      = True
                current_cols = []
                S[:]         = 0.0
                filled[:]    = False
                continue

            if not in_ovlp:
                continue

            if not stripped or stripped.startswith("*"):
                continue

            if stripped.startswith("---"):
                break
            
            if stripped[0].isalpha():
                break

            
            tokens = stripped.split()
            if all(t.lstrip("-").isdigit() for t in tokens):
                try:
                    current_cols = [int(t) - 1 for t in tokens]  # 1-based → 0-based
                    continue
                except ValueError:
                    pass

            if not current_cols:
                continue

            parts = stripped.split()
            if len(parts) < 2:
                continue

            try:
                row_idx = int(parts[0]) - 1  # 1-based → 0-based
            except ValueError:
                continue

            if not (0 <= row_idx < nbas):
                continue

            try:
                vals = [
                    float(x.replace("D", "E").replace("d", "e"))
                    for x in parts[1:]
                ]
            except ValueError:
                continue

            n_vals = min(len(vals), len(current_cols))
            for i in range(n_vals):
                col_idx = current_cols[i]
                if 0 <= col_idx < nbas:
                    S[row_idx, col_idx] = vals[i]
                    filled[row_idx, col_idx] = True

    if not filled.any():
        print("[DEBUG] _parse_overlap_from_output: no data found — "
              "add 'Overlap' keyword to &SEWARD")
        return None

    diag = np.diag(S)
    if not np.allclose(diag, 1.0, atol=1e-3):
        print(f"[DEBUG] Diagonal not ~1.0: min={diag.min():.4f}, max={diag.max():.4f}")
        return None

    S = 0.5 * (S + S.T)

    print(f"[DEBUG] Parsed S_ao from text output, shape {S.shape}")
    return S

def occ_tuple_to_alpha_beta(occ_tuple, active_space):
    """
    Convert a CASSCF occupation tuple to lists of occupied
    alpha and beta 1-based MO indices (active orbitals only).

    Occupation convention from OpenMolcas:
        2  → doubly occupied  (alpha + beta)
        1  → singly occupied, alpha
       -1  → singly occupied, beta
        0  → unoccupied

    The determinant phase convention is based on interleaved spin orbitals
    ``(1a, 1b, 2a, 2b, ...)``. This function returns separate alpha and beta
    lists for determinant factorization; :func:`slater_det_overlap` restores
    the corresponding interleaved-order permutation phase.
    """
    alpha_orbs = []
    beta_orbs  = []

    for i, occ in enumerate(occ_tuple):
        orb_idx = active_space[i]   
        if occ == 2:
            alpha_orbs.append(orb_idx)
            beta_orbs.append(orb_idx)
        elif occ == 1:
            alpha_orbs.append(orb_idx)
        elif occ == -1:
            beta_orbs.append(orb_idx)
        # occ == 0 → unoccupied, skip

    return alpha_orbs, beta_orbs

def build_full_orbital_lists(occ_tuple, active_space, inactive_orbs):
    """
    Build full alpha and beta occupied MO index lists including
    inactive (always doubly occupied) orbitals.

    Parameters
    ----------
    occ_tuple     : occupation tuple for active orbitals only
    active_space  : list of 1-based MO indices for active orbitals
    inactive_orbs : list of 1-based MO indices for inactive orbitals

    Returns
    -------
    alpha_orbs, beta_orbs : sorted lists of 1-based MO indices
    """
    alpha_orbs = list(inactive_orbs)
    beta_orbs  = list(inactive_orbs)

    act_alpha, act_beta = occ_tuple_to_alpha_beta(occ_tuple, active_space)

    alpha_orbs += act_alpha
    beta_orbs  += act_beta

    return sorted(alpha_orbs), sorted(beta_orbs)

def slater_det_overlap(alpha_K, beta_K, alpha_L, beta_L, S_mo):
    """
    Compute the overlap between two Slater determinants K and L:

        <Phi_K | Phi_L> = det(S_alpha_KL) * det(S_beta_KL)

    Parameters
    ----------
    alpha_K, beta_K : 1-based occupied MO index lists for bra determinant K
    alpha_L, beta_L : 1-based occupied MO index lists for ket determinant L
    S_mo            : (nmo x nmo) MO overlap matrix, 0-based indexing

    Returns
    -------
    overlap : complex scalar

    Notes
    -----
    Alpha and beta blocks are factorized for efficiency, but OpenMolcas
    determinants are defined in interleaved spin-orbital order
    ``(1a, 1b, 2a, 2b, ...)``. The permutation phases converting between that
    convention and grouped ``(all alpha, all beta)`` order are included here.
    Determinants with different numbers of alpha or beta electrons have zero
    overlap. Consequently, different M_S projections are exactly orthogonal in
    this spin-free collinear representation.
    """
    aK = [a - 1 for a in alpha_K]
    aL = [a - 1 for a in alpha_L]
    bK = [b - 1 for b in beta_K]
    bL = [b - 1 for b in beta_L]

    if len(aK) != len(aL) or len(bK) != len(bL):
        return 0.0 + 0.0j

    det_alpha = (
        1.0 + 0.0j if len(aK) == 0
        else np.linalg.det(S_mo[np.ix_(aK, aL)])
    )
    det_beta = (
        1.0 + 0.0j if len(bK) == 0
        else np.linalg.det(S_mo[np.ix_(bK, bL)])
    )

    # The factorized determinants above use grouped (all-alpha, all-beta)
    # ordering. OpenMolcas occupations use interleaved (1a,1b,2a,2b,...)
    # ordering, so restore the permutation phase for each determinant.
    phase_K = -1.0 if sum(b < a for b in beta_K for a in alpha_K) % 2 else 1.0
    phase_L = -1.0 if sum(b < a for b in beta_L for a in alpha_L) % 2 else 1.0
    return phase_K * phase_L * det_alpha * det_beta


def ci_overlap_general(data_bra, data_ket, S_mo, active_space, inactive_orbs, nstates,
    coeff_thresh=1e-6,
    verbose=False,
):
    """
    Compute the full MRCI overlap matrix <Psi_I(bra) | Psi_J(ket)>:

        S_CI[I,J] = sum_{K,L} C_I^K * C_J^L * det(S_alpha_KL) * det(S_beta_KL)

    Parameters
    ----------
    data_bra, data_ket : [energies, confs_list, CI_list]
        confs_list[istate] = list of occupation tuples
        CI_list[istate]    = list of CI coefficients
    S_mo               : (nmo x nmo) MO overlap matrix
                         Pass MO_prev.conj().T @ S_ao @ MO_curr for time-overlap
                         Pass MO_curr.conj().T @ S_ao @ MO_curr for same-time
    active_space       : list of 1-based MO indices (active orbitals)
    inactive_orbs      : list of 1-based MO indices (inactive orbitals)
    nstates            : number of electronic states
    coeff_thresh       : skip determinant pairs where |C_K * C_L| < threshold
    verbose            : print screening statistics

    Returns
    -------
    S_ci : np.ndarray, shape (nstates, nstates), dtype complex128

    Notes
    -----
    State ordering is exactly the outer-list ordering in ``data_bra`` and
    ``data_ket``. Data returned directly by :func:`read_ci_vectors` is ordered
    by spin-free root. Data returned by :func:`expand_ms_projections` is ordered
    first by root and then by descending M_S.

    In a spin-free calculation, matrix elements between different M_S values
    vanish because those states contain different alpha/beta electron counts.
    Equal-M_S components belonging to different roots need not have zero time
    overlap. The different-M_S rule does not apply after spin-orbit or other
    non-collinear interactions mix spin projections.
    """
    S_ci = np.zeros((nstates, nstates), dtype=np.complex128)

    cache_bra = {}
    cache_ket = {}

    for I in range(nstates):
        for K, det_K in enumerate(data_bra[1][I]):
            cache_bra[(I, K)] = build_full_orbital_lists(
                det_K, active_space, inactive_orbs
            )
    for J in range(nstates):
        for L, det_L in enumerate(data_ket[1][J]):
            cache_ket[(J, L)] = build_full_orbital_lists(
                det_L, active_space, inactive_orbs
            )

    total_pairs   = 0
    skipped_pairs = 0

    for I in range(nstates):
        coeffs_I = data_bra[2][I]
        for J in range(nstates):
            coeffs_J  = data_ket[2][J]
            overlap_IJ = 0.0 + 0.0j

            for K in range(len(data_bra[1][I])):
                c_K = coeffs_I[K]
                for L in range(len(data_ket[1][J])):
                    c_L = coeffs_J[L]
                    total_pairs += 1

                    if abs(c_K * c_L) < coeff_thresh:
                        skipped_pairs += 1
                        continue

                    aK, bK = cache_bra[(I, K)]
                    aL, bL = cache_ket[(J, L)]

                    overlap_IJ += c_K * c_L * slater_det_overlap(
                        aK, bK, aL, bL, S_mo
                    )

            S_ci[I, J] = overlap_IJ

    if verbose:
        pct = 100.0 * skipped_pairs / max(total_pairs, 1)
        print(
            f"[ci_overlap_general] Total det pairs : {total_pairs} | "
            f"Skipped (|CK*CL| < {coeff_thresh:.0e}): {skipped_pairs} ({pct:.1f}%)"
        )

    return S_ci


def _infer_nstates_from_ciroot(ciroot):
    """
    Extract the number of states from ciroot parameter.
    
    Args:
        ciroot: Number of roots in various formats
        
    Returns:
        int or None: The number of states to compute, or None if cannot infer
    """
    if ciroot is None:
        return None

    if isinstance(ciroot, int):
        return ciroot

    if isinstance(ciroot, (list, tuple)):
        return int(ciroot[0]) if len(ciroot) > 0 else None

    if isinstance(ciroot, str):
        parts = ciroot.replace(",", " ").split()
        ints = [int(x) for x in parts if x.isdigit()]
        return ints[0] if ints else None

    return None

class tmp:
    pass


def _compute_spin_manifold(coords, params, itraj, manifold, manifold_index,
                           multiple_manifolds):
    """Run and evaluate one fixed-multiplicity OpenMolcas manifold."""
    atom_labels = params["atom_labels"]
    timestep = params.get("timestep", 0)
    base_run_params = copy.deepcopy(params.get("molcas_run_params", {
        "basis": "ANO-RCC-VDZP",
        "charge": 0,
        "spin": 1,
        "nactel": "6 0 0",
        "inactive": 5,
        "ras2": 6,
        "ciroot": "2 2 1",
        "prwf": 0.0,
        "thre": "1.0e-10",
    }))
    molcas_run_params = copy.deepcopy(base_run_params)
    molcas_run_params.update(copy.deepcopy(manifold))

    multiplicity = int(molcas_run_params.get("spin", 1))
    nroots = int(manifold.get(
        "nroots", params.get("nroots", params.get(
            "nstates", _infer_nstates_from_ciroot(
                molcas_run_params.get("ciroot")) or 2))
    ))
    include_ms = bool(manifold.get(
        "include_ms_projections",
        params.get("include_ms_projections", multiple_manifolds),
    ))
    nstates = nroots * multiplicity if include_ms else nroots

    gradient_states, nac_pairs = _normalize_property_requests(
        nroots,
        manifold.get("gradient_states", params.get(
            "gradient_states", molcas_run_params.get("gradient_states"))),
        manifold.get("nac_pairs", params.get(
            "nac_pairs", molcas_run_params.get("nac_pairs"))),
    )
    if not nac_pairs and molcas_run_params.get("nac_states") is not None:
        legacy_i, legacy_j = molcas_run_params["nac_states"]
        nac_pairs = [(int(legacy_i) - 1, int(legacy_j) - 1)]
    molcas_run_params["nstates"] = nroots
    molcas_run_params["gradient_states"] = gradient_states
    molcas_run_params["nac_pairs"] = nac_pairs
    if "nac_nocsf" in params and "nac_nocsf" not in manifold:
        molcas_run_params["nac_nocsf"] = params["nac_nocsf"]

    prwf = float(str(molcas_run_params.get("prwf", 0.0)).replace(
        "D", "E").replace("d", "e"))
    if prwf > 0.0:
        warnings.warn(
            f"PRWF={prwf:g} truncates CSFs before determinant expansion; use "
            "PRWF=0.0 for rigorous post-expansion ci_threshold filtering.",
            UserWarning,
            stacklevel=3,
        )

    history_key = (itraj, manifold_index) if multiple_manifolds else itraj
    suffix = f"_spin{multiplicity}_manifold{manifold_index}" if multiple_manifolds else ""
    wd = f"{params.get('working_directory_prefix', 'wd')}_itraj{itraj}{suffix}"
    jobid = f"_timestep_{timestep}_traj_{itraj}{suffix}"
    run_params = {
        "atom_labels": atom_labels,
        "exe": params.get("exe", "pymolcas"),
        "molcas_run_params": molcas_run_params,
        "working_directory": wd,
        "molcas_jobid": jobid,
        "input_prefix": params.get("molcas_input_prefix", "input_"),
        "output_prefix": params.get("molcas_output_prefix", "output_"),
        "fileorb": params["rasorb_prev"].get(history_key),
    }
    out_path, rasorb_file = run_molcas(coords, run_params)
    info, MO_curr, data_curr = read_molcas_orbital_info({
        "filename": out_path,
        "rasorb_file": rasorb_file,
        "nstates": nroots,
        "h5_path": rasorb_file.replace(".RasOrb", ".rasscf.h5"),
        "ci_threshold": params.get(
            "ci_threshold", params.get("ci_coeff_thresh", 1e-6)),
    })
    if isinstance(data_curr, dict):
        energies = data_curr.get("energies", [])
        confs = data_curr.get("confs", [])
        coefficients = data_curr.get("CI", [])
    else:
        energies, confs, coefficients = data_curr
    if len(energies) < nroots:
        raise ValueError(
            f"Spin multiplicity {multiplicity}: found {len(energies)} energies, "
            f"expected {nroots}."
        )
    root_energies = list(energies[:nroots])
    check_energy_continuity(
        root_energies, params["energies_prev"].get(history_key),
        itraj=itraj, timestep=timestep,
        thresh_eV=params.get("energy_continuity_thresh_eV", 0.3),
    )

    ms_labels = [(root, multiplicity - 1) for root in range(nroots)]
    if include_ms:
        confs, coefficients, ms_labels = expand_ms_projections(
            confs, coefficients, multiplicity)
    energies = [root_energies[root] for root, _ in ms_labels]
    if not any(confs):
        raise ValueError(
            f"Spin multiplicity {multiplicity}: all CI states are empty; "
            "check PRWF, ci_threshold, and RASSCF convergence."
        )

    MO_curr = np.asarray(MO_curr, dtype=np.complex128)
    nbas = info.get("nao", MO_curr.shape[0])
    S_ao = np.asarray(read_ao_overlap(out_path, nbas), dtype=np.complex128)
    active_space = info.get(
        "actual_orbital_space",
        list(range(info["min_active"], info["max_active"] + 1)),
    )
    inactive_orbs = list(range(1, info["nocc"] + 1))
    for state_idx, state_confs in enumerate(confs):
        for conf_idx, conf in enumerate(state_confs):
            if len(conf) != len(active_space):
                raise ValueError(
                    f"Occupation tuple length ({len(conf)}) != active_space "
                    f"length ({len(active_space)}) for spin {multiplicity}, "
                    f"state {state_idx + 1}, configuration {conf_idx + 1}."
                )

    current_data = (energies, confs, coefficients)
    is_first_time = params["is_first_time"].get(history_key, True)
    if is_first_time:
        MO_prev = copy.deepcopy(MO_curr)
        data_prev = copy.deepcopy(current_data)
    else:
        MO_prev = params["MO_prev"].get(history_key, MO_curr)
        data_prev = params["data_prev"].get(history_key, current_data)
    st_mo = MO_prev.conj().T @ S_ao @ MO_curr
    s_mo = MO_curr.conj().T @ S_ao @ MO_curr
    overlap_kwargs = {
        "active_space": active_space,
        "inactive_orbs": inactive_orbs,
        "nstates": nstates,
        "coeff_thresh": params.get("overlap_pair_threshold", 0.0),
        "verbose": params.get("verbose", False),
    }
    st_ci = ci_overlap_general(data_prev, current_data, st_mo, **overlap_kwargs)
    s_ci = ci_overlap_general(current_data, current_data, s_mo, **overlap_kwargs)
    gradients, nac_vectors = read_alaska_properties(
        out_path, len(atom_labels), gradient_states=gradient_states,
        nac_pairs=nac_pairs, nstates=nroots,
    )

    params["MO_prev"][history_key] = copy.deepcopy(MO_curr)
    params["data_prev"][history_key] = copy.deepcopy(current_data)
    params["is_first_time"][history_key] = False
    params["rasorb_prev"][history_key] = rasorb_file
    params["energies_prev"][history_key] = root_energies
    return {
        "multiplicity": multiplicity,
        "nroots": nroots,
        "nstates": nstates,
        "energies": energies,
        "ms_labels": ms_labels,
        "time_overlap": st_ci,
        "overlap": s_ci,
        "gradients": gradients,
        "nac_vectors": nac_vectors,
    }

def molcas_compute_adi(q, params, full_id):
    """
    Perform a single-time-step electronic structure evaluation using OpenMolcas
    for trajectory-based nonadiabatic dynamics. Computes molecular orbitals (MOs),
    SA-CASSCF CI states, and their overlaps between consecutive time steps,
    constructing the adiabatic Hamiltonian, vibronic Hamiltonian, and derivative
    couplings.

    The function is designed for trajectory-based nonadiabatic methods, such as:
        - FSSH (Fewest Switches Surface Hopping)
        - Ehrenfest dynamics
        - Mapping-based methods
        - Exact factorization / quantum trajectory approaches

    Workflow
    --------
    1. Extract nuclear coordinates for the trajectory from `q`.
    2. Normalize either the legacy single-multiplicity input or a list of
       fixed-spin manifolds.
    3. For each manifold, write and run an independent OpenMolcas SA-CASSCF
       calculation in a trajectory- and spin-specific directory.
    4. Parse the output files:
        - `job.out` → CASSCF energies, CI vectors, orbital space metadata
        - `job.RasOrb` → MO coefficient matrix
    5. Build a determinant cache from CI vectors (truncated by `ci_coeff_thresh`).
    6. Compute CI overlaps between the previous and current time step using
       `ci_overlap_general()` (MO-transformed Slater determinant overlaps).
    7. Assemble:
        - Time-overlap matrix between consecutive CI states
        - Adiabatic Hamiltonian (diagonal = CASSCF state energies)
        - Vibronic Hamiltonian (including approximate derivative couplings)
    8. Update trajectory-specific previous-state data in `params`.

    Parameters
    ----------
    q : MATRIX
        Nuclear coordinates for all trajectories.
        Shape: (3 * N_atoms, N_trajectories)
        Units: Bohr
        Column `itraj` corresponds to trajectory `itraj`.

    params : dict
        Dictionary of simulation parameters and trajectory state information.
        Keys used include:

        **Required:**
        atom_labels : list of str
            Atomic symbols, e.g., ["O", "H", "H"].
        molcas_run_params : dict
            Parameters for OpenMolcas SA-CASSCF:
                - basis : str, e.g. "ANO-RCC-VDZP"
                - nactel : int, number of active electrons
                - ras2 : list[int], active orbital indices (1-based)
                - ciroot : list[list[int]], e.g. [[2,2],[2,1]] for state averaging
            See OpenMolcas documentation for all available keywords.
            These values are also used as common defaults for every entry in
            ``spin_manifolds``.

        **Optional / Internal (updated in-place):**
        dt : float, default=41.0
            Nuclear time step in atomic units (1 fs ≈ 41.341 a.u.).
        molcas_exe : str, default="pymolcas"
            Command to invoke OpenMolcas.
        working_directory_prefix : str, default="wd"
            Prefix for trajectory-specific directories.
        ci_threshold : float, default=1e-6
            Discard determinants whose final combined weight ``|C|^2`` is below
            this value. ``ci_coeff_thresh`` is accepted as a legacy alias.
        include_ms_projections : bool, default=False
            Represent every spin-free root by all ``2S+1`` distinct M_S
            projections, ordered root-major from M_S=S through M_S=-S.
            Here ``spin`` in ``molcas_run_params`` is the multiplicity ``2S+1``.
            It defaults to ``True`` when ``spin_manifolds`` is supplied, while
            retaining ``False`` for the legacy single-manifold calculation.
        spin_manifolds : list[dict], optional
            Fixed-spin calculations to combine into one Libra state space.
            Each dictionary must define a unique positive ``spin`` multiplicity,
            overrides ``molcas_run_params``, and may additionally define
            ``nroots``, ``include_ms_projections``, ``gradient_states``, and
            ``nac_pairs``. List order defines manifold order; roots and M_S
            projections define the ordering within each manifold. Omitting this
            key selects the backward-compatible single-multiplicity path.
        nroots : int, optional
            Number of spin-free CASSCF roots. If omitted, ``nstates`` is used,
            followed by the number inferred from ``ciroot``. When M_S expansion
            is enabled, the returned Libra dimension is
            ``nroots * spin_multiplicity``.
        overlap_pair_threshold : float, default=0.0
            Optional screening threshold for products of determinant
            coefficients inside the overlap double sum.
        gradient_states : "all", int, or iterable[int], optional
            Zero-based states whose analytical energy gradients are requested
            from ALASKA and returned in ``d1ham_adi``.
        nac_pairs : "all", pair[int, int], or iterable[pair[int, int]], optional
            Zero-based state pairs whose analytical spatial derivative-coupling
            vectors are requested from ALASKA and returned in ``dc1_adi``.
        nac_nocsf : bool, default=False
            Add ALASKA's NOCSF keyword to NAC calculations.
        scf_method : {"auto", "rhf", "uhf", "rohf"}, default="auto"
            Preliminary orbital method. ``auto`` selects UHF for multiplicity
            greater than one. OpenMolcas has no SCF-level ROHF implementation;
            ``rohf`` therefore proceeds directly to its spin-adapted RASSCF
            module, which is the documented OpenMolcas route for ROHF.
        uhf_orbital_set : {"alpha", "beta"}, default="beta"
            UHF orbital set supplied to the first RASSCF calculation.
        is_first_time : dict
            Dictionary keyed by trajectory index (`itraj`) with boolean values.
            True indicates that the current step is the first step of this trajectory.
        act_state : dict
            Dictionary keyed by trajectory index (`itraj`) with integer values
            indicating the active electronic state for this trajectory.
        MO_prev : dict
            Previous MO coefficients per trajectory (updated in-place).
            Shape: (nbas, nbas) per entry.
        data_prev : dict
            Previous CI data (energies, CI vectors) per trajectory (updated in-place).
        coordinates_prev : dict
            Previous nuclear coordinates per trajectory (updated in-place).
        verbose : bool, default=False
            Print detailed progress information during execution.
        timestep : int
            Current time-step index (used for file naming and diagnostics).

    full_id : int or object
        Encoded trajectory identifier (decoded to extract `itraj`).

    Returns
    -------
    obj : tmp (Libra temporary object)
        Object containing adiabatic electronic properties for this trajectory.
        Attributes include:

        ham_adi : CMATRIX (nstates, nstates)
            Adiabatic Hamiltonian matrix.
            Diagonal entries are SA-CASSCF state energies (in Hartree).

        hvib_adi : CMATRIX (nstates, nstates)
            Vibronic Hamiltonian including nonadiabatic coupling:
                Hvib_ij = E_i δ_ij - i d_ij
            where d_ij is the approximate derivative coupling.

        time_overlap_adi : CMATRIX (nstates, nstates)
            Time-overlap matrix S_ij(t, t+dt) = ⟨Ψ_i(t) | Ψ_j(t+dt)⟩.
            Computed via `ci_overlap_general()` using the determinant cache
            and MO-transformed Slater determinant overlaps.

        basis_transform : CMATRIX (nstates, nstates)
            Basis transformation matrix (currently set to identity).

        d1ham_adi : CMATRIXList
            One derivative Hamiltonian per nuclear degree of freedom. Requested
            energy gradients populate the corresponding diagonal elements.

        dc1_adi : CMATRIXList
            One spatial derivative-coupling matrix per nuclear degree of freedom.
            Requested ALASKA NAC vectors populate anti-Hermitian off-diagonals.

        ms_labels : list[tuple[int, int]]
            State labels ``(spin_free_root, 2*M_S)`` in the same ordering as
            every returned electronic matrix. Retained for compatibility; root
            numbers restart at zero in each manifold.

        spin_labels : list[tuple[int, int, int]]
            Unambiguous global labels
            ``(spin_multiplicity, spin_free_root, 2*M_S)``.

    Notes
    -----
    - All computations are performed in **trajectory-specific directories**
      to ensure thread safety when running multiple trajectories in parallel.
      Mixed-spin calculations add ``_spin<M>_manifold<K>`` to keep their files
      and orbital restarts independent.
    - The SA-CASSCF calculation uses a **state-averaged** formalism; energies
      are printed for each root included in the averaging.
    - Without ``include_ms_projections``, each OpenMolcas root appears once and
      ``ms_labels`` identifies its highest-weight component. With the option
      enabled, every root is expanded into all ``2S+1`` projections. Thus two
      doublet roots give ``[(0, 1), (0, -1), (1, 1), (1, -1)]``, while two
      quartet roots give twelve states ordered as root 0 projections
      ``3/2, 1/2, -1/2, -3/2`` followed by the same ordering for root 1.
    - Spin-free energies and gradients are identical for every M_S component
      of a root. They are therefore repeated on the corresponding diagonal
      entries. A root-to-root ALASKA NAC is copied only between components with
      equal M_S; different-M_S elements remain zero.
    - Different M_S components have zero same-time and time overlap in the
      present spin-free, collinear formalism. Spin-orbit-coupled states require
      a spin-mixed representation and do not obey this block structure.
    - Different multiplicities are assembled as exact zero-coupled blocks for
      the Hamiltonian, overlaps, derivative couplings, and time couplings.
      Therefore this basis is ready to receive later SOC matrix elements, but
      it cannot produce intersystem crossing until such spin-dependent terms
      are supplied.
    - CI overlaps are computed using the method of Plasser et al. (JCP 2016):
      the Slater determinant overlap is factorised into MO overlap contributions,
      and the CI overlap is assembled as:

        ⟨Ψ_I | Ψ_J⟩ = Σ_{pq} C_I^p * C_J^q * det( MO_prev^T * S_AO * MO_curr )

      where S_AO is the atomic orbital overlap matrix (approximated as identity
      in the MO basis following an orthonormalisation step).

    - Derivative couplings are **approximated** from the anti-symmetric part of
      the time-overlap matrix divided by 2 dt:

          d_ij ≈ [ S_ij(t, t+dt) - S_ji(t, t+dt) ] / (2 * dt)

      This is the **finite-difference overlap-based** approximation (the
      "Hammes-Schiffer–Tully" approach), valid when dt is small.

    - Energies are in **Hartree**, time in **atomic units**, coordinates in
      **Bohr**, and overlaps are **dimensionless**.

    - `is_first_time` and `act_state` are dictionaries keyed by trajectory index,
      enabling simultaneous tracking of multiple trajectories.

    - If MO/CI data is missing for the previous step (first timestep), the
      overlap is set to the identity matrix.

    Example
    -------
    >>> params = {
    ...     "atom_labels": ["O", "H", "H"],
    ...     "molcas_run_params": {
    ...         "basis": "ANO-RCC-VDZP",
    ...         "nactel": 4,
    ...         "ras2": [2, 3, 4, 5, 6, 7],
    ...         "ciroot": [[2, 2], [2, 1]],
    ...     },
    ...     "dt": 41.0,
    ...     "ci_coeff_thresh": 0.001,
    ...     "verbose": True,
    ... }
    >>> obj = molcas_compute_adi(q, params, full_id)
    >>> print(obj.ham_adi)
    >>> print(obj.time_overlap_adi)

    A triplet calculation with two spin-free roots can request explicit spin
    projections as follows::

        params["nroots"] = 2
        params["molcas_run_params"]["spin"] = 3
        params["include_ms_projections"] = True
        obj = molcas_compute_adi(q, params, full_id)
        # obj.ham_adi is 6 x 6
        # obj.ms_labels == [(0, 2), (0, 0), (0, -2),
        #                   (1, 2), (1, 0), (1, -2)]

    One singlet root and one triplet root form a four-state spin-diabatic basis
    using independent OpenMolcas calculations::

        params["spin_manifolds"] = [
            {"spin": 1, "nroots": 1, "ciroot": "1 1 1"},
            {"spin": 3, "nroots": 1, "ciroot": "1 1 1"},
        ]
        obj = molcas_compute_adi(q, params, full_id)
        # obj.spin_labels == [(1, 0, 0),
        #                     (3, 0, 2), (3, 0, 0), (3, 0, -2)]
        # All singlet-triplet matrix blocks are zero without SOC.
    """
    Id = Cpp2Py(full_id)
    itraj = Id[-1]
    coords = q.col(itraj)
    for name in ("MO_prev", "data_prev", "is_first_time", "rasorb_prev",
                 "energies_prev"):
        params.setdefault(name, {})

    requested_manifolds = params.get("spin_manifolds")
    multiple_manifolds = requested_manifolds is not None
    if multiple_manifolds:
        if not isinstance(requested_manifolds, (list, tuple)) or not requested_manifolds:
            raise ValueError("spin_manifolds must be a non-empty list of dictionaries")
        if any(not isinstance(item, dict) for item in requested_manifolds):
            raise TypeError("Each spin_manifolds entry must be a dictionary")
        manifolds = list(requested_manifolds)
        if any("spin" not in item for item in manifolds):
            raise ValueError("Each spin_manifolds entry must define spin multiplicity")
        multiplicities = [int(item["spin"]) for item in manifolds]
        if any(value < 1 for value in multiplicities):
            raise ValueError("Spin multiplicities must be positive integers")
        if len(set(multiplicities)) != len(multiplicities):
            raise ValueError("spin_manifolds must not repeat a spin multiplicity")
    else:
        manifolds = [{}]

    results = [
        _compute_spin_manifold(
            coords, params, itraj, manifold, index, multiple_manifolds
        )
        for index, manifold in enumerate(manifolds)
    ]
    nstates = sum(result["nstates"] for result in results)
    dt = params.get("dt", 1.0 * units.fs2au)
    atom_labels = params["atom_labels"]

    obj = SimpleNamespace()
    obj.spin_labels = [
        (result["multiplicity"], root, twice_ms)
        for result in results for root, twice_ms in result["ms_labels"]
    ]
    obj.ms_labels = [label[1:] for label in obj.spin_labels]
    obj.ham_adi          = CMATRIX(nstates, nstates)
    obj.nac_adi          = CMATRIX(nstates, nstates)
    obj.hvib_adi         = CMATRIX(nstates, nstates)
    obj.time_overlap_adi = CMATRIX(nstates, nstates)
    obj.overlap_adi      = CMATRIX(nstates, nstates)
    obj.basis_transform  = CMATRIX(nstates, nstates)
    obj.d1ham_adi        = CMATRIXList()
    obj.dc1_adi           = CMATRIXList()

    for _ in range(3 * len(atom_labels)):
        obj.d1ham_adi.append(CMATRIX(nstates, nstates))
        obj.dc1_adi.append(CMATRIX(nstates, nstates))

    offset = 0
    for result in results:
        local_nstates = result["nstates"]
        for i in range(local_nstates):
            gi = offset + i
            energy = complex(result["energies"][i])
            obj.ham_adi.set(gi, gi, energy)
            obj.hvib_adi.set(gi, gi, energy)
            obj.basis_transform.set(gi, gi, 1.0 + 0.0j)
            for j in range(local_nstates):
                gj = offset + j
                obj.time_overlap_adi.set(
                    gi, gj, complex(result["time_overlap"][i, j]))
                obj.overlap_adi.set(gi, gj, complex(result["overlap"][i, j]))

        for root, gradient in result["gradients"].items():
            local_states = [
                state for state, (state_root, _) in enumerate(result["ms_labels"])
                if state_root == root
            ]
            for atom in range(len(atom_labels)):
                for xyz in range(3):
                    for state in local_states:
                        obj.d1ham_adi[3 * atom + xyz].set(
                            offset + state, offset + state,
                            complex(gradient[atom, xyz]),
                        )

        # ALASKA NAC=i,j prints <Psi_j|nabla Psi_i>. Only equal-M_S
        # components within this fixed-spin manifold are coupled.
        for (root_i, root_j), vector in result["nac_vectors"].items():
            local_pairs = [
                (i, j) for i, (ri, msi) in enumerate(result["ms_labels"])
                for j, (rj, msj) in enumerate(result["ms_labels"])
                if ri == root_i and rj == root_j and msi == msj
            ]
            for atom in range(len(atom_labels)):
                for xyz in range(3):
                    value = complex(vector[atom, xyz])
                    for i, j in local_pairs:
                        gi, gj = offset + i, offset + j
                        obj.dc1_adi[3 * atom + xyz].set(gj, gi, value)
                        obj.dc1_adi[3 * atom + xyz].set(
                            gi, gj, -value.conjugate())
        offset += local_nstates

    for i in range(nstates):
        for j in range(i + 1, nstates):
            dij = (obj.time_overlap_adi.get(i, j) - obj.time_overlap_adi.get(j, i)) / (2.0 * dt)
            obj.nac_adi.set(i, j,  dij)
            obj.nac_adi.set(j, i, -dij.conjugate())
            obj.hvib_adi.set(i, j, -1.0j * dij)
            obj.hvib_adi.set(j, i,  1.0j * dij.conjugate())

    return obj
