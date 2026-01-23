# *********************************************************************************
# * Copyright (C) 2025-2026 Alexey V. Akimov
# *
# * This file is distributed under the terms of the GNU General Public License
# * as published by the Free Software Foundation, either version 3 of
# * the License, or (at your option) any later version.
# * See the file LICENSE in the root directory of this distribution
# * or <http://www.gnu.org/licenses/>.
# *
# *********************************************************************************/

"""
.. module:: interfaces
   :platform: Unix, Windows
   :synopsis: this module implements functions for using the citools utilities in 
   connection with particular codes
.. moduleauthor:: Alexey V. Akimov, Kevin Walsh, ChatGPT

"""

from collections import Counter
from typing import List, Tuple, Iterable, Dict, Any
import numpy as np
from scipy.sparse import coo_matrix

from . import csf
from . import slatdet as sd

#from liblibra_core import MATRIX, CMATRIX

def find_matches(
    config: List[int],
    possible_configs: List[Tuple[List[int], float]],
    max_unpaired: int
) -> List[Tuple[List[int], float]]:
    """
    Find configurations matching a reference Slater determinant in orbital composition
    and order, while allowing different spin sign assignments.

    The `max_unpaired` parameter here represents **2 × Mₛ**, the total spin polarization:
    - Singlet: 2Mₛ = 0
    - Doublet: 2Mₛ = ±1
    - Triplet: 2Mₛ = ±2
    etc.

    Parameters
    ----------
    config : list[int]
        Reference configuration (single determinant) in Libra format,
        where positive orbital indices correspond to α-spin and negative to β-spin.
        Example: [6, -6, 7, -7, 9, -8]
    possible_configs : list[tuple[list[int], float]]
        List of all possible consistently ordered Slater determinants (SDs)
        for the given active space, represented as tuples `(det, weight)` or `(det, coeff)`.
    max_unpaired : int
        Target 2Mₛ spin polarization (sum of spin signs).
        For singlet configurations use 0.

    Returns
    -------
    list[tuple[list[int], float]]
        A filtered list of configurations from `possible_configs` that:
        - Have the same absolute orbital composition as `config`
        - Differ only by spin sign assignments
        - Have total spin polarization `Σ sign(orbital) == max_unpaired`

    Examples
    --------
    >>> ref = [6, -6, 7, -7, 9, -8]
    >>> configs = [
    ...     ([6, -6, 7, -7, 8, 9], 1),
    ...     ([6, -6, 7, -7, 8, -9], 1),
    ...     ([6, -6, 7, -7, -8, 9], 1),
    ...     ([6, -6, 7, -7, -8, -9], 1)
    ... ]
    >>> find_matches(ref, configs, max_unpaired=0)
    [([6, -6, 7, -7, 8, -9], 1), ([6, -6, 7, -7, -8, 9], 1)]
    """

    # Reference absolute orbital composition
    ref_counter = Counter(abs(x) for x in config)

    # Keep configurations with same orbital composition and given spin polarization
    matches = [
        entry for entry in possible_configs
        if Counter(abs(x) for x in entry[0]) == ref_counter
        and int(np.sum(np.sign(entry[0]))) == max_unpaired
    ]

    return matches


def build_minimal_csf_basis(
    configs: List[Tuple[int, ...]],
    active_space: List[int],
    nelec: int,
    max_unpaired: int
) -> Tuple[List[Tuple[Any, ...]], List[Tuple[Any, ...]]]:
    """
    Construct a minimal CSF (Configuration State Function) basis
    from raw MOPAC/Libra configurations.

    Parameters
    ----------
    configs : list[tuple[int]]
        Configurations extracted from Libra or MOPAC output (raw orbital numbers).
    active_space : list[int]
        List of active orbital numbers defining the active space.
    nelec : int
        Total number of electrons.
    max_unpaired : int
        Twice the target spin projection (2*Ms).
        For example, 0 corresponds to singlet-only configurations.

    Returns
    -------
    min_basis : list[tuple]
        List of matching determinants (CSFs) as (configuration, phase) tuples.
    transposed_basis : list[tuple]
        Transposed representation, i.e. `list(zip(*min_basis))`.
        Each element groups together all configurations or all phases across the basis.

    Notes
    -----
    - Each entry in `min_basis` corresponds to a determinant that matches
      one of the given MOPAC configurations, within the spin constraint `max_unpaired`.
    - `transposed_basis` is convenient for splitting the basis into
      separate lists of determinants and phases.
    - If no matches are found, `transposed_basis` is returned as an empty list.

    Example
    -------
    >>> configs0_raw = [(6, -6, 7, -7, 9, -8)]
    >>> active_space = [6, 7, 8, 9, 10, 11]
    >>> nelec = 6
    >>> max_unpaired = 0   # singlet configurations only
    >>>
    >>> min_basis, (all_confs, all_phases) = build_minimal_csf_basis(
    ...     configs0_raw, active_space, nelec, max_unpaired
    ... )
    >>>
    >>> print(min_basis)
    [((6, -6, 7, -7, 8, -9), 1),
     ((6, -6, 7, -7, -8, 9), 1)]
    """

    # Generate all possible determinants (with parity) for the given active space
    dets: List[Tuple[Any, ...]] = list(sd.generate_determinants_with_parity(active_space, nelec))

    min_basis: List[Tuple[Any, ...]] = []
    for config in configs:
        matches = find_matches(config, dets, max_unpaired)
        min_basis.extend(matches)

    # Avoid error when basis is empty
    transposed_basis: List[Tuple[Any, ...]] = list(zip(*min_basis)) if min_basis else []

    return min_basis, transposed_basis

def build_minimal_csf_basis_singlet(
    configs: List[Tuple[int, ...]],
    active_space: List[int],
    nelec: int,
    max_unpaired: int
) -> Tuple[List[Tuple[Any, ...]], List[Tuple[Any, ...]]]:
    """
    Construct a minimal CSF (Configuration State Function) basis
    from raw MOPAC/Libra configurations.

    Parameters
    ----------
    configs : list[tuple[int]]
        Configurations extracted from Libra or MOPAC output (raw orbital numbers).
    active_space : list[int]
        List of active orbital numbers defining the active space.
    nelec : int
        Total number of electrons.
    max_unpaired : int
        Twice the target spin projection (2*Ms).
        For example, 0 corresponds to singlet-only configurations.

    Returns
    -------
    min_basis : list[tuple]
        List of matching determinants (CSFs) as (configuration, phase) tuples.
    transposed_basis : list[tuple]
        Transposed representation, i.e. `list(zip(*min_basis))`.
        Each element groups together all configurations or all phases across the basis.

    Notes
    -----
    - Each entry in `min_basis` corresponds to a determinant that matches
      one of the given MOPAC configurations, within the spin constraint `max_unpaired`.
    - `transposed_basis` is convenient for splitting the basis into
      separate lists of determinants and phases.
    - If no matches are found, `transposed_basis` is returned as an empty list.

    Example
    -------
    >>> configs0_raw = [(6, -6, 7, -7, 9, -8)]
    >>> active_space = [6, 7, 8, 9, 10, 11]
    >>> nelec = 6
    >>> max_unpaired = 0   # singlet configurations only
    >>>
    >>> min_basis, (all_confs, all_phases) = build_minimal_csf_basis(
    ...     configs0_raw, active_space, nelec, max_unpaired
    ... )
    >>>
    >>> print(min_basis)
    [((6, -6, 7, -7, 8, -9), 1),
     ((6, -6, 7, -7, -8, 9), 1)]
    """

    # Generate all possible determinants (with parity) for the given active space
    dets: List[Tuple[Any, ...]] = list(sd.generate_single_excitations(active_space, nelec))

    min_basis: List[Tuple[Any, ...]] = []
    for config in configs:
        matches = find_matches(config, dets, max_unpaired)
        min_basis.extend(matches)

    # Avoid error when basis is empty
    transposed_basis: List[Tuple[Any, ...]] = list(zip(*min_basis)) if min_basis else []

    return min_basis, transposed_basis



def map_to_active_indices(
    configs: List[Tuple[int, ...]],
    active_space: List[int]
) -> List[Tuple[int, ...]]:
    """
    Map orbital configurations from absolute orbital numbers to
    indices relative to the given active space (starting from 1).

    The signs of the orbitals (±) are preserved, meaning that
    positive/negative indices still indicate α/β spins.

    Parameters
    ----------
    configs : list[tuple[int, ...]]
        List of configurations, where each configuration is a tuple
        of signed orbital numbers (e.g., (6, -6, 7, -7, 9, -8)).
        Positive values correspond to α-spin, and negative to β-spin.
    active_space : list[int]
        The ordered list of orbital numbers defining the active space.
        Example: [6, 7, 8, 9, 10, 11]

    Returns
    -------
    list[tuple[int, ...]]
        List of configurations with each orbital mapped to its
        corresponding index in the active space (starting from 1).
        Signs are preserved.

    Raises
    ------
    ValueError
        If any orbital in `configs` is not present in the `active_space`.

    Example
    -------
    >>> configs = [
    ...     (6, -6, 7, -7, 8, -8),
    ...     (6, -6, 7, -7, 9, -8)
    ... ]
    >>> active_space = [6, 7, 8, 9, 10, 11]
    >>> map_to_active_indices(configs, active_space)
    [(1, -1, 2, -2, 3, -3), (1, -1, 2, -2, 4, -3)]

    Notes
    -----
    - The first orbital in `active_space` is assigned index 1, the second index 2, etc.
    - This function is useful for transforming determinants into active-space-based
      indexing schemes used in spin-adapted CSFs or CI vector construction.
    """

    # Precompute mapping: orbital_number → 1-based index in active space
    index_map = {orb: i + 1 for i, orb in enumerate(active_space)}

    mapped_configs: List[Tuple[int, ...]] = []

    for conf in configs:
        try:
            mapped_conf = tuple(int(np.sign(orb)) * index_map[abs(orb)] for orb in conf)
        except KeyError as e:
            raise ValueError(f"Orbital {e.args[0]} not found in active space {active_space}")
        mapped_configs.append(mapped_conf)

    return mapped_configs



def conf2csf_matrix(
    min_basis: List[Tuple[Tuple[int, ...], int]],
    all_confs: List[Tuple[int, ...]],
    all_phases: List[int],
    S: int = 0,
    Ms: int = 0
) -> coo_matrix:
    """
    Build the configuration-to-CSF (Configuration State Function) transformation matrix.

    This matrix `T` transforms the determinant (Slater determinant) basis
    into the spin-adapted CSF basis for a given total spin `S` and projection `Ms`.

    Parameters
    ----------
    min_basis : list[tuple[tuple[int], int]]
        Minimal determinant basis (e.g., output from `build_minimal_csf_basis`),
        where each element is a tuple `(configuration, phase)`.
    all_confs : list[tuple[int]]
        List of configurations (determinants) corresponding to rows of the matrix.
        Typically obtained as `all_confs, all_phases = zip(*min_basis)`.
    all_phases : list[int]
        Parity or phase factors corresponding to each determinant.
    S : int, optional
        Total spin quantum number (default is 0, singlet).
    Ms : int, optional
        Spin projection quantum number (default is 0).

    Returns
    -------
    T_sparse : scipy.sparse.coo_matrix
        Sparse configuration-to-CSF transformation matrix.
        Rows correspond to determinants in `all_confs`,
        columns correspond to CSFs in `min_basis`.

    Notes
    -----
    - The matrix `T` satisfies the relation:
      \[
      |\text{CSF}_j\rangle = \sum_i T_{ij} |\text{Det}_i\rangle
      \]
    - Each column corresponds to a spin-adapted CSF with given `(S, Ms)`.
    - The function assumes `csf.generate_CSFs_grouped()` produces
      a dictionary `csfs[(S, Ms)]`, where each entry is a list of
      `(determinant, coefficient)` pairs.

    Example
    -------
    >>> # 1. Build minimal basis from MOPAC/Libra configs
    >>> min_basis, (all_confs, all_phases) = build_minimal_csf_basis(
    ...     configs0_raw, active_space=[6,7,8,9,10,11], nelec=6, max_unpaired=0
    ... )
    >>>
    >>> # 2. Build configuration–CSF transformation matrix for singlet states
    >>> T = conf2csf_matrix(min_basis, all_confs, all_phases, S=0, Ms=0)
    >>>
    >>> print("Shape of T:", T.num_of_rows, "x", T.num_of_cols)
    >>> print("First few elements:")
    >>> T.show_matrix()
    """

    # Generate CSFs for the given basis, grouped by (S, Ms)
    csfs: Dict[Tuple[int, int], List[List[Tuple[List[int], float]]]] = \
        csf.generate_CSFs_grouped(min_basis)

    #csf.print_csfs(csfs)  # optional diagnostic output

    nconfigs = len(min_basis)
    if (S, Ms) not in csfs:
        raise ValueError(f"No CSFs found for (S={S}, Ms={Ms}).")

    ncsfs = len(csfs[(S, Ms)])  # number of CSFs for the given spin manifold

    """
    # Initialize complex transformation matrix
    T = CMATRIX(nconfigs, ncsfs)

    # Populate the matrix
    for j, csf_group in enumerate(csfs[(S, Ms)]):
        for det, coeff in csf_group:
            det_tuple = tuple(det)
            if det_tuple not in all_confs:
                raise ValueError(f"Determinant {det_tuple} not found in all_confs.")
            i = all_confs.index(det_tuple)
            T.set(i, j, coeff * all_phases[i] * (1.0 + 0j))
    """

    # List to store non-zero entries
    rows, cols, data = [], [], []

    # Build a dict for faster lookup
    conf_index = {tuple(conf): i for i, conf in enumerate(all_confs)}

    # Populate sparse T
    for j, csf_group in enumerate(csfs[(S, Ms)]):
        for det, coeff in csf_group:
            det_tuple = tuple(det)
            i = conf_index.get(det_tuple)
            if i is None:
                raise ValueError(f"Determinant {det_tuple} not found in all_confs.")
            rows.append(i)
            cols.append(j)
            data.append(coeff * all_phases[i] * (1.0 + 0j))

    # Convert to a scipy sparse matrix (CSR)
    T_sparse = coo_matrix((data, (rows, cols)), shape=(nconfigs, ncsfs), dtype=np.complex128)

    return T_sparse




def configs_and_T_matrix(
    configs0_raw: List[Tuple[int, ...]],
    active_space: List[int],
    orbital_space: List[int],
    nelec: int,
    S: int,
    Ms: int
) -> Tuple[List[Tuple[int, ...]], coo_matrix]:
    """
    Generate the minimal active-space configurations mapped to a given orbital space
    and the configuration-to-CSF transformation matrix for a CAS with given spin.

    Parameters
    ----------
    configs0_raw : list[tuple[int]]
        List of raw configurations from Libra/MOPAC (signed orbital indices).
    active_space : list[int]
        Orbitals defining the active space used to generate the minimal determinant basis.
    orbital_space : list[int]
        Orbital indices used for mapping configurations (output will be relative to this space).
    nelec : int
        Number of active electrons.
    S : int
        Total spin quantum number.
    Ms : int
        Spin projection quantum number.

    Returns
    -------
    mapped_basis : list[tuple[int, ...]]
        List of minimal configurations mapped to the specified `orbital_space`,
        with signs preserved (positive = α-spin, negative = β-spin).
    T : coo_matrix
        Complex-valued configuration-to-CSF transformation matrix.

    Example
    -------
    >>> # Build configurations and T matrix for singlet CAS
    >>> mapped_basis, T = configs_and_T_matrix(
    ...     configs0_raw,
    ...     active_space=[6, 7, 8, 9, 10, 11],
    ...     orbital_space=[2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
    ...     nelec=6,
    ...     S=0,
    ...     Ms=0
    ... )
    >>> # mapped_basis shows minimal determinants mapped to orbital_space indices
    >>> print(mapped_basis[:5])
    [(5, -5, 6, -6, 7, -7),
     (5, -5, 6, -6, 7, -8),
     (5, -5, 6, -6, -7, 8),
     (5, -5, 6, -6, 7, -9),
     (5, -5, 6, -6, -7, 9)]
    >>> # T is the configuration-to-CSF transformation matrix
    >>> print("Shape of T:", T.num_of_rows, "x", T.num_of_cols)
    """

    # 1. Build minimal SD basis with spin constraint 2*Ms
    min_basis, (all_confs, all_phases) = build_minimal_csf_basis(
        configs0_raw, active_space, nelec, 2*Ms
    )

    # 2. Map configurations to the specified orbital space
    mapped_basis = map_to_active_indices(all_confs, orbital_space)

    # 3. Compute configuration-to-CSF transformation matrix
    T = conf2csf_matrix(min_basis, all_confs, all_phases, S, Ms)

    return mapped_basis, T


def configs_and_T_matrix_singlet(
    configs0_raw: List[Tuple[int, ...]],
    active_space: List[int],
    orbital_space: List[int],
    nelec: int,
    S: int,
    Ms: int
) -> Tuple[List[Tuple[int, ...]], coo_matrix]:
    """
    Generate the minimal active-space configurations mapped to a given orbital space
    and the configuration-to-CSF transformation matrix for a CAS with given spin.

    Parameters
    ----------
    configs0_raw : list[tuple[int]]
        List of raw configurations from Libra/MOPAC (signed orbital indices).
    active_space : list[int]
        Orbitals defining the active space used to generate the minimal determinant basis.
    orbital_space : list[int]
        Orbital indices used for mapping configurations (output will be relative to this space).
    nelec : int
        Number of active electrons.
    S : int
        Total spin quantum number.
    Ms : int
        Spin projection quantum number.

    Returns
    -------
    mapped_basis : list[tuple[int, ...]]
        List of minimal configurations mapped to the specified `orbital_space`,
        with signs preserved (positive = α-spin, negative = β-spin).
    T : coo_matrix
        Complex-valued configuration-to-CSF transformation matrix.

    Example
    -------
    >>> # Build configurations and T matrix for singlet CAS
    >>> mapped_basis, T = configs_and_T_matrix(
    ...     configs0_raw,
    ...     active_space=[6, 7, 8, 9, 10, 11],
    ...     orbital_space=[2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
    ...     nelec=6,
    ...     S=0,
    ...     Ms=0
    ... )
    >>> # mapped_basis shows minimal determinants mapped to orbital_space indices
    >>> print(mapped_basis[:5])
    [(5, -5, 6, -6, 7, -7),
     (5, -5, 6, -6, 7, -8),
     (5, -5, 6, -6, -7, 8),
     (5, -5, 6, -6, 7, -9),
     (5, -5, 6, -6, -7, 9)]
    >>> # T is the configuration-to-CSF transformation matrix
    >>> print("Shape of T:", T.num_of_rows, "x", T.num_of_cols)
    """

    # 1. Build minimal SD basis with spin constraint 2*Ms
    min_basis, (all_confs, all_phases) = build_minimal_csf_basis_singlet(
        configs0_raw, active_space, nelec, 2*Ms
    )

    # 2. Map configurations to the specified orbital space
    mapped_basis = map_to_active_indices(all_confs, orbital_space)

    # 3. Compute configuration-to-CSF transformation matrix
    T = conf2csf_matrix(min_basis, all_confs, all_phases, S, Ms)

    return mapped_basis, T




def freeze(x):
    """
    Recursively convert a nested list structure into an immutable,
    hashable tuple representation.

    This function is intended to create a stable key for deep value
    comparison of configurations that may contain arbitrarily nested
    Python lists.

    Parameters
    ----------
    x : list or any
        A (possibly nested) list structure. Non-list elements are
        returned unchanged.

    Returns
    -------
    tuple or any
        An immutable representation of `x` suitable for use as a key
        in dictionaries or sets.

    Notes
    -----
    - Only `list` containers are converted recursively.
    - Elements are assumed to be hashable once lists are removed.
    - The original structure is *not* modified.
    """
    return tuple(freeze(i) if isinstance(i, list) else i for i in x)


def unique_confs(configs1, configs2, n_excited_states):
    """
    Collect unique configurations from two sets of *excited-state*
    configuration lists using deep value comparison, while preserving
    first-occurrence order.

    The expected input hierarchy is::

        configsX : list (over excited electronic states)
          └── ex : list (configurations in a given excited state)
                └── sd : list (possibly nested; e.g., Slater determinant)

    Only the first `n_excited_states` excited states from each input
    are considered. The ground state is *not* included and is assumed
    to be handled separately.

    Uniqueness is defined by the *deep contents* of each configuration
    (`sd`), not by object identity. Internally, configurations are
    converted into immutable tuple representations using `freeze()`
    to allow efficient hashing and comparison.

    Parameters
    ----------
    configs1, configs2 : list of list of list
        State-resolved configuration lists for excited states.
        `configsX[i]` contains the configurations contributing to
        excited state `i + 1`.

    n_excited_states : int
        Number of excited electronic states to include (starting from
        index 0 of `configsX`).

    Returns
    -------
    list of list
        A list of unique configurations (`sd`) in the order of their
        first appearance across
        `configs1[:n_excited_states]` followed by
        `configs2[:n_excited_states]`.

    Raises
    ------
    ValueError
        If `n_excited_states` is non-positive.

    Notes
    -----
    - The ground-state reference determinant is intentionally excluded.
    - Order is preserved across both inputs.
    - Time complexity is O(N) with respect to the total number of
      configurations examined.
    - This function is suitable for CI/SD-style workflows where
      configurations (e.g., Slater determinants) may be deeply nested.
    - The returned configurations are the *original objects*, not their
      frozen representations.
    """
    if n_excited_states <= 0:
        raise ValueError("n_excited_states must be a positive integer")

    unique_sds = []
    seen = set()

    # Iterate over the first n_excited_states excited states
    for configs in (
        configs1[:n_excited_states],
        configs2[:n_excited_states],
    ):
        for ex in configs:
            for sd in ex:
                key = freeze(sd)
                if key not in seen:
                    seen.add(key)
                    unique_sds.append(sd)

    return unique_sds



def ci_amplitudes_mtx(nstates, common_sd_basis, configs, ci_amplitudes):
    """
    Construct the CI coefficient matrix in a common Slater-determinant basis.

    The resulting matrix C has the structure::

        C[sd_index, state_index]

    where:
      - state_index = 0 corresponds to the ground state
      - state_index >= 1 corresponds to excited states
      - sd_index = 0 corresponds to the reference determinant
      - sd_index >= 1 corresponds to determinants in `common_sd_basis`

    Parameters
    ----------
    nstates : int
        Total number of electronic states, including the ground state.

    common_sd_basis : list
        List of unique Slater determinants forming the global CI basis,
        excluding the reference determinant. The ordering of this list
        defines the row ordering of the CI matrix.

    configs : list of list
        `configs[i]` contains the Slater determinants contributing to
        excited state `i+1`.

    ci_amplitudes : list of list
        `ci_amplitudes[i][j]` is the CI coefficient corresponding to
        `configs[i][j]` in excited state `i+1`.

    Returns
    -------
    numpy.ndarray
        CI coefficient matrix of shape `(nsd, nstates)`, where::

            nsd = len(common_sd_basis) + 1

        The matrix uses extended precision (`np.float128`) and is
        initialized such that the ground state is the reference
        determinant with unit coefficient.

    Notes
    -----
    - The reference determinant is assumed to contribute only to the
      ground state with coefficient 1.0.
    - Excited-state CI expansions do not include the reference
      determinant explicitly.
    - This function preserves the ordering of `common_sd_basis`.
    """
    # Number of Slater determinants including the reference determinant
    nsd = len(common_sd_basis) + 1

    # Allocate CI matrix
    C = np.zeros((nsd, nstates), dtype=np.float128)

    # Ground state = reference determinant
    C[0, 0] = 1.0

    # Precompute determinant -> global index mapping (O(1) lookup)
    #sd_to_index = {sd: i for i, sd in enumerate(common_sd_basis)}

    # Build hashable determinant → index mapping
    sd_to_index = {
        freeze(sd): i for i, sd in enumerate(common_sd_basis)
    }


    # Fill excited-state CI coefficients
    for istate in range(nstates - 1):
        for sd, amp in zip(configs[istate], ci_amplitudes[istate]):
            key = freeze(sd)
            j_glob = sd_to_index[key]
            C[j_glob + 1, istate + 1] = amp

    return C




def sd_and_csf_overlaps_singlet(
    st_mo,
    lowest_orbital,
    highest_orbital,
    nelec,
    homo_indx,
    common_sd_basis,
    _active_space=None,
    S=0,
    Ms=0,
    max_unpaired=0,
):
    """
    Compute Slater-determinant (SD) and configuration-state-function (CSF)
    overlap matrices for singlet excitations using molecular-orbital
    time-overlaps.

    This function constructs a reference determinant and a set of singly
    excited determinants defined by `common_sd_basis`, maps them into a
    spin-adapted singlet CSF basis, and computes both SD and CSF overlap
    matrices.

    Parameters
    ----------
    st_mo : sparse matrix or array-like, shape (2*norb, 2*norb)
        Molecular-orbital time-overlap matrix in the spin–orbital basis.
        Spin-up and spin-down blocks are assumed to be included explicitly.

    lowest_orbital : int
        Lowest spatial orbital index (1-based) included in `st_mo`.

    highest_orbital : int
        Highest spatial orbital index (1-based) included in `st_mo`.

    nelec : int
        Total number of electrons.

    homo_indx : int
        Index of the highest occupied molecular orbital (HOMO), 1-based.

    common_sd_basis : list
        List of excitation descriptors defining the singly excited
        determinants. Each element is assumed to encode an excitation
        `(occ, vir)` in spin–orbital notation.

    _active_space : iterable of int, optional
        Spatial orbital indices (1-based) defining the active space for
        spin adaptation. If None, all available orbitals are used.

    S : int, optional
        Total spin quantum number (default: 0, singlet).

    Ms : int, optional
        Spin projection quantum number (default: 0).

    max_unpaired : int, optional
        Maximum number of unpaired electrons allowed (currently unused;
        included for interface consistency).

    Returns
    -------
    st_csf : ndarray
        CSF overlap matrix.

    st_sd : ndarray
        Slater-determinant overlap matrix.

    Raises
    ------
    ValueError
        If the dimensions of `st_mo` are inconsistent with the declared
        orbital range.

    Notes
    -----
    - The reference determinant is assumed to be a closed-shell singlet.
    - Only singly excited determinants relative to the reference are
      constructed.
    - Determinants are mapped to spin-adapted CSFs via a transformation
      matrix `T`, such that::

          S_CSF = Tᵀ · S_SD · T

    - The ordering of determinants is preserved throughout.
    """

    # ==================================================================
    # Basic consistency checks
    # ==================================================================
    if nelec % 2 != 0:
        raise ValueError("Closed-shell singlet requires an even number of electrons")

    if S != 0 or Ms != 0:
        raise ValueError("This routine is restricted to singlet states (S=0, Ms=0)")

    if lowest_orbital > highest_orbital:
        raise ValueError("lowest_orbital must be <= highest_orbital")

    if homo_indx < lowest_orbital or homo_indx > highest_orbital:
        raise ValueError("homo_indx must lie within [lowest_orbital, highest_orbital]")

    # ------------------------------------------------------------------
    # Validate MO overlap matrix size
    # ------------------------------------------------------------------
    norb = highest_orbital + 1 - lowest_orbital
    if st_mo.shape[0] != 2 * norb or st_mo.shape[1] != 2 * norb:
        raise ValueError(
            f"MO overlap matrix shape {st_mo.shape} is inconsistent with "
            f"orbital range [{lowest_orbital}, {highest_orbital}]"
        )

    # ------------------------------------------------------------------
    # Define orbital spaces (1-based spatial indices)
    # ------------------------------------------------------------------
    orbital_space = list(range(lowest_orbital, highest_orbital + 1))

    if _active_space is None:
        active_space = orbital_space
    else:
        active_space = list(_active_space)
        if not set(active_space).issubset(orbital_space):
            raise ValueError("Active space orbitals must be within orbital_space")

    # ==================================================================
    # Build reference determinant
    # ==================================================================
    gs = sd.make_ref_det(nelec, homo_indx)

    if len(gs) != nelec:
        raise ValueError("Reference determinant does not contain nelec electrons")

    if len(set(gs)) != len(gs):
        raise ValueError("Reference determinant contains duplicate spin orbitals")

    if not all(lowest_orbital <= abs(o) <= highest_orbital for o in gs):
        raise ValueError("Reference determinant orbitals outside MO overlap range")

    # ==================================================================
    # Build singly excited determinants
    # ==================================================================
    configs0_raw = [tuple(gs)]

    for det in common_sd_basis:
        if len(det) != 2:
            raise ValueError(
                "Each element of common_sd_basis must encode a single excitation (occ, vir)"
            )

        occ, vir = det

        if occ not in gs:
            raise ValueError(f"Occupied orbital {occ} not present in reference determinant")

        if vir in gs:
            raise ValueError(f"Virtual orbital {vir} already occupied in reference determinant")

        ex = sd.make_excitation(gs, occ, vir)
        configs0_raw.append(tuple(ex))

    # ==================================================================
    # Spin adaptation: SD → CSF
    # ==================================================================
    mapped_basis, T = configs_and_T_matrix_singlet(
        configs0_raw,
        active_space,
        orbital_space,
        nelec,
        S,
        Ms,
    )

    if T.shape[0] != 2*(len(configs0_raw) - 1) + 1 :
        raise ValueError("Transformation matrix T has inconsistent dimensions")

    dets = list(mapped_basis)

    # ==================================================================
    # Overlap matrices
    # ==================================================================
    st_mo_dense = st_mo.todense()

    st_sd = sd.slater_overlap_matrix(
        dets, dets, st_mo_dense, complex_valued=False
    )

    if st_sd.shape[0] != st_sd.shape[1]:
        raise ValueError("SD overlap matrix is not square")

    st_csf = T.T @ st_sd @ T

    return st_csf, st_sd


