# *********************************************************************************
# * Copyright (C) 2025 Alexey V. Akimov
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
.. moduleauthor:: Alexey V. Akimov, ChatGPT

"""

from collections import Counter
from typing import List, Tuple, Iterable, Dict, Any
import numpy as np

from . import csf
from . import slatdet as sd

from liblibra_core import *

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


from typing import List, Tuple, Any

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
    >>> import numpy as np
    >>> # Example input from MOPAC/Libra
    >>> configs0_raw = [(6, -6, 7, -7, 9, -8)]
    >>> active_space = [6, 7, 8, 9, 10, 11]
    >>> nelec = 6
    >>> max_unpaired = 0   # singlet configurations only
    >>>
    >>> min_basis, (all_confs, all_phases) = build_minimal_csf_basis(
    ...     configs0_raw, active_space, nelec, max_unpaired
    ... )
    >>>
    >>> print(len(min_basis))
    19
    >>> print(all_confs[0])
    (6, -6, 7, -7, 8, -8)
    >>> print(all_phases[0])
    1
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
) -> "CMATRIX":
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
    T : CMATRIX
        Complex-valued configuration–CSF transformation matrix of shape
        `(n_determinants, n_CSFS)`.

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

    return T


def configs_and_T_matrix(
    configs0_raw: List[Tuple[int, ...]],
    active_space: List[int],
    nelec: int,
    S: int,
    Ms: int
) -> Tuple[List[Tuple[int, ...]], "CMATRIX"]:
    """
    Generate the minimal active-space configurations and the configuration-to-CSF
    transformation matrix for a given CAS (Complete Active Space) and spin.

    Parameters
    ----------
    configs0_raw : list[tuple[int]]
        List of raw configurations from Libra/MOPAC (signed orbital indices).
    active_space : list[int]
        List of orbitals defining the active space.
    nelec : int
        Number of active electrons.
    S : int
        Total spin quantum number.
    Ms : int
        Spin projection quantum number (Ms).

    Returns
    -------
    mapped_basis : list[tuple[int, ...]]
        List of minimal configurations mapped to active-space indices (starting from 1),
        with signs preserved (positive = alpha, negative = beta).
    T : CMATRIX
        Complex-valued configuration-to-CSF transformation matrix.

    Example
    -------
    >>> mapped_basis, T = configs_and_T_matrix(
    ...     configs0_raw,
    ...     active_space=[6, 7, 8, 9, 10, 11],
    ...     nelec=6,
    ...     S=0,
    ...     Ms=0
    ... )
    >>> print("Number of configurations:", len(mapped_basis))
    >>> print("First configuration:", mapped_basis[0])
    >>> print("Shape of T:", T.num_of_rows, "x", T.num_of_cols)
    """

    # 1. Create minimal SD basis with spin constraint 2*Ms
    min_basis, (all_confs, all_phases) = build_minimal_csf_basis(
        configs0_raw, active_space, nelec, 2*Ms
    )

    # 2. Map minimal basis to active-space indices
    mapped_basis = map_to_active_indices(all_confs, active_space)

    # 3. Compute configuration-to-CSF transformation matrix
    T = conf2csf_matrix(min_basis, all_confs, all_phases, S, Ms)

    return mapped_basis, T


