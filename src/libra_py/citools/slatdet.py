from __future__ import annotations
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
.. module:: slatdet
   :platform: Unix, Windows
   :synopsis: this module implements auxiliary functions for working with Slater determiants
.. moduleauthor:: Alexey V. Akimov, ChatGPT

"""

import itertools
import math
import numpy as np
from typing import Generator, Tuple, List, Sequence, Any

# ---------- helpers ----------

def permutation_parity(original: Sequence[Any], sorted_target: Sequence[Any]) -> int:
    """
    Compute the permutation parity (+1 or -1) that maps one sequence to another.

    Given two sequences containing the same unique elements, this function determines 
    whether the permutation required to reorder `original` into `sorted_target` is 
    even or odd. It does so by constructing the index mapping between the two sequences 
    and counting the number of inversions.

    Parameters
    ----------
    original : sequence of Any
        The starting ordering of elements. Typically represents a specific configuration
        of spin-orbitals or basis functions.
    sorted_target : sequence of Any
        The target ordering containing the same elements as `original`, but arranged
        in the desired canonical order (e.g., from `canonical_sort_key`).

    Returns
    -------
    int
        The parity of the permutation:
        
        - `+1` if the permutation is even (can be achieved with an even number of swaps)  
        - `-1` if the permutation is odd (requires an odd number of swaps)

    Examples
    --------
    >>> permutation_parity((2, 1), (1, 2))
    -1

    >>> permutation_parity((3, 1, 2), (1, 2, 3))
    +1

    >>> permutation_parity((1, 2, 3), (1, 2, 3))
    +1

    Notes
    -----
    - Both `original` and `sorted_target` must contain the same elements, without duplicates.
    - The algorithm counts inversions in the index mapping derived from the correspondence 
      between the two sequences.
    - This function is typically used in the construction of antisymmetrized states such as 
      Slater determinants and configuration state functions (CSFs), where permutation parity 
      determines the sign of a term under electron exchange.
    """
    # Build list of target indices for each element in original,
    # taking care in case values repeat (not expected here).
    target_positions = []
    used = {}
    for idx, val in enumerate(sorted_target):
        used.setdefault(val, []).append(idx)
    # For each original element, take the first unused index from used[val]
    taken = {val: 0 for val in used}
    for val in original:
        poslist = used[val]
        take_idx = taken[val]
        target_positions.append(poslist[take_idx])
        taken[val] += 1

    # Count inversions in target_positions
    inv = 0
    for i in range(len(target_positions)):
        for j in range(i + 1, len(target_positions)):
            if target_positions[i] > target_positions[j]:
                inv += 1

    return -1 if (inv % 2 == 1) else 1




def canonical_sort_key(x: int) -> tuple[int, int]:
    """
    Canonical sorting key for spin-orbitals.

    This function defines a consistent ordering of spin-orbitals used to 
    canonicalize Slater determinants. It ensures that spin-orbitals are 
    first ordered by the absolute value of their orbital index (spatial orbital),
    and then by spin: α-spin (positive index) precedes β-spin (negative index).

    Parameters
    ----------
    x : int
        Spin-orbital label, where positive values (e.g., +1, +2, +3) represent 
        α-spin orbitals and negative values (e.g., -1, -2, -3) represent β-spin orbitals.

    Returns
    -------
    tuple of (int, int)
        Sorting key `(abs(x), spin_order)` where:
        
        - `abs(x)` ensures grouping by spatial orbital index.  
        - `spin_order` is 0 for α-spin (x > 0) and 1 for β-spin (x < 0), 
          so that α-spin comes before β-spin in canonical ordering.

    Examples
    --------
    >>> sorted([+2, -1, +1, -2], key=canonical_sort_key)
    [1, -1, 2, -2]

    >>> canonical_sort_key(-3)
    (3, 1)

    Notes
    -----
    This key is typically used together with functions such as 
    `generate_determinants_with_parity` to define a canonical 
    determinant ordering and compute permutation parities.
    """
    return (abs(x), 0 if x > 0 else 1)



# ---------- determinant generator (returns sorted det + parity) ----------
def generate_determinants_with_parity(
    active_orbitals: List[int],
    N: int,
    allow_double_occupancy: bool = True
) -> Generator[Tuple[Tuple[int, ...], int], None, None]:
    """
    Generate all possible Slater determinants (configurations) from a given set of 
    active orbitals, including their permutation parities.

    Each spatial orbital in `active_orbitals` contributes two spin-orbitals: 
    α (represented by +i) and β (represented by -i). The function yields all 
    unique combinations of `N` spin-orbitals consistent with the Pauli principle 
    and, optionally, with or without double occupancy of the same spatial orbital.

    The permutation parity is computed relative to the canonically sorted order 
    of spin-orbitals, allowing subsequent antisymmetrization when constructing 
    configuration state functions (CSFs).

    Parameters
    ----------
    active_orbitals : list of int
        List of spatial orbital indices that define the active space, e.g. [1, 2, 3].
    N : int
        Number of electrons to distribute among the spin-orbitals.
    allow_double_occupancy : bool, optional
        If True (default), both α and β spins of the same orbital may be occupied.
        If False, configurations containing both +i and -i for the same orbital are excluded.

    Yields
    ------
    tuple of (tuple of int, int)
        A pair `(determinant, parity)` where:
        
        - **determinant** (`tuple[int, ...]`):  
          Spin-orbital occupations represented as integers, 
          with positive values for α-spin and negative for β-spin,
          sorted canonically (e.g., `(-2, -1, +1, +3)`).
        
        - **parity** (`int`):  
          The parity (+1 or -1) of the permutation needed to transform
          the original ordering into canonical order.

    Examples
    --------
    >>> list(generate_determinants_with_parity([1, 2], 2))
    [((-2, -1), 1), ((-2, 1), -1), ((-2, 2), 1), ((-1, 1), 1), ((-1, 2), -1), ((1, 2), 1)]

    >>> list(generate_determinants_with_parity([1, 2], 2, allow_double_occupancy=False))
    [((-2, -1), 1), ((-2, 1), -1), ((-1, 2), -1), ((1, 2), 1)]

    Notes
    -----
    This function assumes the following helper functions are defined elsewhere:
    
    - `canonical_sort_key(spin_orbital: int) -> Any`  
      Defines the canonical ordering of spin-orbitals.

    - `permutation_parity(original: tuple, sorted_tuple: tuple) -> int`  
      Returns +1 for even permutations and -1 for odd ones.
    """
    spin_orbitals: List[int] = []
    for i in active_orbitals:
        spin_orbitals.append(+i)   # α
        spin_orbitals.append(-i)   # β

    for combo in itertools.combinations(spin_orbitals, N):
        if not allow_double_occupancy:
            # forbid both spins of same orbital
            if any((+i in combo) and (-i in combo) for i in active_orbitals):
                continue
        original: Tuple[int, ...] = tuple(combo)
        sorted_combo: Tuple[int, ...] = tuple(sorted(original, key=canonical_sort_key))
        parity: int = permutation_parity(original, sorted_combo)
        yield sorted_combo, parity

def generate_single_excitations(active_orbitals: List[int], nelec: int):
    """
    Generate the ground-state determinant first, followed by all single-excitation
    determinants (sorted + parity), for a closed-shell reference inside an active space.

    Parameters
    ----------
    active_orbitals : list[int]
        Ordered list of spatial orbital numbers in the active space.
        Example: [1,2,3,4,5,6,7,8,9,10]

    nelec : int
        Number of electrons in the active space. Must be even for closed-shell.
        Example: 10 -> 5 occupied spatial orbitals (paired).

    Yields
    ------
    (det_sorted, parity)
        det_sorted: tuple of signed spin-orbitals (canonical order)
        parity: +1 or -1
    """
    if nelec % 2 != 0:
        raise ValueError("Closed-shell reference requires an even number of electrons.")

    n_occ = nelec // 2
    occ = active_orbitals[:n_occ]
    virt = active_orbitals[n_occ:]

    # -----------------------------
    # 1. Ground state determinant
    # -----------------------------
    ref_raw = []
    for o in occ:
        ref_raw.append(+o)  # alpha
        ref_raw.append(-o)  # beta

    ref_sorted = tuple(sorted(ref_raw, key=canonical_sort_key))
    ref_parity = permutation_parity(tuple(ref_raw), ref_sorted)

    # yield ground-state first
    yield ref_sorted, ref_parity

    # -----------------------------
    # 2. Single excitations
    # -----------------------------
    for o in occ:
        for v in virt:

            # α-spin excitation: +o → +v
            det_raw = list(ref_raw)
            det_raw.remove(+o)
            det_raw.append(+v)
            det_sorted = tuple(sorted(det_raw, key=canonical_sort_key))
            parity = permutation_parity(tuple(det_raw), det_sorted)
            yield det_sorted, parity

            # β-spin excitation: -o → -v
            det_raw = list(ref_raw)
            det_raw.remove(-o)
            det_raw.append(-v)
            det_sorted = tuple(sorted(det_raw, key=canonical_sort_key))
            parity = permutation_parity(tuple(det_raw), det_sorted)
            yield det_sorted, parity

def slater_overlap_matrix(dets_A, dets_B, S_orb, complex_valued=False):
    """
    Compute the matrix of overlaps between two possibly distinct sets of
    Slater determinants (α/β spins orthogonal).

    Parameters
    ----------
    dets_A : list[tuple[int]]
        Bra determinants, tuples of signed orbital indices.
        Positive = alpha, negative = beta.

    dets_B : list[tuple[int]]
        Ket determinants, tuples of signed orbital indices.

    S_orb : np.ndarray
        Spatial orbital overlap matrix (n_orb, n_orb).

    complex_valued : bool, optional
        If True, results are complex-valued. Default is False.

    Returns
    -------
    S_AB : np.ndarray
        Overlap matrix of shape (len(dets_A), len(dets_B)),
        with elements <D_I^(A) | D_J^(B)>.
    """
    dtype = complex if complex_valued else float
    nA, nB = len(dets_A), len(dets_B)
    S_AB = np.zeros((nA, nB), dtype=dtype)

    # Precompute α and β orbital indices (0-based)
    alpha_A = [np.array([abs(o) - 1 for o in d if o > 0], dtype=int) for d in dets_A]
    beta_A  = [np.array([abs(o) - 1 for o in d if o < 0], dtype=int) for d in dets_A]
    alpha_B = [np.array([abs(o) - 1 for o in d if o > 0], dtype=int) for d in dets_B]
    beta_B  = [np.array([abs(o) - 1 for o in d if o < 0], dtype=int) for d in dets_B]

    for i in range(nA):
        a_i, b_i = alpha_A[i], beta_A[i]
        for j in range(nB):
            a_j, b_j = alpha_B[j], beta_B[j]

            # Build α and β submatrices
            S_a = S_orb[np.ix_(a_i, a_j)]
            S_b = S_orb[np.ix_(b_i, b_j)]

            # Product of determinants
            S_AB[i, j] = np.linalg.det(S_a) * np.linalg.det(S_b)

    return S_AB


def make_ref_det(nelec, homo_indx):
    """
    Construct a closed-shell reference determinant in spin-orbital notation.

    Spin orbitals are labeled by integers:
      +i  → spin-up orbital i
      -i  → spin-down orbital i

    Parameters
    ----------
    nelec : int
        Total number of electrons (must be even).
    homo_indx : int
        HOMO spatial orbital index (1-based).

    Returns
    -------
    list[int]
        Spin-orbital occupation list, using +i (alpha) and -i (beta).
    """
    if nelec % 2 != 0:
        raise ValueError("Closed-shell reference requires an even number of electrons")

    nocc = nelec // 2
    lowest_occ = homo_indx - nocc + 1

    if lowest_occ < 1:
        raise ValueError("Not enough orbitals below HOMO to place all electrons")

    return [
        spin
        for i in range(lowest_occ, homo_indx + 1)
        for spin in (i, -i)
    ]



def make_excitation(ref_det, occ, vir):
    """
    Generate a single excitation from a reference Slater determinant.

    This function replaces one occupied spin orbital (`occ`) in the
    reference determinant with a virtual spin orbital (`vir`).

    Parameters
    ----------
    ref_det : list of int
        Reference Slater determinant represented as a list of occupied
        spin orbitals.

    occ : int
        Occupied spin orbital to be removed.

    vir : int
        Virtual spin orbital to be inserted.

    Returns
    -------
    list of int
        New Slater determinant corresponding to the excitation.

    Raises
    ------
    ValueError
        If `occ` is not present in `ref_det`.

    Notes
    -----
    - The returned determinant is a new list; the reference determinant
      is not modified.
    - No checks are performed for duplicate occupations or Pauli
      violations.
    - Orbital ordering is preserved except for the replaced index.
    """

    if vir in ref_det:
        raise ValueError(F"Orbital {vir} is already present in reference determinant, is not valid virtual orbital")

    if occ not in ref_det:
        raise ValueError(F"Orbital {occ} is not present in reference determinant")


    res = list(ref_det)

    # O(N) lookup; acceptable for CI-size determinants
    idx = res.index(occ)
    res[idx] = vir

    return res

