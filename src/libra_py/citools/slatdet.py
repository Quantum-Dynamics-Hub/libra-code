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

