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
.. module:: csf
   :platform: Unix, Windows
   :synopsis: this module implements functions for constructing configuration state functions (CSFs)
.. moduleauthor:: Alexey V. Akimov, ChatGPT

"""

from collections import defaultdict
import numpy as np
from typing import Iterable, Tuple, Dict, List

from . import clebsch_gordan


def group_key_from_det(det):
    """
    Compute a grouping key for a determinant based on orbital occupation.

    The function counts the number of electrons occupying each spatial orbital 
    in the determinant and separates orbitals into:
      - doubly occupied orbitals
      - singly occupied orbitals

    Parameters
    ----------
    det : tuple[int]
        Determinant represented as a tuple of spin-orbitals.
        Positive values indicate α-spin, negative values β-spin, and 
        the absolute value corresponds to the spatial orbital index.
        Example: `(1, -2, 2, -1)`.

    Returns
    -------
    tuple[tuple[int, ...], tuple[int, ...]]
        A tuple `(double_orbs, single_orbs)` where:
        - `double_orbs` : tuple of orbital indices that are doubly occupied (2 electrons)
        - `single_orbs` : tuple of orbital indices that are singly occupied (1 electron)
        
        The orbital indices in each tuple are sorted in ascending order.

    Example
    -------
    >>> group_key_from_det((1, -1, 2, -3))
    ((1,), (2, 3))

    >>> group_key_from_det((1, 2, -2, 3, -3))
    ((2, 3), (1,))
    
    Notes
    -----
    - This key is typically used to group determinants that share the same spatial 
      occupation pattern when building Configuration State Functions (CSFs).
    """
    occ = {}
    for x in det:
        i = abs(x)
        occ[i] = occ.get(i, 0) + 1
    double_orbs = tuple(sorted([i for i, c in occ.items() if c == 2]))
    single_orbs = tuple(sorted([i for i, c in occ.items() if c == 1]))
    return (double_orbs, single_orbs)



def build_unpaired_spinlist(det, single_orbs):
    """
    Construct a spin list for the unpaired electrons in a determinant.

    Given a determinant and a list of singly-occupied orbitals (`single_orbs`),
    this function returns a tuple of spin projections for the unpaired electrons
    in the specified orbital order. The spins are represented as integers:
        +1 → α-spin (spin up)
        -1 → β-spin (spin down)

    Parameters
    ----------
    det : tuple[int]
        Determinant represented as a tuple of spin-orbitals. 
        Positive values indicate α-spin, negative values β-spin.
        Example: `(1, -2, 3, -4)`.
    single_orbs : iterable[int]
        Ordered list of spatial orbital indices that are singly occupied in this determinant.
        The order determines the order of spins in the returned tuple.

    Returns
    -------
    tuple[int, ...]
        Tuple of +1/-1 values corresponding to the spin projections of the unpaired 
        electrons in the order specified by `single_orbs`.

    Raises
    ------
    RuntimeError
        If a singly-occupied orbital listed in `single_orbs` is missing from the determinant.

    Examples
    --------
    >>> build_unpaired_spinlist((1, -2, 3), [2, 3])
    (-1, +1)
    >>> build_unpaired_spinlist((1, 2, -3), [2, 3])
    (+1, -1)
    """
    spins = []
    for orb in single_orbs:
        if orb in det:
            spins.append(+1)
        elif -orb in det:
            spins.append(-1)
        else:
            raise RuntimeError("Unpaired orbital missing from determinant")
    return tuple(spins)



def generate_CSFs_grouped(
    dets_with_parity: Iterable[Tuple[Tuple[int, ...], int]],
    tol: float = 1e-12
) -> Dict[Tuple[float, float], List[List[Tuple[Tuple[int, ...], float]]]]:
    """
    Build Configuration State Functions (CSFs) grouped by spatial-occupation patterns.

    This function constructs spin-adapted CSFs from a list of determinants 
    with their permutation parities. Determinants are first grouped by their 
    spatial orbital occupation patterns (e.g., which orbitals are doubly or singly 
    occupied), and then spin couplings are applied to unpaired electrons within 
    each group using Clebsch–Gordan coefficients.

    Parameters
    ----------
    dets_with_parity : iterable of (tuple[int, ...], int)
        Iterable of determinants and their canonicalization parities.
        - det: tuple of integers representing spin-orbitals, e.g., `(1, -2, 3, -4)`
        - parity: +1 or -1, the permutation parity used when canonicalizing the determinant
    tol : float, optional
        Tolerance for comparing M values when selecting coupled spin states.
        Default is 1e-12.

    Returns
    -------
    csf_dict : dict[tuple[float, float], list[list[tuple[tuple[int, ...], float]]]]
        Dictionary keyed by `(S, M)` total spin and projection:
        - `S` (float): total spin quantum number
        - `M` (float): total spin projection
        Each value is a list of CSFs for that `(S, M)` sector.  
        Each CSF is a list of `(det, coeff)` pairs, where:
        - `det`: determinant tuple from `dets_with_parity`
        - `coeff`: float coefficient of that determinant in the CSF

    Notes
    -----
    - Determinants are first grouped by spatial occupation pattern (doubly and singly 
      occupied orbitals). Each group is processed independently.
    - Closed-shell groups (no unpaired electrons) produce pure singlet CSFs `(S=0, M=0)`.
    - Open-shell groups use `recursive_couple_spins_int` to generate all allowed 
      `(S2, M2)` couplings (where S2=2*S, M2=2*M) for the unpaired electrons.
    - CSFs are normalized and phase-consistent: the element of largest absolute value 
      is forced positive.
    - Vectors with all zero coefficients are skipped.
    - The output dictionary converts the internal `defaultdict` to a regular `dict`.

    Examples
    --------
    >>> dets_with_parity = [((1, -2, 3), 1), ((1, -3, 2), -1)]
    >>> csf_dict = generate_CSFs_grouped(dets_with_parity)
    >>> for key, csfs in csf_dict.items():
    ...     print(key, csfs)
    (0.5, 0.5) [[((1, -2, 3), 0.70710678), ((1, -3, 2), 0.70710678)]]
    (1.5, 1.5) [[((1, -2, 3), 1.0)]]

    Dependencies
    ------------
    - `group_key_from_det(det)`: function that computes a hashable key from the determinant's
      spatial occupation pattern.
    - `build_unpaired_spinlist(det, single_orbs)`: function that extracts +1/-1 spin 
      assignments for unpaired electrons.
    - `clebsch_gordan.recursive_couple_spins_int(spins_2)`: function that recursively couples 
      spins to produce all `(S2, M2)` multiplets.
    - NumPy (`np`) for array handling and vector normalization.

    """
    # 1) group determinants by pattern (double_orbs, single_orbs)
    groups = defaultdict(list)
    for det, parity in dets_with_parity:
        key = group_key_from_det(det)
        groups[key].append((tuple(det), float(parity)))

    # temporary storage: for each (S,M,group_key) keep list of (det, coeff)
    temp = defaultdict(list)

    # 2) process each group independently
    for group_key, entries in groups.items():
        double_orbs, single_orbs = group_key
        nd = len(entries)
        det_list = [e[0] for e in entries]
        parity_list = [e[1] for e in entries]

        # handle closed-shell group (no unpaired electrons)
        if len(single_orbs) == 0:
            # Each determinant is a pure singlet CSF (S=0, M=0)
            for det, parity in entries:
                temp[(0.0, 0.0, group_key)].append([(det, parity * 1.0)])  # wrap here
            continue

        # For this open-shell pattern we will build CSF vectors for each (S,M)
        SM_candidates = set()
        per_det_couplings = []
        for det, parity in entries:
            spins_2 = build_unpaired_spinlist(det, single_orbs)
            Ms2_det = sum(spins_2)
            coupled = clebsch_gordan.recursive_couple_spins_int(spins_2)
            mapping = {}
            for (S2, M2), state_list in coupled.items():
                if abs(M2 - Ms2_det) > tol:
                    continue
                coeff_for_config = 0.0
                for spin_config, coeff in state_list:
                    if tuple(spin_config) == tuple(spins_2):
                        coeff_for_config += coeff
                if abs(coeff_for_config) > 0.0:
                    mapping[(S2, M2)] = coeff_for_config
                    SM_candidates.add((S2, M2))
            per_det_couplings.append((det, parity, mapping))

        for (S2, M2) in sorted(SM_candidates):
            vec = np.zeros(nd, dtype=float)
            for i, (det, parity, mapping) in enumerate(per_det_couplings):
                coeff = mapping.get((S2, M2), 0.0)
                if coeff != 0.0:
                    vec[i] = float(coeff * parity)
            if np.linalg.norm(vec) < 1e-14:
                continue
            vec /= np.linalg.norm(vec)
            imax = int(np.argmax(np.abs(vec)))
            if vec[imax] < 0:
                vec *= -1.0
            csf = []
            for i, val in enumerate(vec):
                if abs(val) > 1e-14:
                    csf.append((det_list[i], float(val)))
            S = S2 / 2.0
            M = M2 / 2.0
            temp[(S, M, group_key)].append(csf)

    # 3) flatten per-group CSFs into final csf_dict keyed by (S, M)
    csf_dict = defaultdict(list)
    for (S, M, group_key), csf_list in temp.items():
        for csf in csf_list:
            if isinstance(csf, tuple) and isinstance(csf[0], (tuple, list)) and isinstance(csf[1], (int, float)):
                csf_dict[(S, M)].append([csf])
            elif isinstance(csf, list):
                csf_dict[(S, M)].append(csf)
            else:
                raise ValueError(f"Unexpected CSF entry format: {csf}")

    return dict(csf_dict)



def print_csfs(csfs):
    """
    Nicely print Configuration State Functions (CSFs) returned by `generate_CSFs_grouped`.

    This function iterates over all `(S, M)` sectors in the `csfs` dictionary
    and prints each CSF with its determinant(s) and corresponding coefficients.
    It handles both standard CSF lists `[(det, coeff), ...]` and accidental 
    tuple entries `(det, coeff)` safely.

    Parameters
    ----------
    csfs : dict[tuple[float, float], list[list[tuple[tuple[int, ...], float]]]]
        Dictionary of CSFs as returned by `generate_CSFs_grouped`.
        Keys are `(S, M)` total spin and projection.
        Values are lists of CSFs, each CSF is a list of `(det, coeff)` pairs.

    Prints
    ------
    For each `(S, M)` sector:
      - Line indicating the spin: "S,Ms = (S,M)"
      - Each CSF with index: "  CSF #i:"
      - Each determinant and its coefficient: "    det  coeff" with coeff formatted as `+/-0.######`

    Example
    -------
    >>> csfs = {
    ...     (0.5, 0.5): [[((1, -2, 3), 0.707107), ((1, -3, 2), 0.707107)]],
    ...     (1.5, 1.5): [[((1, -2, 3), 1.0)]]
    ... }
    >>> print_csfs(csfs)
    S,Ms = (0.5,0.5)
      CSF #1:
        (1, -2, 3)  +0.707107
        (1, -3, 2)  +0.707107
    S,Ms = (1.5,1.5)
      CSF #1:
        (1, -2, 3)  +1.000000
    """
    for key in sorted(csfs.keys()):
        S, M = key
        print(f"S,Ms = ({S},{M})")
        for i, csf in enumerate(csfs[key], start=1):
            print(f"  CSF #{i}:")
            # handle both [(det,coeff)] and accidental tuples
            if isinstance(csf, tuple) and len(csf) == 2 and isinstance(csf[0], (tuple, list)):
                csf = [csf]
            for det, coeff in csf:
                print(f"    {det}  {coeff:+.6f}")



