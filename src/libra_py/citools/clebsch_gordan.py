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
.. module:: clebsch_gordan
   :platform: Unix, Windows
   :synopsis: this module implements functions for computing Clebsch-Gordan coefficients
.. moduleauthor:: Alexey V. Akimov, ChatGPT

"""

import math
from typing import Union, List, Dict, Tuple

# ---- helpers ----
def fact(n):
    """Safe factorial: returns 0 for negative ints, factorial for >=0."""
    if n < 0:
        return 0
    return math.factorial(n)

def clebsch_gordan(
    j1: Union[int, float],
    j2: Union[int, float],
    J: Union[int, float],
    m1: Union[int, float],
    m2: Union[int, float],
    M: Union[int, float],
    tol: float = 1e-12
) -> float:
    """
    Compute the Clebsch–Gordan coefficient ⟨j₁ m₁, j₂ m₂ | J M⟩ using the Racah formula.

    This function implements a purely numerical (no `sympy`) version of the Racah formula
    for Clebsch–Gordan coefficients, which describe the coupling of two angular momenta
    `j₁` and `j₂` to form a total angular momentum `J`. The result is a real-valued
    coefficient corresponding to the overlap between the uncoupled basis
    |j₁ m₁⟩ |j₂ m₂⟩ and the coupled basis |J M⟩.

    Parameters
    ----------
    j1, j2, J : int or float
        Angular momentum quantum numbers (can be integers or half-integers).
        Must satisfy the triangular condition:
        |j₁ − j₂| ≤ J ≤ j₁ + j₂.
    m1, m2, M : int or float
        Magnetic quantum numbers corresponding to `j₁`, `j₂`, and `J`.
        Must satisfy m₁ + m₂ = M.
    tol : float, optional
        Numerical tolerance for selection rule enforcement. Default is 1e-12.

    Returns
    -------
    float
        The Clebsch–Gordan coefficient ⟨j₁ m₁, j₂ m₂ | J M⟩.
        Returns 0.0 if the input violates selection rules or factorial arguments
        are negative.

    Notes
    -----
    The implementation follows the Racah formula (see, e.g., Varshalovich *et al.*, 
    *Quantum Theory of Angular Momentum*, 1988):

    $$
    \\langle j_1 m_1, j_2 m_2 | J M \\rangle =
    \\delta_{M, m_1 + m_2} (-1)^{j_1 - j_2 + M}
    \\sqrt{2J + 1}
    \\sqrt{\\frac{(J + j_1 - j_2)! (J - j_1 + j_2)! (j_1 + j_2 - J)!}
                 {(J + j_1 + j_2 + 1)!}} \\\\
    \\times
    \\sqrt{(J + M)! (J - M)! (j_1 + m_1)! (j_1 - m_1)! (j_2 + m_2)! (j_2 - m_2)!}
    \\sum_k
    \\frac{(-1)^k}{
      k! (j_1 + j_2 - J - k)! (j_1 - m_1 - k)! (j_2 + m_2 - k)!
      (J - j_2 + m_1 + k)! (J - j_1 - m_2 + k)! }.
    $$

    All angular momentum values are internally converted to integer "twice-values"
    (2j) to allow half-integer arithmetic using standard integer factorials.

    Selection Rules
    ---------------
    - |m₁| ≤ j₁, |m₂| ≤ j₂, |M| ≤ J  
    - m₁ + m₂ = M  
    - |j₁ − j₂| ≤ J ≤ j₁ + j₂  

    Examples
    --------
    >>> clebsch_gordan(1, 1, 2, 1, 1, 2)
    1.0

    >>> clebsch_gordan(1, 1, 1, 1, 0, 1)
    -0.7071067811865476  # = -1/√2

    >>> clebsch_gordan(1, 1, 1, 0, 1, 1)
    0.7071067811865476   # = +1/√2

    >>> clebsch_gordan(1, 1, 0, 1, -1, 0)
    0.5773502691896258   # = 1/√3

    References
    ----------
    - D. A. Varshalovich, A. N. Moskalev, and V. K. Khersonskii,
      *Quantum Theory of Angular Momentum*, World Scientific (1988)
    - Edmonds, A. R., *Angular Momentum in Quantum Mechanics*, Princeton (1996)

    """
    # Selection rules
    if abs(m1) > j1 + tol or abs(m2) > j2 + tol or abs(M) > J + tol:
        return 0.0
    if abs((m1 + m2) - M) > tol:
        return 0.0
    if (J < abs(j1 - j2) - tol) or (J > (j1 + j2) + tol):
        return 0.0

    # convert to integer twice-values
    j1_2 = int(round(2 * j1))
    j2_2 = int(round(2 * j2))
    J_2  = int(round(2 * J))
    m1_2 = int(round(2 * m1))
    m2_2 = int(round(2 * m2))
    M_2  = int(round(2 * M))

    # factorial arguments (all integers, must be >=0)
    A = (J_2 + j1_2 - j2_2) // 2
    B = (J_2 - j1_2 + j2_2) // 2
    C = (j1_2 + j2_2 - J_2) // 2
    D = (j1_2 + j2_2 + J_2) // 2 + 1

    # check nonnegative
    if min(A, B, C) < 0 or D <= 0:
        return 0.0

    # prefactor (integer exponent for sign)
    sign_exp = (j1_2 - j2_2 + M_2) // 2
    pref1 = (-1) ** sign_exp
    # sqrt term 1
    try:
        sqrt1 = math.sqrt((J_2 + 1) * fact(A) * fact(B) * fact(C) / fact(D))
    except ValueError:
        return 0.0

    # sqrt term 2
    arg1 = (J_2 + M_2) // 2
    arg2 = (J_2 - M_2) // 2
    arg3 = (j1_2 + m1_2) // 2
    arg4 = (j1_2 - m1_2) // 2
    arg5 = (j2_2 + m2_2) // 2
    arg6 = (j2_2 - m2_2) // 2
    if min(arg1, arg2, arg3, arg4, arg5, arg6) < 0:
        return 0.0

    sqrt2 = math.sqrt(
        fact(arg1) * fact(arg2) * fact(arg3) * fact(arg4) * fact(arg5) * fact(arg6)
    )

    pref = pref1 * sqrt1 * sqrt2

    # summation limits k
    k_min = max(
        0,
        ((j2_2 - J_2 - m1_2) // 2),
        ((j1_2 - J_2 + m2_2) // 2)
    )
    k_max = min(
        (j1_2 + j2_2 - J_2) // 2,
        (j1_2 - m1_2) // 2,
        (j2_2 + m2_2) // 2
    )

    total = 0.0
    for k in range(k_min, k_max + 1):
        a = (j1_2 + j2_2 - J_2) // 2 - k
        b = (j1_2 - m1_2) // 2 - k
        c = (j2_2 + m2_2) // 2 - k
        d = (J_2 - j2_2 + m1_2) // 2 + k
        e = (J_2 - j1_2 - m2_2) // 2 + k
        if min(a, b, c, d, e, k) < 0:
            continue
        term = ((-1) ** k) / (fact(k) * fact(a) * fact(b) * fact(c) * fact(d) * fact(e))
        total += term

    return pref * total



# ---- fixed recursive coupling ----

def recursive_couple_spins_int(spin_list_2: List[int]) -> Dict[Tuple[int, int], List[Tuple[Tuple[int, ...], float]]]:
    """
    Recursively couple a list of spin-½ particles using Clebsch–Gordan coefficients (integer 2*M convention).

    This function builds the total spin states for a system of N spin-½ particles 
    (each represented by +1 for α-spin and −1 for β-spin, corresponding to 2*M_i values).  
    The coupling is performed recursively: each additional spin-½ is combined with the 
    previously coupled subsystem using the Clebsch–Gordan coefficients to produce 
    total spin multiplets characterized by total spin S and total projection M.

    Results are expressed in the integer convention for twice-values:
    - `S2 = 2 * S`
    - `M2 = 2 * M`

    Each entry of the returned dictionary corresponds to a particular total spin 
    manifold (S2, M2) and contains all spin configurations contributing to it 
    with their corresponding coefficients in the coupled basis.

    Parameters
    ----------
    spin_list_2 : list of int
        List of twice spin projections for each individual spin-½ particle:
        +1 for α-spin (m = +½), −1 for β-spin (m = −½).  
        For example, `[+1, -1, +1]` represents a three-electron spin configuration
        (αβ α).

    Returns
    -------
    dict[tuple[int, int], list[tuple[tuple[int, ...], float]]]
        A dictionary keyed by `(S2, M2)` representing total spin and projection
        (both doubled integers).  
        Each value is a list of tuples `(spin_config, coeff)`, where:
        
        - **spin_config** (`tuple[int, ...]`): individual spin projections (each ±1)
          for the N spins, in the order they were coupled.
        - **coeff** (`float`): Clebsch–Gordan coefficient weight for that configuration
          in the total spin state.

    Examples
    --------
    >>> recursive_couple_spins_int([+1])
    {(1, 1): [((1,), 1.0)]}

    >>> recursive_couple_spins_int([+1, -1])
    {
        (2, 0): [((1, -1), 0.70710678), ((-1, 1), 0.70710678)],
        (0, 0): [((1, -1), 0.70710678), ((-1, 1), -0.70710678)]
    }

    >>> recursive_couple_spins_int([+1, +1, -1])
    # Returns multiple coupled spin states for S=3/2 and S=1/2 sectors.

    Notes
    -----
    - This function uses a **recursive coupling scheme**:
      it first couples the last N−1 spins, then adds the remaining one.
    - The Clebsch–Gordan coefficients are evaluated using an external function
      `clebsch_gordan_fn(j1, j2, J, m1, m2, M)`, which must be defined in scope.
    - Spin values are handled in the integer “2*M” convention for numerical stability.
    - The resulting dictionary can be used to construct spin-adapted configuration
      state functions (CSFs) or project determinants onto total-spin eigenstates.

    Implementation Details
    ----------------------
    - For N = 1, the result is trivially S = ½ (S2 = 1) with M2 = ±1.
    - For N > 1, the function iteratively couples spin-½ (j₂ = ½) with the 
      previously obtained spin states (S_rest₂ / 2), yielding new total spins 
      `S_total₂ = |S_rest₂ − 1|` and `S_rest₂ + 1`.
    - Terms with vanishing Clebsch–Gordan coefficients (|CG| < 1e−14) are discarded.
    """
    N = len(spin_list_2)
    if N == 0:
        return {}  # no spins
    if N == 1:
        s = spin_list_2[0]
        return {(1, s): [((s,), 1.0)]}  # S2=1, M2=s

    first = spin_list_2[0]
    rest = spin_list_2[1:]
    rest_coupled = recursive_couple_spins_int(rest)
    coupled = {}

    for (S_rest2, M_rest2), states in rest_coupled.items():
        M_total2 = M_rest2 + first
        # allowed total S2 values when coupling S_rest (S_rest2/2) with spin-1/2 (1)
        S_total2_options = [abs(S_rest2 - 1), S_rest2 + 1]
        for S_total2 in S_total2_options:
            j1 = S_rest2 / 2.0
            j2 = 0.5
            J = S_total2 / 2.0
            m1 = M_rest2 / 2.0
            m2 = first / 2.0
            M = M_total2 / 2.0

            cg = clebsch_gordan(j1, j2, J, m1, m2, M)
            if abs(cg) < 1e-14:
                continue
            for spin_config, coeff in states:
                new_conf = (first,) + spin_config
                coupled.setdefault((S_total2, M_total2), []).append((new_conf, coeff * cg))
    return coupled


