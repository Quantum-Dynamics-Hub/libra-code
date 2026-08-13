"""Spin-adapted configuration state functions (CSFs).

A spin-free single excitation from orbital 1 to orbital 2 corresponds to two
open-shell determinants:

    |-1,  2>  beta electron remains in orbital 1, alpha electron in orbital 2
    | 1, -2>  alpha electron remains in orbital 1, beta electron in orbital 2

The singlet CSF has total spin S = 0 and M_s = 0.  It is the antisymmetric
spin combination.  The triplet M_s = 0 CSF has S = 1 and is the symmetric
combination.  The spatial part carries the opposite exchange symmetry so the
total electronic wavefunction remains antisymmetric.
"""

import numpy as np

from libra_py.citools import csf


def print_csf(label, terms):
    print(label)
    norm = sum(coeff * coeff for _, coeff in terms)
    for det, coeff in terms:
        print(f"  {coeff:+.8f} * {det}")
    print(f"  norm = {norm:.8f}\n")


dets_with_parity = [
    ((-1, 2), 1),
    ((1, -2), 1),
]

csfs = csf.generate_CSFs_grouped(dets_with_parity)

singlet = csfs[(0.0, 0.0)][0]
triplet_ms0 = csfs[(1.0, 0.0)][0]

print_csf("Singlet CSF, S = 0, M_s = 0", singlet)
print_csf("Triplet CSF, S = 1, M_s = 0", triplet_ms0)

overlap = sum(
    coeff_a * dict(triplet_ms0).get(det, 0.0)
    for det, coeff_a in singlet
)
print(f"<singlet | triplet M_s=0> = {overlap:.8f}")
print("Expected orthogonality:", np.isclose(overlap, 0.0))

