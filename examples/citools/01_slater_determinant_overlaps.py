"""Slater-determinant overlaps from molecular-orbital time-overlaps.

This example uses two electrons in two spatial orbitals.  The determinant
labels are:

    (1, -1)  closed-shell reference
    (-1, 2)  alpha excitation 1 -> 2, in canonical determinant notation
    (1, -2)  beta excitation 1 -> 2

The doubled spin-orbital overlap matrix has alpha and beta blocks and zero
off-diagonal spin blocks.  This encodes the orthogonality of alpha and beta
spin functions.
"""

import numpy as np

from libra_py.citools import slatdet as sd


def spin_orbital_indices(det, nspatial):
    """Map signed orbital labels to doubled spin-orbital matrix indices."""
    return [abs(o) - 1 if o > 0 else nspatial + abs(o) - 1 for o in det]


def full_spin_orbital_overlap(dets, spin_overlap):
    """Reference implementation using one full determinant per SD pair."""
    nspatial = spin_overlap.shape[0] // 2
    out = np.zeros((len(dets), len(dets)))

    for i, bra in enumerate(dets):
        bra_idx = spin_orbital_indices(bra, nspatial)
        for j, ket in enumerate(dets):
            ket_idx = spin_orbital_indices(ket, nspatial)
            out[i, j] = np.linalg.det(spin_overlap[np.ix_(bra_idx, ket_idx)])

    return out


alpha_overlap = np.array([
    [1.00, 0.20],
    [0.30, 1.10],
])
beta_overlap = np.array([
    [0.95, 0.10],
    [0.15, 1.05],
])
spin_overlap = np.block([
    [alpha_overlap, np.zeros_like(alpha_overlap)],
    [np.zeros_like(beta_overlap), beta_overlap],
])

dets = [(1, -1), (-1, 2), (1, -2)]
ordering_phases = [sd.alpha_beta_ordering_phase(det) for det in dets]

factorized = sd.slater_overlap_matrix(
    dets,
    dets,
    spin_overlap,
    phases_A=ordering_phases,
    phases_B=ordering_phases,
    spin_orbital_matrix=True,
)
full = full_spin_orbital_overlap(dets, spin_overlap)

print("Doubled spin-orbital overlap matrix:")
print(spin_overlap)
print("\nDeterminants:")
for det, phase in zip(dets, ordering_phases):
    print(f"  {str(det):>12}  canonical -> alpha/beta phase = {phase:+d}")

print("\nFactorized alpha/beta SD overlap:")
print(factorized)
print("\nFull spin-orbital determinant reference:")
print(full)
print("\nAgreement:", np.allclose(factorized, full))
