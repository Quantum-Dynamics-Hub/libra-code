"""Construct and apply phase-corrected and LD projectors for three states."""

import numpy as np

from libra_py.dyn.control_params import DynControlParams
from libra_py.dyn.transformations.projectors import (
    compute_projector,
    update_projection,
)


permutation = np.array([1, 2, 0])
phases = np.array([1.0j, -1.0, np.exp(0.25j)])
time_overlap = np.zeros((3, 3), dtype=complex)
time_overlap[np.arange(3), permutation] = phases
energies = np.diag([-0.20, 0.00, 0.30])

# Algorithm 21 uses the Hungarian assignment. Phase correction makes the
# diagonal of S @ P positive and real after the new states have been reordered.
params = DynControlParams(
    state_tracking_algo=21,
    do_phase_correction=1,
    phase_correction_tol=1.0e-10,
)
projector = compute_projector(params, energies, time_overlap)
print("permutation-and-phase projector:\n", projector)
print("corrected overlap S @ P:\n", time_overlap @ projector)

# Algorithm -1 is local diabatization. update_projection imports
# orthogonalized_T from local_diabatization.py and returns a unitary transform.
nonunitary_overlap = np.array(
    [
        [0.91, 0.15, -0.02],
        [-0.11, 0.86, 0.18],
        [0.04, -0.14, 0.88],
    ],
    dtype=complex,
)
ld_projector = update_projection(
    {"state_tracking_algo": -1},
    nonunitary_overlap,
    energies_current=energies,
)
print("LD projector:\n", ld_projector)
print("LD unitarity error:", np.max(np.abs(ld_projector.conj().T @ ld_projector - np.eye(3))))
