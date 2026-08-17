"""CI-state time-overlaps from spin-adapted determinant overlaps.

Package interfaces such as CP2K, DFTB+, and MOPAC ultimately provide two
kinds of data to `citools`:

1. one-electron molecular-orbital time-overlaps, S_mo(t, t+dt)
2. CI expansions of excited states in terms of spin-free single excitations

This example builds a small model by hand.  The active orbital space contains
two occupied orbitals and one virtual orbital.  The two spin-free excitations
are:

    [2, 3]  HOMO -> LUMO
    [1, 3]  lower occupied orbital -> LUMO

The ground state is implicit and is represented by the closed-shell reference.
"""

import numpy as np

from libra_py.citools import ci


# A one-electron time-overlap matrix between spatial MOs at two consecutive
# nuclear geometries.  It is not exactly symmetric because it connects two
# different times: <phi_i(t) | phi_j(t+dt)>.
spatial_time_overlap = np.array([
    [0.99, 0.04, 0.02],
    [0.03, 0.98, 0.05],
    [0.01, 0.06, 0.97],
])

# Restricted alpha/beta spin functions are orthogonal, so the doubled
# spin-orbital matrix has zero off-diagonal spin blocks.
spin_time_overlap = np.kron(np.eye(2), spatial_time_overlap)

# data[1] is a list over excited states; each excited state is a list of
# spin-free single excitations [occupied, virtual].
#
# data[2] has the corresponding CI amplitudes.  These are deliberately simple
# normalized mixtures in the two-configuration space.
data_t = [
    [0.10, 0.18],
    [
        [[2, 3], [1, 3]],
        [[2, 3], [1, 3]],
    ],
    [
        [0.96, 0.28],
        [-0.28, 0.96],
    ],
]
data_t_dt = [
    [0.11, 0.19],
    [
        [[2, 3], [1, 3]],
        [[2, 3], [1, 3]],
    ],
    [
        [0.95, 0.31],
        [-0.31, 0.95],
    ],
]

params = {
    "homo_indx": 2,
    "nocc": 1,
    "nvirt": 1,
    "nelec": 4,
    "nstates": 3,
    "active_space": [1, 2, 3],
}

st_ci = ci.overlap(spin_time_overlap, data_t, data_t_dt, params)

print("Spatial MO time-overlap:")
print(spatial_time_overlap)
print("\nCI state time-overlap, including ground state as state 0:")
print(st_ci)
print("\nDiagonal elements near 1 indicate state continuity over the time step.")
print("Off-diagonal elements are finite-time approximations to nonadiabatic mixing.")

