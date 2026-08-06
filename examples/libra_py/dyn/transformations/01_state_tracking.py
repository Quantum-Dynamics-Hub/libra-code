"""Track the identities of three states across one electronic-structure step."""

import numpy as np

from libra_py.dyn.transformations.state_tracking import (
    get_reordering,
    get_stochastic_reordering2,
    hungarian_algorithm,
    make_cost_matrix,
    permutation_matrix,
    permute_states,
)


# Rows are old states and columns are newly labelled states. The largest
# overlaps describe the cyclic identity map old 0->new 1, 1->2, and 2->0.
time_overlap = np.array(
    [
        [0.05, 0.98, 0.03],
        [0.02, 0.04, 0.97],
        [0.99, 0.01, 0.02],
    ],
    dtype=complex,
)
energies = np.diag([-0.20, 0.00, 0.30])

# alpha=0 gives pure overlap tracking. A positive alpha and selectors 1--3
# suppress assignments between energetically distant states.
cost = make_cost_matrix(
    time_overlap, energies, alpha=0.5, scaling_function=2
)
greedy = get_reordering(time_overlap)
optimal = hungarian_algorithm(
    time_overlap, energies, alpha=0.5, scaling_function=2
)
stochastic = get_stochastic_reordering2(
    time_overlap, rng=np.random.default_rng(4)
)

print("energy-aware overlap score:\n", cost)
print("greedy permutation:", greedy)
print("optimal permutation:", optimal)
print("sampled permutation:", stochastic)
print("optimal permutation matrix:\n", permutation_matrix(optimal))

# Each row is a trajectory-specific permutation. This maps active-state labels
# from the old basis into the newly labelled basis.
active_states = np.array([0, 2])
trajectory_permutations = np.asarray([optimal, [2, 0, 1]])
print("old active states:", active_states)
print("new active states:", permute_states(trajectory_permutations, active_states))
