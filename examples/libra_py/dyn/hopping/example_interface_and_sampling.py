"""Common ``tsh_method`` interface, batches, and stochastic proposals."""

from _example_setup import use_repo_sources

use_repo_sources()

import numpy as np

from libra_py.dyn.control_params import DynControlParams
from libra_py.dyn.hopping import (
    TSH_METHODS,
    hop,
    hop_proposal_probabilities,
    propose_hops,
)


def main():
    print("C++-compatible tsh_method values:")
    for number, name in TSH_METHODS.items():
        print(f"  {number:2d}: {name}")

    density = np.diag([0.2, 0.7, 0.1]).astype(complex)
    hvib = np.diag([0.0, 0.01, 0.02]).astype(complex)
    mash = hop_proposal_probabilities(
        DynControlParams(tsh_method=6), density, hvib, active_states=0
    )
    print("\nMASH through common interface:", mash)

    no_hop = hop_proposal_probabilities(
        DynControlParams(tsh_method=-1), density, hvib, active_states=2
    )
    print("No-hop option from state 2:    ", no_hop)

    current = {
        "ham_dia": np.array([[0.1, 0.01], [0.01, -0.1]]),
        "d1ham_dia": np.array([[[1.0, 0.0], [0.0, -1.0]]]),
    }
    previous = {"ham_dia": np.diag([-0.1, 0.1])}
    batch = hop_proposal_probabilities(
        DynControlParams(tsh_method=3, rep_lz=0),
        None,
        None,
        active_states=[0, 0, 0],
        ham=[current, current, current],
        ham_prev=[previous, previous, previous],
        momentum=np.array([[0.5, 1.0, 4.0]]),
        inverse_mass=[1.0],
    )
    print("\nBatched LZ probabilities at velocities 0.5, 1, and 4:")
    print(batch)

    print("\nDeterministic sampling with specified random numbers:")
    probabilities = [0.2, 0.3, 0.5]
    for random_number in (0.05, 0.25, 0.75):
        print(f"  ksi={random_number:.2f} -> state {hop(0, probabilities, random_number)}")

    rng = np.random.default_rng(2)
    proposed = propose_hops(
        [[0.2, 0.8], [0.9, 0.1], [0.4, 0.6]], [0, 1, 0], rng
    )
    print("Seeded batch proposal:        ", proposed)


if __name__ == "__main__":
    main()
