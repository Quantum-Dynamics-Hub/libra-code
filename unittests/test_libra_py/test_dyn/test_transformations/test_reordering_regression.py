"""NumPy regressions adapted from the legacy ``test_reordering.py`` suite."""

from itertools import permutations

import numpy as np
import pytest

from libra_py.dyn.transformations.state_tracking import (
    get_reordering,
    hungarian_algorithm,
    munkres_kuhn,
    permutation_matrix,
)


def _matrix(rows):
    """Construct a real overlap matrix from a readable nested sequence."""

    return np.asarray(rows, dtype=complex)


# These identity, multiply mixed, and ambiguous overlap matrices reproduce the
# cases in unittests/test_reordering.py. The second value is the optimal
# assignment expected by that suite; the third is the consistently oriented
# greedy result. They differ for ambiguous mixed-state matrices.
MIXED_OVERLAP_CASES = [
    pytest.param(np.eye(4, dtype=complex), [0, 1, 2, 3], [0, 1, 2, 3], id="identity-4"),
    pytest.param(
        _matrix([
            [1, 0, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1],
            [0, 1, 0, 0],
        ]),
        [0, 2, 3, 1],
        [0, 2, 3, 1],
        id="cycle-3-of-4",
    ),
    pytest.param(
        _matrix([
            [0, 1, 0, 0],
            [0, 0, 1, 0],
            [1, 0, 0, 0],
            [0, 0, 0, 1],
        ]),
        [1, 2, 0, 3],
        [1, 2, 0, 3],
        id="cycle-first-3",
    ),
    pytest.param(
        _matrix([
            [1, 0, 0, 0],
            [0, 0.76, 0.79, 0],
            [0, 0.79, 0.80, 0],
            [0, 0, 0, 1],
        ]),
        [0, 2, 1, 3],
        [0, 1, 2, 3],
        id="double-mixing",
    ),
    pytest.param(
        _matrix([
            [1, 0, 0, 0, 0],
            [0, 0.71, 0.75, 0.77, 0],
            [0, 0.72, 0.75, 0.78, 0],
            [0, 0.73, 0.76, 0.80, 0],
            [0, 0, 0, 0, 1],
        ]),
        [0, 2, 1, 3, 4],
        [0, 3, 1, 2, 4],
        id="triple-mixing-a",
    ),
    pytest.param(
        _matrix([
            [1, 0, 0, 0, 0, 0, 0],
            [0, 0.73, 0.75, 0, 0, 0, 0],
            [0, 0.76, 0.76, 0, 0, 0, 0],
            [0, 0, 0, 0.71, 0.72, 0, 0],
            [0, 0, 0, 0.74, 0.78, 0, 0],
            [0, 0, 0, 0, 0, 0.77, 0.78],
            [0, 0, 0, 0, 0, 0.79, 0.81],
        ]),
        [0, 2, 1, 3, 4, 5, 6],
        [0, 2, 1, 4, 3, 6, 5],
        id="three-double-mixing-blocks",
    ),
    pytest.param(
        _matrix([
            [1, 0, 0, 0, 0],
            [0, 0.77, 0.78, 0.78, 0],
            [0, 0.76, 0.75, 0.89, 0],
            [0, 0.75, 0.73, 0.77, 0],
            [0, 0, 0, 0, 1],
        ]),
        [0, 2, 3, 1, 4],
        [0, 2, 3, 1, 4],
        id="triple-mixing-b",
    ),
    pytest.param(
        _matrix([
            [1, 0, 0, 0, 0],
            [0, 0.65, 0.78, 0.77, 0],
            [0, 0.89, 0.75, 0.76, 0],
            [0, 0.77, 0.54, 0.65, 0],
            [0, 0, 0, 0, 1],
        ]),
        [0, 2, 1, 3, 4],
        [0, 3, 1, 2, 4],
        id="triple-mixing-c",
    ),
]


@pytest.mark.parametrize("overlap, expected, greedy_expected", MIXED_OVERLAP_CASES)
@pytest.mark.parametrize("algorithm", [munkres_kuhn, hungarian_algorithm])
def test_legacy_mixed_overlap_optimal_reorderings(
    overlap, expected, greedy_expected, algorithm
):
    """Optimal algorithms reproduce the assignments in the legacy suite."""

    del greedy_expected
    actual = algorithm(overlap, np.zeros_like(overlap), 0.0)
    np.testing.assert_array_equal(actual, expected)


@pytest.mark.parametrize("overlap, expected, greedy_expected", MIXED_OVERLAP_CASES)
def test_mixed_overlap_greedy_reordering_uses_consistent_orientation(
    overlap, expected, greedy_expected
):
    """Greedy outputs retain old-to-new orientation in ambiguous cases."""

    del expected
    np.testing.assert_array_equal(get_reordering(overlap), greedy_expected)


@pytest.mark.parametrize("nstates", [2, 3, 4])
def test_every_exact_permutation_for_all_deterministic_algorithms(nstates):
    """Cover every exact permutation through four states.

    Inputs are normalized to the module convention
    ``S[old_state, new_state]`` before being passed to each algorithm.
    """

    energies = np.zeros((nstates, nstates))
    for expected_tuple in permutations(range(nstates)):
        expected = np.asarray(expected_tuple)
        old_to_new = permutation_matrix(expected)
        munkres_overlap = old_to_new.T
        greedy_overlap = old_to_new.T
        hungarian_overlap = old_to_new.T

        np.testing.assert_array_equal(
            munkres_kuhn(munkres_overlap, energies, 0.0, 0), expected
        )
        np.testing.assert_array_equal(get_reordering(greedy_overlap), expected)
        np.testing.assert_array_equal(
            hungarian_algorithm(hungarian_overlap, energies, 0.0), expected
        )
