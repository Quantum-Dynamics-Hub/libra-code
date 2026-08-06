"""End-to-end projector conventions used by state tracking and propagation."""

import numpy as np
import pytest

from libra_py.dyn.control_params import DynControlParams
from libra_py.dyn.propagation.electronic import propagate_electronic_method
from libra_py.dyn.transformations.projectors import compute_projector


PERMUTATION = np.array([1, 2, 0])


def _phase_permuted_overlap():
    overlap = np.zeros((3, 3), dtype=complex)
    overlap[np.arange(3), PERMUTATION] = [1.0j, -1.0, np.exp(0.3j)]
    return overlap


def test_projector_columns_define_state_reassignment_and_phase_correction():
    """P[perm[i],i] maps old label i to new raw label perm[i]."""

    overlap = _phase_permuted_overlap()
    projector = compute_projector(
        DynControlParams(state_tracking_algo=21, do_phase_correction=1),
        np.diag([-0.2, 0.0, 0.3]),
        overlap,
    )

    np.testing.assert_array_equal(np.argmax(np.abs(projector), axis=0), PERMUTATION)
    np.testing.assert_allclose(overlap @ projector, np.eye(3), atol=1.0e-14)


def test_projector_transforms_new_raw_hamiltonian_to_old_state_labels():
    """P† H_new P orders diagonal energies by the tracked old labels."""

    projector = compute_projector(
        {"state_tracking_algo": 21},
        np.diag([10.0, 20.0, 30.0]),
        np.abs(_phase_permuted_overlap()),
    )
    hamiltonian_new = np.diag([10.0, 20.0, 30.0])
    dynamically_consistent = projector.conj().T @ hamiltonian_new @ projector

    np.testing.assert_allclose(
        np.diag(dynamically_consistent),
        np.diag(hamiltonian_new)[PERMUTATION],
    )


@pytest.mark.parametrize("method", [0, 1, 2, 3, 4, 5, 10, 11, 12, 13, 14, 15])
def test_cpp_projector_methods_map_coefficients_with_T_at_zero_time(method):
    """Dynamics.cpp options 0--5 and 10--15 end with C_new = T C_old at dt=0."""

    projector = compute_projector(
        {"state_tracking_algo": 21},
        np.diag([0.0, 0.1, 0.2]),
        np.abs(_phase_permuted_overlap()),
    )[None, ...]
    coefficients = np.array([[1.0, 2.0j, -0.5]])
    zero = np.zeros((1, 3, 3), dtype=complex)

    result = propagate_electronic_method(
        coefficients,
        ham_current=zero,
        ham_previous=zero,
        hvib_current=zero,
        hvib_previous=zero,
        dt=0.0,
        method=method,
        projector=projector,
    )

    expected = np.einsum("tij,tj->ti", projector, coefficients)
    np.testing.assert_allclose(result, expected)


def test_cpp_special_methods_6_7_8_preserve_their_distinct_basis_behavior():
    projector = compute_projector(
        {"state_tracking_algo": 21},
        np.diag([0.0, 0.1, 0.2]),
        np.abs(_phase_permuted_overlap()),
    )[None, ...]
    coefficients = np.array([[1.0, 2.0j, -0.5]])
    zero = np.zeros((1, 3, 3), dtype=complex)
    common = dict(
        ham_current=zero,
        ham_previous=zero,
        hvib_current=zero,
        hvib_previous=zero,
        dt=0.0,
        projector=projector,
    )

    np.testing.assert_allclose(
        propagate_electronic_method(coefficients, method=6, **common), coefficients
    )
    np.testing.assert_allclose(
        propagate_electronic_method(coefficients, method=7, **common),
        np.einsum("tji,tj->ti", projector.conj(), coefficients),
    )
    np.testing.assert_allclose(
        propagate_electronic_method(coefficients, method=8, **common), coefficients
    )


@pytest.mark.parametrize(
    "method, expected_energy",
    [
        (0, 2.0), (1, 2.0), (2, 2.0),
        (3, 2.0), (4, 3.0), (5, 3.0),
        (6, 1.0), (7, 3.0), (8, 3.0),
        (10, 2.0), (11, 2.0), (12, 2.0),
        (13, 2.0), (14, 3.0), (15, 3.0),
    ],
)
def test_electronic_integrator_selector_uses_cpp_time_point_formula(
    method, expected_energy
):
    """Check accumulated phases for scalar old/current Hamiltonians."""

    coefficients = np.ones((1, 1), dtype=complex)
    old_h = np.array([[[1.0]]])
    new_h = np.array([[[3.0]]])
    old_v = np.array([[[2.0]]])
    new_v = np.array([[[4.0]]])
    result = propagate_electronic_method(
        coefficients,
        ham_current=new_h,
        ham_previous=old_h,
        hvib_current=new_v,
        hvib_previous=old_v,
        dt=0.2,
        method=method,
    )
    np.testing.assert_allclose(result[0, 0], np.exp(-0.2j * expected_energy))
