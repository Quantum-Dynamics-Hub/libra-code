import numpy as np

from libra_py.dyn.models import PhenolModel
from libra_py.units import Angst

from _test_utils import assert_engine_builds_adiabatic


def test_phenol_shapes_symmetry_and_overlap():
    model = PhenolModel()
    q = np.asarray(
        [
            [0.95 * Angst, 1.05 * Angst],
            [0.20, 0.65],
        ]
    )

    result = model.evaluate(q)

    assert result["H_dia"].shape == (2, 3, 3)
    assert result["dH_dia"].shape == (2, 2, 3, 3)
    assert np.allclose(result["H_dia"], np.swapaxes(result["H_dia"], -1, -2).conj())
    assert np.allclose(result["dH_dia"], np.swapaxes(result["dH_dia"], -1, -2).conj())
    assert np.allclose(result["S_dia"], np.broadcast_to(np.eye(3), (2, 3, 3)))


def test_phenol_derivatives_match_finite_difference():
    model = PhenolModel()
    points = [
        np.asarray([0.98 * Angst, 0.25]),
        np.asarray([1.10 * Angst, 0.80]),
    ]

    for point in points:
        result = model.evaluate(point.reshape(2, 1))
        expected = _finite_difference_derivative(model, point)
        assert np.allclose(result["dH_dia"][0], expected, rtol=1.0e-4, atol=1.0e-6)


def test_phenol_works_through_hamiltonian_engine():
    q = np.asarray(
        [
            [0.95 * Angst, 1.05 * Angst, 1.15 * Angst],
            [0.20, 0.65, 1.10],
        ]
    )
    assert_engine_builds_adiabatic(PhenolModel(), q)


def _finite_difference_derivative(model, point, step=1.0e-5):
    derivatives = np.zeros((2, 3, 3), dtype=complex)
    for dof in range(2):
        dq = np.zeros_like(point)
        dq[dof] = step
        plus = model.diabatic((point + dq).reshape(2, 1))[0]
        minus = model.diabatic((point - dq).reshape(2, 1))[0]
        derivatives[dof] = (plus - minus) / (2.0 * step)
    return derivatives
