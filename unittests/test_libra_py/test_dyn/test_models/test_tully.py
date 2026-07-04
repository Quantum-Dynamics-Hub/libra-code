import numpy as np

from libra_py.dyn.backends import backend
from libra_py.dyn.core.storage import TensorStorage
from libra_py.dyn.core.trajectory import Trajectory
from libra_py.dyn.hamiltonians import (
    HamiltonianEngine,
    compute_adiabatic_from_diabatic,
    compute_diabatic,
)
from libra_py.dyn.models import TullyModel1, TullyModel2, TullyModel3


def test_tully_diabatic_hamiltonians_match_reference_formulas():
    x = np.array([-2.0, 0.0, 1.5])
    cases = [
        (TullyModel1(), _reference_model1),
        (TullyModel2(), _reference_model2),
        (TullyModel3(), _reference_model3),
    ]

    for model, reference_fn in cases:
        result = model.evaluate(np.asarray([x]))
        expected_h, expected_dh = reference_fn(x)

        assert np.allclose(result["H_dia"], expected_h)
        assert np.allclose(result["dH_dia"], expected_dh)
        assert np.allclose(result["DC1_dia"], np.zeros_like(expected_dh))
        assert set(result) == {"H_dia", "dH_dia", "DC1_dia", "S_dia"}


def test_tully_adiabatic_energies_match_closed_form_two_state_values():
    x = np.array([-2.0, -0.25, 0.25, 2.0])

    for model in (TullyModel1(), TullyModel2(), TullyModel3()):
        storage, traj = _compute_from_hamiltonian_machinery(model, x)
        expected = _two_state_eigenvalues(storage.ham_dia[0, traj.tbf_ids])

        assert np.allclose(
            np.diagonal(storage.ham_adi[0, traj.tbf_ids], axis1=-2, axis2=-1).real,
            expected,
        )
        assert np.allclose(
            storage.ham_adi[0, traj.tbf_ids],
            np.swapaxes(storage.ham_adi[0, traj.tbf_ids].conj(), -1, -2),
        )


def test_tully_adiabatic_derivative_couplings_match_finite_difference():
    x = np.array([-2.0, -0.5, 0.5, 2.0])

    for model in (TullyModel1(), TullyModel2(), TullyModel3()):
        storage, traj = _compute_from_hamiltonian_machinery(model, x)
        expected = np.array([_finite_difference_dc01(model, point) for point in x])

        dc1_adi = storage.dc1_adi[0, traj.tbf_ids]
        assert np.allclose(dc1_adi[:, 0, 0, 1].real, expected, atol=2.0e-5)
        assert np.allclose(dc1_adi[:, 0, 1, 0], -dc1_adi[:, 0, 0, 1])


def test_tully_model_writes_expected_values_through_hamiltonian_engine():
    model = TullyModel1()
    q = np.array([-2.0, 0.5, 2.0])
    momentum = np.array([10.0, 12.0, 14.0])
    inverse_mass = 1.0 / 2000.0

    storage = TensorStorage(
        backend=backend,
        ntraj=1,
        ndof=1,
        nstates=2,
        ntbf_initial=3,
        ntbf_capacity=3,
    )
    storage.allocate_hamiltonian_derivatives(der_lvl=1)
    storage.q[0, :, 0] = q
    storage.p[0, :, 0] = momentum
    storage.iM[0, :, 0] = inverse_mass

    traj = Trajectory(0)
    traj.tbf_ids = [0, 1, 2]
    HamiltonianEngine(backend).evaluate(traj, storage, model, rep="adiabatic")

    reference_storage, reference_traj = _compute_from_hamiltonian_machinery(model, q)
    reference_idx = reference_traj.tbf_ids
    expected_nac01 = momentum * inverse_mass * reference_storage.dc1_adi[0, reference_idx, 0, 0, 1].real

    assert np.allclose(storage.ham_dia[0, :3], reference_storage.ham_dia[0, reference_idx])
    assert np.allclose(storage.ham_adi[0, :3], reference_storage.ham_adi[0, reference_idx])
    assert np.allclose(storage.dc1_adi[0, :3], reference_storage.dc1_adi[0, reference_idx])
    assert np.allclose(storage.nac_adi[0, :3, 0, 1].real, expected_nac01)
    assert np.allclose(storage.hvib_adi[0, :3], storage.ham_adi[0, :3] - 1j * storage.nac_adi[0, :3])


def _reference_model1(x, A=0.010, B=1.600, C=0.005, D=1.000):
    exp_pos = np.exp(-B * x)
    exp_neg = np.exp(B * x)
    v11 = np.where(x > 0, A * (1.0 - exp_pos), -A * (1.0 - exp_neg))
    dv11 = np.where(x > 0, A * B * exp_pos, A * B * exp_neg)
    z = np.exp(-D * x * x)
    v12 = C * z
    dv12 = -2.0 * x * C * D * z
    return _reference_two_state(x, v11, -v11, v12, dv11, -dv11, dv12)


def _reference_model2(x, A=0.100, B=0.280, C=0.015, D=0.060, E=0.050):
    zh = np.exp(-B * x * x)
    zc = np.exp(-D * x * x)
    v11 = np.zeros_like(x)
    v22 = E - A * zh
    v12 = C * zc
    dv11 = np.zeros_like(x)
    dv22 = 2.0 * A * B * x * zh
    dv12 = -2.0 * C * D * x * zc
    return _reference_two_state(x, v11, v22, v12, dv11, dv22, dv12)


def _reference_model3(x, A=0.0006, B=0.1000, C=0.9000):
    v11 = np.zeros_like(x) + A
    v22 = -v11
    exp_left = np.exp(C * x)
    exp_right = np.exp(-C * x)
    v12 = np.where(x <= 0, B * exp_left, B * (2.0 - exp_right))
    dv12 = np.where(x <= 0, B * C * exp_left, B * C * exp_right)
    return _reference_two_state(x, v11, v22, v12, np.zeros_like(x), np.zeros_like(x), dv12)


def _reference_two_state(x, v11, v22, v12, dv11, dv22, dv12):
    h = np.zeros((len(x), 2, 2), dtype=complex)
    dh = np.zeros((len(x), 1, 2, 2), dtype=complex)
    h[:, 0, 0] = v11
    h[:, 1, 1] = v22
    h[:, 0, 1] = v12
    h[:, 1, 0] = v12
    dh[:, 0, 0, 0] = dv11
    dh[:, 0, 1, 1] = dv22
    dh[:, 0, 0, 1] = dv12
    dh[:, 0, 1, 0] = dv12
    return h, dh


def _two_state_eigenvalues(h):
    a = h[..., 0, 0].real
    b = h[..., 1, 1].real
    c = h[..., 0, 1].real
    center = 0.5 * (a + b)
    radius = np.sqrt((0.5 * (a - b)) ** 2 + c * c)
    return np.stack((center - radius, center + radius), axis=-1)


def _compute_from_hamiltonian_machinery(model, x):
    storage = TensorStorage(
        backend=backend,
        ntraj=1,
        ndof=1,
        nstates=2,
        ntbf_initial=len(x),
        ntbf_capacity=len(x),
    )
    storage.allocate_hamiltonian_derivatives(der_lvl=1)
    traj = Trajectory(0)
    traj.tbf_ids = list(range(len(x)))
    storage.q[0, traj.tbf_ids, 0] = x

    compute_diabatic(storage, traj, model, der_lvl=1)
    compute_adiabatic_from_diabatic(storage, traj, der_lvl=1)
    return storage, traj


def _finite_difference_dc01(model, x, step=1.0e-5):
    center = _basis_transform_at(model, x).real
    plus = _basis_transform_at(model, x + step).real
    minus = _basis_transform_at(model, x - step).real

    plus = _align_eigenvector_phases(center, plus)
    minus = _align_eigenvector_phases(center, minus)
    derivative = (plus - minus) / (2.0 * step)
    return (center.T @ derivative)[0, 1]


def _basis_transform_at(model, x):
    storage, traj = _compute_from_hamiltonian_machinery(model, np.asarray([x]))
    return storage.basis_transform[0, traj.tbf_ids[0]]


def _align_eigenvector_phases(reference, vectors):
    aligned = vectors.copy()
    for state in range(reference.shape[1]):
        if np.dot(reference[:, state], aligned[:, state]) < 0.0:
            aligned[:, state] *= -1.0
    return aligned
