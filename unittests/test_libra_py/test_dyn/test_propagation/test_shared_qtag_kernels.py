"""Shared propagation kernels used by classical and coherent Gaussian dynamics."""
import numpy as np
import pytest
from scipy.linalg import expm

from libra_py.dyn.backends.arrays import to_numpy
from libra_py.dyn.propagation import hermitian_propagator, drift_coordinates


@pytest.mark.parametrize("backend", ["numpy", "torch"])
@pytest.mark.parametrize("columns", [False, True])
def test_hermitian_evolution_matches_matrix_exponential(backend, columns):
    """Real Hamiltonians must also evolve complex vectors on either backend."""
    h = np.array([[.2, .13], [.13, .7]])
    c = np.array([.6+.2j, -.3+.7j])
    if columns:
        c = np.column_stack([c, c.conj()])
    expected = expm(-.4j*h) @ c
    if backend == "torch":
        torch = pytest.importorskip("torch")
        h, c = torch.as_tensor(h), torch.as_tensor(c)
        c.requires_grad_()
    actual = hermitian_propagator(c, h, .4)
    np.testing.assert_allclose(to_numpy(actual), expected, atol=1e-14)
    if backend == "torch":
        actual.real.sum().backward()
        assert c.grad is not None
        assert actual.dtype == c.dtype
        assert actual.device == c.device


@pytest.mark.parametrize("backend", ["numpy", "torch"])
def test_drift_accepts_both_storage_and_packet_layouts(backend):
    """The pure drift preserves inputs and broadcasts mass along the DOF axis."""
    q = np.array([[.1, .2], [.3, .4], [.5, .6]])
    p = np.array([[1., 2.], [3., 4.], [5., 6.]])
    mass = np.array([.2, .3])
    expected = q + .4*mass*p
    if backend == "torch":
        torch = pytest.importorskip("torch")
        q, p, mass = [torch.as_tensor(x) for x in (q, p, mass)]
    before = to_numpy(q).copy()
    np.testing.assert_allclose(to_numpy(drift_coordinates(q, p, mass, .4)), expected)
    np.testing.assert_allclose(
        to_numpy(drift_coordinates(q.T, p.T, mass[:, None], .4)), expected.T)
    np.testing.assert_array_equal(to_numpy(q), before)
