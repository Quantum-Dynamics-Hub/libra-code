"""Quantum-momentum estimates retain the independent QTAG reference semantics."""
import numpy as np
import pytest

from libra_py.dyn.algorithms.quantum_potential import quantum_momentum
from libra_py.dyn.backends.arrays import to_numpy
from qtag_reference.compute import qtag_momentum


@pytest.mark.parametrize("backend", ["numpy", "torch"])
@pytest.mark.parametrize("options", [
    {"mom_calc_type": 0},
    {"mom_calc_type": 2},
    {"mom_calc_type": 1, "momentum_fit": "direct"},
    {"mom_calc_type": 1, "momentum_fit": "legacy", "linfit_beta": 1e-3},
])
def test_quantum_momentum_matches_reference(backend, options):
    """Check raw fields, analytic slopes, and both fits on unequal-width packets."""
    rng = np.random.default_rng(41)
    q, p, s = [rng.normal(size=(2, 5)) for _ in range(3)]
    a = rng.uniform(.7, 2.1, size=q.shape)
    c = rng.normal(size=5) + 1j*rng.normal(size=5)
    args = [q, p, a, s, c]
    if backend == "torch":
        torch = pytest.importorskip("torch")
        args = [torch.as_tensor(value) for value in args]
    actual = quantum_momentum(*args, options)
    expected = qtag_momentum(*args, options)
    for value, reference in zip(actual, expected):
        assert value.shape == q.shape
        np.testing.assert_allclose(to_numpy(value), to_numpy(reference), atol=1e-12)
        if backend == "torch":
            assert value.dtype == args[0].dtype
            assert value.device == args[0].device
