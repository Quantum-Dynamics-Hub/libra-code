"""Regression tests for nonorthogonal LDR projection and propagation."""

from pathlib import Path
import sys

import numpy as np
import pytest
import torch

SRC = Path(__file__).resolve().parents[4] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from libra_py.dynamics.ldr_torch.compute import ldr_solver


def assert_close(actual, expected, *, atol=1e-10, rtol=1e-10):
    np.testing.assert_allclose(actual.detach().cpu().numpy(),
                               expected.detach().cpu().numpy(), atol=atol, rtol=rtol)

def make_solver(**overrides):
    params = dict(
        device="cpu", nstates=1, qgrid=[[0.0], [0.8], [1.7]],
        q0=[0.0], p0=[0.0], alpha=[1.0], mass=[1.0], k=[4.0],
        dt=0.17, E=[[0.1, 0.3, 0.8]],
    )
    params.update(overrides)
    solver = ldr_solver(params)
    solver.buildSH()
    return solver


def electronic_basis():
    """Different orthonormal electronic frames at three nuclear centers."""
    angles = torch.tensor([0.0, 0.4, 0.9], dtype=torch.float64)
    c, s = angles.cos(), angles.sin()
    frames = torch.stack([c, -s, s, c], dim=-1).reshape(3, 2, 2)
    # Columns are electronic eigenvectors; rows of vectors are (state, grid).
    vectors = frames.permute(2, 0, 1).reshape(6, 2).to(torch.cdouble)
    return vectors.conj() @ vectors.T


def test_initial_gaussian_equal_to_basis_function_recovers_unit_coefficient():
    solver = make_solver()
    solver.initialize_C()
    assert_close(
        solver.C0, torch.tensor([1, 0, 0], dtype=torch.cdouble), atol=1e-13, rtol=0
    )
    # The default electronic frame must retain inter-center nuclear overlaps.
    assert_close(solver.S, solver.s_nucl.to(torch.cdouble))


def test_initialization_uses_all_electronic_overlaps_and_resets_coefficients():
    se = electronic_basis()
    solver = make_solver(nstates=2, E=torch.zeros(2, 3),
                         s_elec=se, elec_ampl=se[:, 0].reshape(2, 3))
    assert solver.elec_ampl[1].abs().max() > 0.5
    solver.initialize_C()
    expected = torch.zeros(6, dtype=torch.cdouble)
    expected[0] = 1
    assert_close(solver.C0, expected, atol=1e-13, rtol=0)
    solver.elec_ampl = se[:, 3].reshape(2, 3)
    solver.initialize_C()
    expected[0], expected[3] = 0, 1
    assert_close(solver.C0, expected, atol=1e-13, rtol=0)


def test_gaussian_projection_matches_independent_two_dimensional_quadrature():
    solver = make_solver(
        qgrid=[[-0.7, 0.1], [0.4, -0.6], [1.3, 0.8]],
        q0=[0.2, -0.3], p0=[1.1, -0.7], alpha=[1.3, 0.8],
        mass=[1.4, 2.0], k=[0.9, 0.6],
        elec_ampl=[1.0, 0.8+0.2j, 0.7-0.1j],
    )
    x = torch.linspace(-12, 12, 30001, dtype=torch.float64)
    alpha0 = torch.sqrt(solver.k * solver.mass) / 2
    overlaps = torch.ones(3, dtype=torch.cdouble)
    for f in range(2):
        chi = (2 * solver.alpha[f] / torch.pi)**0.25 * torch.exp(
            -solver.alpha[f] * (x[:, None] - solver.qgrid[:, f])**2
        )
        target = (2 * alpha0[f] / torch.pi)**0.25 * torch.exp(
            -alpha0[f] * (x - solver.q0[f])**2
            + 1j * solver.p0[f] * (x - solver.q0[f])
        )
        overlaps *= torch.trapezoid(chi * target[:, None], x, dim=0)
    b = overlaps * solver.elec_ampl[0]
    expected = torch.linalg.solve(solver.S, b)
    expected /= torch.vdot(expected, solver.S @ expected).real.sqrt()
    solver.initialize_C()
    assert_close(solver.C0, expected, atol=1e-11, rtol=1e-11)


def test_legacy_overlap_vector_matches_explicit_single_state_block():
    params = dict(nstates=2, istate=1, E=torch.zeros(2, 3), s_elec=electronic_basis())
    legacy = make_solver(**params, elec_ampl=[0.8, 0.9, 1.0])
    full = make_solver(**params, elec_ampl=[[0, 0, 0], [0.8, 0.9, 1.0]])
    legacy.initialize_C()
    full.initialize_C()
    assert_close(legacy.C0, full.C0)


def test_complex_overlap_square_root_recovers_known_positive_hermitian_matrix():
    # Construct S=A^2 with an explicitly positive Hermitian A. This checks
    # reconstruction in the original basis without using eigh as an oracle.
    root = torch.tensor([
        [2.0, 0.3+0.4j, -0.2+0.1j],
        [0.3-0.4j, 3.0, 0.6-0.2j],
        [-0.2-0.1j, 0.6+0.2j, 4.0],
    ], dtype=torch.cdouble)
    solver = make_solver()
    solver.S = root @ root
    solver.compute_propagator()
    assert_close(solver.S_half, root, atol=1e-12, rtol=1e-12)


@pytest.mark.parametrize("alpha", [[1.0], [[0.7], [1.3], [2.1]]])
def test_complex_gauge_covariance_projection_and_propagation(alpha):
    se = electronic_basis()
    energies = [[0.1, 0.3, 0.8], [0.5, 0.9, 1.4]]
    ref = make_solver(nstates=2, E=energies, s_elec=se, alpha=alpha,
                      elec_ampl=se[:, 0].reshape(2, 3), k=[0.6], p0=[0.7])
    phase = torch.exp(1j * torch.tensor([0.1, 0.7, -1.2, 2.1, -0.4, 1.5],
                                      dtype=torch.float64))
    phased = make_solver(
        nstates=2, E=energies, alpha=alpha, s_elec=phase.conj()[:, None] * se * phase[None, :],
        elec_ampl=(phase.conj() * se[:, 0]).reshape(2, 3), k=[0.6], p0=[0.7],
    )
    for solver in [ref, phased]:
        solver.compute_propagator()
        solver.initialize_C()
        assert_close(solver.S_half @ solver.S_half, solver.S,
                                   atol=1e-12, rtol=1e-12)
        assert_close(solver.U.conj().T @ solver.S @ solver.U, solver.S,
                                   atol=1e-12, rtol=1e-12)
        # Independent reference: exponential of the generalized generator.
        exact_u = torch.linalg.matrix_exp(-1j * solver.dt *
                                          torch.linalg.solve(solver.S, solver.H))
        assert_close(solver.U, exact_u, atol=1e-12, rtol=1e-12)
    assert_close(phased.C0, phase.conj() * ref.C0,
                               atol=1e-12, rtol=1e-12)
    assert_close(phased.U, phase.conj()[:, None] * ref.U * phase[None, :],
                               atol=1e-12, rtol=1e-12)
    ref.C_curr, phased.C_curr = ref.C0.clone(), phased.C0.clone()
    initial_energy = ref.compute_total_energy()
    for _ in range(30):
        ref.C_curr = ref.U @ ref.C_curr
        phased.C_curr = phased.U @ phased.C_curr
    assert_close(phased.C_curr, phase.conj() * ref.C_curr,
                               atol=1e-11, rtol=1e-11)
    assert_close(ref.compute_total_energy(), initial_energy)
    assert_close(phased.compute_total_energy(), initial_energy)
    assert_close(phased.compute_denmat().diagonal(),
                               ref.compute_denmat().diagonal())


@pytest.mark.parametrize("kind", ["list", "numpy", "float32", "float64", "complex64"])
def test_input_dtypes_are_normalized_and_real_overlaps_propagate(kind):
    se = np.ones((3, 3), dtype=np.float64)
    e = np.array([[0.1, 0.3, 0.8]], dtype=np.float64)
    ampl = np.ones(3)
    if kind == "list":
        se, e, ampl = se.tolist(), e.tolist(), ampl.tolist()
    elif kind != "numpy":
        dtype = getattr(torch, kind)
        se, e, ampl = (torch.as_tensor(v, dtype=dtype) for v in (se, e, ampl))
    solver = make_solver(s_elec=se, E=e, elec_ampl=ampl)
    assert solver.s_elec.dtype == solver.elec_ampl.dtype == torch.cdouble
    assert solver.E.dtype == torch.float64
    solver.compute_propagator()
    solver.initialize_C()
    c = solver.U @ solver.C0
    assert_close(torch.vdot(c, solver.S @ c).real,
                               torch.tensor(1.0, dtype=torch.float64))


@pytest.mark.skipif(not torch.cuda.is_available(), reason="CUDA unavailable")
def test_cpu_input_tensors_are_moved_to_cuda():
    solver = make_solver(device="cuda", s_elec=torch.ones(3, 3),
                         E=torch.zeros(1, 3), elec_ampl=torch.ones(3))
    solver.compute_propagator()
    solver.initialize_C()
    for value in [solver.E, solver.s_elec, solver.elec_ampl, solver.S_half, solver.C0]:
        assert value.device.type == "cuda"
    assert torch.isfinite(solver.U @ solver.C0).all()


@pytest.mark.parametrize("params, message", [
    ({"elec_ampl": [[1, 1]]}, "elec_ampl"),
    ({"elec_ampl": [0, 0, 0]}, "nonzero projected norm"),
    ({"s_elec": [[1]]}, "s_elec"),
    ({"E": [[1j, 0, 0]]}, "real electronic energies"),
    ({"E": [[1]]}, "E must have shape"),
    ({"istate": 1}, "istate"),
    ({"alpha": [-1.0]}, "widths must be positive"),
])
def test_invalid_initialization_inputs_are_rejected(params, message):
    with pytest.raises(ValueError, match=message):
        make_solver(**params).initialize_C()


def test_singular_overlap_is_rejected_before_building_propagator():
    solver = make_solver(qgrid=[[0.0], [0.0]], E=[[0.0, 0.0]])
    with pytest.raises(ValueError, match="positive definite"):
        solver.compute_propagator()
