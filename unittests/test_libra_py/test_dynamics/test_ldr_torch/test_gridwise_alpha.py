"""Independent Gaussian integral checks for gridwise LDR exponents."""
from pathlib import Path
import sys

import numpy as np
import pytest
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "src"))
from libra_py.dynamics.ldr_torch.compute import ldr_solver


def make_solver(alpha, ndim=2, **overrides):
    params = dict(device="cpu", nstates=1,
                  qgrid=np.array([[-0.7, 0.1], [0.4, -0.6], [1.3, 0.8]])[:, :ndim],
                  alpha=alpha, q0=[0.2, -0.3][:ndim], p0=[1.1, -0.7][:ndim],
                  mass=[1.4, 2.0][:ndim], k=[0.9, 0.6][:ndim],
                  E=[[0.1, 0.3, 0.8]], dt=0.17)
    params.update(overrides)
    solver = ldr_solver(params)
    solver.buildSH()
    return solver


@pytest.mark.parametrize("ndim", [1, 2])
def test_gridwise_integrals_and_projection_against_quadrature(ndim):
    alpha = torch.tensor([[0.7, 1.4], [1.6, 0.6], [1.1, 2.0]], dtype=torch.float64)[:, :ndim]
    solver = make_solver(alpha, ndim)
    x = torch.linspace(-12, 12, 30001, dtype=torch.float64)
    overlap, kinetic, position = [], [], []
    b = torch.ones(3, dtype=torch.cdouble)
    a0 = (solver.k * solver.mass).sqrt() / 2
    for f in range(ndim):
        dx = x[:, None] - solver.qgrid[:, f]
        g = (2 * alpha[:, f] / torch.pi)**0.25 * torch.exp(-alpha[:, f] * dx**2)
        grad = -2 * alpha[:, f] * dx * g
        overlap.append(torch.trapezoid(g[:, :, None] * g[:, None, :], x, dim=0))
        kinetic.append(torch.trapezoid(grad[:, :, None] * grad[:, None, :], x, dim=0)
                       / (2 * solver.mass[f]))
        position.append(torch.trapezoid(x[:, None, None] * g[:, :, None] * g[:, None, :], x, dim=0))
        target = (2*a0[f]/torch.pi)**0.25 * torch.exp(
            -a0[f]*(x-solver.q0[f])**2 + 1j*solver.p0[f]*(x-solver.q0[f]))
        b *= torch.trapezoid(g * target[:, None], x, dim=0)
    s = torch.stack(overlap).prod(dim=0)
    t = torch.zeros_like(s)
    c = torch.tensor([0.4+0.2j, -0.1+0.5j, 0.6-0.3j], dtype=torch.cdouble)
    solver.C_curr = c
    expected_position = []
    for f in range(ndim):
        other = torch.ones_like(s)
        for k in range(ndim):
            if k != f:
                other *= overlap[k]
        t += kinetic[f] * other
        expected_position.append((torch.vdot(c, (position[f]*other).to(c.dtype) @ c)
                                  / torch.vdot(c, s.to(c.dtype) @ c)).real)
    torch.testing.assert_close(solver.s_nucl, s, atol=1e-11, rtol=1e-11)
    torch.testing.assert_close(solver.t_nucl, t, atol=1e-11, rtol=1e-11)
    torch.testing.assert_close(torch.stack(solver.compute_average_pos()),
                               torch.stack(expected_position), atol=1e-11, rtol=1e-11)
    expected_c = torch.linalg.solve(s.to(c.dtype), b)
    expected_c /= torch.vdot(expected_c, s.to(c.dtype) @ expected_c).real.sqrt()
    solver.initialize_C()
    torch.testing.assert_close(solver.C0, expected_c, atol=1e-11, rtol=1e-11)


@pytest.mark.parametrize("uniform", [1.5, [1.5], [1.5, 0.75]])
def test_uniform_and_repeated_gridwise_widths_agree(uniform):
    alpha = torch.as_tensor(uniform, dtype=torch.float64).expand(3, 2).clone()
    a, b = make_solver(uniform), make_solver(alpha)
    for solver in (a, b):
        solver.compute_propagator()
        solver.initialize_C()
        solver.C_curr = solver.U @ solver.C0
    for name in ("s_nucl", "t_nucl", "S", "H", "U", "C0", "C_curr"):
        torch.testing.assert_close(getattr(a, name), getattr(b, name), atol=1e-12, rtol=1e-12)
    torch.testing.assert_close(torch.stack(a.compute_average_pos()), torch.stack(b.compute_average_pos()))


def test_projection_recovers_one_gridwise_basis_gaussian():
    alpha = torch.tensor([[0.7, 1.4], [1.6, 0.6], [1.1, 2.0]], dtype=torch.float64)
    mass = torch.tensor([1.4, 2.0], dtype=torch.float64)
    solver = make_solver(alpha, q0=[0.4, -0.6], p0=[0., 0.], k=4*alpha[1]**2/mass)
    solver.initialize_C()
    torch.testing.assert_close(solver.C0, torch.tensor([0., 1., 0.], dtype=torch.cdouble),
                               atol=1e-12, rtol=1e-12)


@pytest.mark.parametrize("alpha", [[], [1., 2., 3.], [[1., 2.]], [[1.], [2.], [3.]],
                                  [[[1., 2.]]], [0., 1.], [float("nan"), 1.],
                                  [[1., 2.], [1., float("inf")], [1., 2.]]])
def test_invalid_gridwise_widths_are_rejected(alpha):
    with pytest.raises(ValueError, match="alpha|widths"):
        make_solver(alpha)


def test_gridwise_widths_saved_with_shape_and_values(tmp_path):
    alpha = [[0.7, 1.4], [1.6, 0.6], [1.1, 2.0]]
    solver = make_solver(alpha, prefix=str(tmp_path / "gridwise"))
    solver.save()
    data = torch.load(tmp_path / "gridwise.pt", weights_only=True)
    torch.testing.assert_close(data["alpha"], torch.tensor(alpha, dtype=torch.float64))
