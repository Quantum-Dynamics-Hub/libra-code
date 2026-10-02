"""Check the numerical conventions used by the manuscript comparison examples."""

import importlib.util
from pathlib import Path

import numpy as np
import pytest
import torch

ROOT = Path(__file__).resolve().parents[4]
spec = importlib.util.spec_from_file_location(
    "ldr_manuscript_examples", ROOT / "examples/libra_py/dynamics/ldr_torch/_common.py"
)
examples = importlib.util.module_from_spec(spec)
spec.loader.exec_module(examples)


def settings_1d():
    return dict(title="Test SAC", model="tully1", ldr_bounds=[[-1.0, 1.0]],
                dvr_bounds=[[-3.0, 3.0]], ldr_spacing=[0.5], dvr_spacing=[0.1],
                q0=[0.0], p0=[0.0], sigma=[1 / np.sqrt(8)], mass=[2000.0],
                istate=0, dt=1.0, duration=3.0, save_every=2, snapshots=[0.0, 3.0])


@pytest.mark.parametrize("ndim", [1, 2])
def test_reconstructed_ldr_density_matches_exact_basis_gaussian(ndim):
    config = settings_1d()
    if ndim == 2:
        config.update(model="flv", gamma=0.08, istate=1,
                      ldr_bounds=[[1.0, 3.0], [-1.0, 1.0]],
                      dvr_bounds=[[0.0, 4.0], [-2.0, 2.0]],
                      ldr_spacing=[0.5, 0.5], dvr_spacing=[0.1, 0.1],
                      q0=[2.0, 0.0], p0=[0.0, 0.0],
                      sigma=[1 / np.sqrt(8)] * 2, mass=[20000.0, 6667.0])
    solver, axes, frames = examples.make_ldr(config, "cpu")
    evaluation_axes = examples.grid_axes(config, "dvr", "cpu")
    density = examples.ldr_density(solver, axes, frames, evaluation_axes)
    expected = examples.gaussian(examples.coordinates(evaluation_axes), config).abs().square()
    np.testing.assert_allclose(density.numpy(), expected.numpy(), atol=1e-12, rtol=1e-10)


@pytest.mark.parametrize("count", [25, 26])
def test_dvr_plane_wave_dispersion_for_odd_and_even_grids(count):
    config = settings_1d()
    dx, mode, mass, dt = 0.25, -4, 2.0, 0.3
    config.update(dvr_bounds=[[0.0, (count - 1) * dx]], dvr_spacing=[dx],
                  q0=[3.0], mass=[mass], dt=dt)
    solver = examples.DVRReference(config, "cpu")
    solver.V.zero_()
    solver.expV_half = torch.eye(2, dtype=torch.cdouble).expand(count, 2, 2)
    momentum = 2 * np.pi * mode / (count * dx)
    psi = torch.exp(1j * momentum * solver.Q[0]) / np.sqrt(count * dx)
    solver.psi_r_dia.zero_()
    solver.psi_r_dia[:, 0] = psi
    solver.step()
    energy = momentum**2 / (2 * mass)
    expected = psi * np.exp(-1j * energy * dt)
    np.testing.assert_allclose(solver.psi_r_dia[:, 0].numpy(), expected.numpy(), atol=1e-12)
    _, kinetic, _, _, _, norm = examples.record_dvr(solver)
    assert float(kinetic) == pytest.approx(energy, abs=1e-12)
    assert float(norm) == pytest.approx(1.0, abs=1e-12)


def test_comparison_outputs_include_final_time_and_matching_initial_state(tmp_path):
    config = settings_1d()
    outputs = [examples.simulate(config, method, tmp_path, "cpu") for method in ["ldr", "dvr"]]
    for data in outputs:
        assert data["time"].tolist() == [0.0, 2.0, 3.0]
        assert data["density_time"].tolist() == [0.0, 3.0]
        assert data["rho"].shape == (3, 2, 2)
        assert data["norm"].dtype == torch.float64
        np.testing.assert_allclose(data["norm"].numpy(), 1.0, atol=1e-12)
        assert (tmp_path / (data["method"] + ".pt")).exists()
        assert (tmp_path / (data["method"] + "_settings.json")).exists()
    np.testing.assert_allclose(outputs[0]["density"][0], outputs[1]["density"][0], atol=1e-12)


def test_1d_phase_alignment_makes_adjacent_overlaps_positive():
    config = settings_1d()
    q = torch.linspace(-5, 5, 51, dtype=torch.float64)[None, :]
    _, frames = torch.linalg.eigh(examples.potential(q, config))
    phase = torch.exp(1j * torch.arange(102, dtype=torch.float64).reshape(51, 1, 2))
    aligned = examples.align_1d(frames * phase)
    overlap = (aligned[:-1].conj() * aligned[1:]).sum(dim=-2)
    assert torch.all(overlap.real > 0)
    np.testing.assert_allclose(overlap.imag.numpy(), 0, atol=1e-14)
