"""Check physical initial states selected in either electronic representation."""

from pathlib import Path
import sys

import pytest
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "src"))
from libra_py.dynamics.ldr_torch.compute import ldr_tdse_solver_multistate


def potential(q, params):
    """Complex, coordinate-dependent mixing distinguishes both choices."""
    x = q[0]
    matrix = torch.zeros(x.numel(), 2, 2, dtype=torch.cdouble, device=q.device)
    matrix[:, 0, 0] = 0.2 + 0.1*x
    matrix[:, 1, 1] = 1.0 - 0.2*x
    matrix[:, 0, 1] = 0.3 + 0.1j*x
    matrix[:, 1, 0] = matrix[:, 0, 1].conj()
    return matrix


def solver(**updates):
    params = dict(device="cpu", Nstates=2, initial_state_index=0,
                  qgrid=[[-0.7], [0.0], [0.8]], alpha=[1.0],
                  q0=[0.0], p0=[0.0], k=[4.0], mass=[1.0], dt=0.02,
                  potential_fn=potential)
    params.update(updates)
    result = ldr_tdse_solver_multistate(params)
    result.initialize_operators()
    return result


def reconstruct(obj, points, coefficients=None):
    """Evaluate the physical wavefunction in the model's diabatic basis."""
    if coefficients is None:
        coefficients = obj.C0
    a = obj._alpha_grid[:, 0]
    gaussian = (2*a/torch.pi)**0.25 * torch.exp(-a*(points[:, None]-obj.qgrid[:, 0])**2)
    return torch.einsum("in,pn,ndi->pd", coefficients.reshape(obj.nstates, obj.ngrids),
                        gaussian.to(torch.cdouble), obj.eigvecs)


@pytest.mark.parametrize("index", [0, 1])
@pytest.mark.parametrize("scheme", ["diagonal", "BAT0", "BAT1", "midpoint"])
def test_default_and_explicit_diabatic_selection_project_the_requested_state(index, scheme):
    default = solver(initial_state_index=index, hamiltonian_scheme=scheme)
    selected = solver(initial_representation="diabatic", initial_state_index=index,
                      hamiltonian_scheme=scheme)
    vector = [0.0, 0.0]
    vector[index] = 1.0
    explicit = solver(initial_electronic_state=vector, hamiltonian_scheme=scheme)
    # The Gaussian is exactly one nuclear basis function, so projection must
    # recover the target vector throughout space, including complex phases.
    points = torch.linspace(-2.0, 2.0, 51, dtype=torch.double)
    expected = ((2/torch.pi)**0.25*torch.exp(-points**2))[:, None] * torch.tensor(vector)
    torch.testing.assert_close(reconstruct(selected, points), expected.to(torch.cdouble),
                               atol=2e-12, rtol=2e-12)
    torch.testing.assert_close(default.C0, selected.C0, atol=0, rtol=0)
    torch.testing.assert_close(selected.C0, explicit.C0, atol=2e-12, rtol=2e-12)
    c1, c2 = selected.C0.clone(), explicit.C0.clone()
    for _ in range(5):
        c1, c2 = selected.U @ c1, explicit.U @ c2
    torch.testing.assert_close(reconstruct(selected, points, c1), reconstruct(explicit, points, c2),
                               atol=2e-12, rtol=2e-12)


@pytest.mark.parametrize("index", [0, 1])
def test_adiabatic_selection_holds_the_state_at_q0_fixed(index):
    selected = solver(initial_representation="adiabatic", initial_state_index=index)
    _, eigenvectors = torch.linalg.eigh(potential(torch.zeros(1, 1, dtype=torch.double), {})[0])
    # Evaluate the same electronic vector at all coordinates, rather than
    # following the coordinate-dependent eigenvectors of the local potential.
    points = torch.linspace(-2.0, 2.0, 51, dtype=torch.double)
    expected = ((2/torch.pi)**0.25*torch.exp(-points**2))[:, None] * eigenvectors[:, index]
    torch.testing.assert_close(reconstruct(selected, points), expected, atol=2e-12, rtol=2e-12)


@pytest.mark.parametrize("representation", ["adiabatic", "diabatic"])
def test_explicit_complex_vector_overrides_representation_and_index(representation):
    vector = torch.tensor([1.0, 2.0j], dtype=torch.cdouble)
    selected = solver(initial_representation=representation, initial_state_index=1,
                      initial_electronic_state=vector)
    points = torch.linspace(-2.0, 2.0, 51, dtype=torch.double)
    expected = ((2/torch.pi)**0.25*torch.exp(-points**2))[:, None] * vector / torch.linalg.vector_norm(vector)
    torch.testing.assert_close(reconstruct(selected, points), expected, atol=2e-12, rtol=2e-12)


@pytest.mark.parametrize("representation", ["invalid", None])
def test_invalid_representation_is_rejected(representation):
    with pytest.raises(ValueError, match="initial_representation must be"):
        solver(initial_representation=representation)


def test_precomputed_inputs_do_not_silently_ignore_representation():
    with pytest.raises(ValueError, match="initial_representation requires potential_fn"):
        ldr_tdse_solver_multistate(dict(initial_representation="diabatic"))
