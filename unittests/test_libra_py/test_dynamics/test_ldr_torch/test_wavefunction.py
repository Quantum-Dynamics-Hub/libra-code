"""Physical spatial projections, global normalization, and saved histories."""

from pathlib import Path
import sys

import pytest
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "src"))
from libra_py.dynamics.ldr_torch.compute import ldr_tdse_solver_multistate, wf_calc_nD


def potential(q, params):
    x = q.sum(dim=0)
    m = params.get("nstates", 2)
    v = torch.zeros(len(x), m, m, dtype=torch.cdouble, device=q.device)
    for i in range(m):
        v[:, i, i] = i + (0.1 + 0.05*i)*x
        for j in range(i + 1, m):
            v[:, i, j] = 0.2/(j-i) + 0.07j*x
            v[:, j, i] = v[:, i, j].conj()
    return v


def make_solver(initialize=True, **updates):
    params = dict(device="cpu", Nstates=2, qgrid=[[-0.8], [0.0], [1.3]],
                  alpha=1.0, q0=[0.0], p0=[0.0], k=[4.0], mass=[1.0],
                  dt=0.13, nsteps=8, potential_fn=potential,
                  initial_representation="diabatic", properties_to_save=[])
    params.update(updates)
    solver = ldr_tdse_solver_multistate(params)
    if initialize:
        solver.initialize_operators()
    return solver


def explicit_sum(solver, points, coefficients):
    """Independent loop over every reference state, center, and electronic component."""
    result = torch.zeros(len(points), solver.nstates, dtype=torch.cdouble)
    for p, q in enumerate(points):
        for n, center in enumerate(solver.qgrid):
            a = solver._alpha_grid[n]
            g = (2*a/torch.pi).pow(0.25).prod() * torch.exp(-(a*(q-center)**2).sum())
            for i in range(solver.nstates):
                for j in range(solver.nstates):
                    result[p, j] += g * coefficients[i*solver.ngrids+n] * solver.eigvecs[n, j, i]
    return result


@pytest.mark.parametrize("ndof,nstates,batch", [(1, 1, 1), (1, 2, 2), (2, 3, 256)])
def test_complex_expansion_matches_explicit_sum(ndof, nstates, batch):
    centers = torch.tensor([[-0.8, 0.2], [0.0, -0.4], [1.3, 0.7]], dtype=torch.double)[:, :ndof]
    widths = torch.tensor([[0.7, 1.2], [1.3, 0.6], [2.1, 0.9]], dtype=torch.double)[:, :ndof]
    points = centers
    s = make_solver(Nstates=nstates, qgrid=centers, alpha=widths, q0=[0.0]*ndof,
                    p0=[0.0]*ndof, k=[4.0]*ndof, mass=[1.0]*ndof,
                    potential_fn_params={"nstates": nstates},
                    psi_batch_size=batch)
    index = torch.arange(s.ndim, dtype=torch.double)
    s.C_curr = torch.cos(index) + 1j*torch.sin(0.7*index)
    s.update_psi_r()
    expected = explicit_sum(s, points, s.C_curr)
    frames = torch.linalg.eigh(potential(points.T, {"nstates": nstates}))[1]
    expected_adi = torch.stack([frames[p].conj().T @ expected[p] for p in range(len(points))])
    assert s.psi_r_dia.shape == s.psi_r_adi.shape == (len(points), nstates)
    torch.testing.assert_close(s.psi_r_dia, expected, atol=2e-12, rtol=2e-12)
    torch.testing.assert_close(s.psi_r_adi, expected_adi, atol=2e-12, rtol=2e-12)
    torch.testing.assert_close(s.psi_r_dia.abs().square().sum(-1),
                               s.psi_r_adi.abs().square().sum(-1))


def test_exact_initial_gaussian_at_qgrid():
    vector = torch.tensor([1.0, 2.0j], dtype=torch.cdouble)
    s = make_solver(initial_electronic_state=vector)
    expected = (2/torch.pi)**0.25 * torch.exp(-s.qgrid.square()) * vector / vector.norm()
    torch.testing.assert_close(s.psi_r_dia, expected, atol=2e-12, rtol=2e-12)


def test_native_projection_and_direct_overlap_density():
    s = make_solver(alpha=[[0.7], [1.3], [2.1]])
    s.C_curr = torch.tensor([0.2j, 0.8, -0.3, 0.4, -0.6j, 0.7j], dtype=torch.cdouble)
    s.update_psi_r()
    for p, q in enumerate(s.qgrid):
        a = s._alpha_grid[:, 0]
        g = (2*a/torch.pi)**0.25 * torch.exp(-a*(q[0]-s.qgrid[:, 0]).square())
        weighted = s.C_curr * g.repeat(s.nstates)
        expected = s.s_elec[torch.arange(s.nstates)*s.ngrids+p] @ weighted
        direct_density = torch.vdot(weighted, s.s_elec @ weighted).real
        torch.testing.assert_close(s.psi_r_adi[p], expected)
        torch.testing.assert_close(s.psi_r_dia[p].abs().square().sum(), direct_density)


def test_reference_phase_changes_preserve_native_densities():
    s = make_solver()
    dia, adi = s.psi_r_dia.clone(), s.psi_r_adi.clone()
    phase = torch.exp(1j*torch.arange(s.ndim, dtype=torch.double).reshape(s.nstates, s.ngrids))
    s.eigvecs *= phase.T[:, None, :]
    s.C_curr *= phase.conj().reshape(-1)
    s.update_psi_r()
    torch.testing.assert_close(s.psi_r_dia, dia)
    torch.testing.assert_close(s.psi_r_adi, phase.T.conj() * adi)
    torch.testing.assert_close(s.psi_r_adi.abs().square(), adi.abs().square())


def test_plain_attributes_update_explicitly_after_coefficient_changes():
    s = make_solver(psi_batch_size=1)
    assert "psi_r_dia" in vars(s) and "psi_r_adi" in vars(s)
    dia = s.psi_r_dia.clone()
    s.C_curr *= 2
    s.update_psi_r()
    torch.testing.assert_close(s.psi_r_dia, 2*dia)
    s.C_curr = 0.5*s.C_curr
    s.update_psi_r()
    torch.testing.assert_close(s.psi_r_dia, dia)


@pytest.mark.parametrize("scheme", ["diagonal", "BAT0", "BAT1", "midpoint"])
def test_history_and_final_snapshot_match_coefficients(tmp_path, scheme):
    s = make_solver(hamiltonian_scheme=scheme, prefix=str(tmp_path / "solution"),
                    qgrid=[[-0.5], [0.2], [1.7], [2.1]], save_every_n_steps=3,
                    properties_to_save=["C_save", "psi_r_dia", "psi_r_adi", "rho_dia", "rho_adi"])
    s.propagate()
    s.save()
    data = torch.load(str(tmp_path / "solution.pt"), weights_only=False)
    assert data["psi_r_dia_all"].shape == data["psi_r_adi_all"].shape == (3, 4, 2)
    assert data["time"] == [0.0, 3*s.dt, 6*s.dt]
    assert s.psi_r_dia_all.shape == (3, 4, 2)
    assert s.rho_dia_all is s.rho_adi_all is None
    assert data["rho_dia_all"] is data["rho_adi_all"] is None
    assert "psi_weights" not in data and not hasattr(s, "psi_weights")
    assert "psi_r_time" not in data
    assert "psi_qgrid" not in data and not hasattr(s, "psi_qgrid")
    for i, c in enumerate(s.C_save):
        torch.testing.assert_close(data["psi_r_dia_all"][i], explicit_sum(s, s.qgrid, c))
        torch.testing.assert_close(data["psi_r_adi_all"][i].abs().square().sum(-1),
                                   data["psi_r_dia_all"][i].abs().square().sum(-1))
    torch.testing.assert_close(data["psi_r_dia"], explicit_sum(s, s.qgrid, s.C_curr))
    assert not torch.allclose(data["psi_r_dia"], data["psi_r_dia_all"][-1])
    saved_initial = s.psi_r_dia_all[0].clone()
    s.C_curr.zero_()
    s.update_psi_r()
    torch.testing.assert_close(s.psi_r_dia, torch.zeros_like(s.psi_r_dia))
    torch.testing.assert_close(s.psi_r_dia_all[0], saved_initial)


def test_default_save_has_final_amplitudes_without_histories(tmp_path):
    s = make_solver(prefix=str(tmp_path / "solution"))
    s.save()
    data = torch.load(str(tmp_path / "solution.pt"), weights_only=False)
    assert data["psi_r_dia_all"] is data["psi_r_adi_all"] is None
    assert data["rho_dia_all"] is data["rho_adi_all"] is None
    torch.testing.assert_close(data["psi_r_dia"], s.psi_r_dia)


def test_legacy_overlaps_support_native_adiabatic_projection(tmp_path):
    full = make_solver()
    params = dict(device="cpu", Nstates=2, qgrid=full.qgrid, alpha=full.alpha,
                  q0=full.q0, p0=full.p0, k=full.k, mass=full.mass,
                  E=full.E, s_elec=full.s_elec, elec_ampl=full.elec_ampl,
                  prefix=str(tmp_path / "legacy"), properties_to_save=["psi_r_adi"])
    legacy = ldr_tdse_solver_multistate(params)
    legacy.initialize_operators()
    assert legacy.psi_r_dia is None
    torch.testing.assert_close(legacy.psi_r_adi, full.psi_r_adi)
    legacy.save_results(0)
    legacy.save()
    with pytest.raises(ValueError, match="require potential_fn"):
        ldr_tdse_solver_multistate(dict(params, properties_to_save=["psi_r_dia"]))


@pytest.mark.parametrize("batch_size", [0, True, 1.5])
def test_invalid_batch_size(batch_size):
    with pytest.raises(ValueError, match="psi_batch_size"):
        make_solver(psi_batch_size=batch_size)


def test_spatial_density_integral_matches_analytic_gaussian_overlap():
    qgrid = torch.linspace(-8, 8, 161, dtype=torch.double)[:, None]
    active = [72, 80, 93]  # Only three centers contribute; all qgrid points are sampled.
    alpha = torch.ones_like(qgrid)
    alpha[active, 0] = torch.tensor([0.7, 1.3, 2.1], dtype=torch.double)
    # This is a quadrature check, not a propagation test. Avoid solving for
    # coefficients in the deliberately dense, linearly dependent basis.
    s = make_solver(initialize=False, qgrid=qgrid, alpha=alpha)
    s.chi_overlap()
    s.build_compound_overlap()
    s.C_curr.reshape(2, s.ngrids)[:, active] = torch.tensor(
        [[1.0, 0.2j, -0.4], [0.1j, 0.3, 0.6j]], dtype=torch.cdouble)
    s.update_psi_r()
    # Analytic integral over all space using nuclear Gaussian overlaps.
    vectors = {n: s.eigvecs[n] @ s.C_curr.reshape(2, s.ngrids)[:, n] for n in active}
    expected = sum(s.s_nucl[n, m] * torch.outer(vectors[n], vectors[m].conj())
                   for n in active for m in active)
    integrated = torch.trapezoid(
        s.psi_r_dia[:, :, None] * s.psi_r_dia[:, None, :].conj(), qgrid[:, 0], dim=0)
    torch.testing.assert_close(integrated, expected, atol=2e-12, rtol=2e-12)
    norm = torch.vdot(s.C_curr, s.S @ s.C_curr)
    for psi in (s.psi_r_dia, s.psi_r_adi):
        raw = psi.abs().square().sum(-1)
        torch.testing.assert_close(torch.trapezoid(raw, qgrid[:, 0]), norm.real)
        torch.testing.assert_close(torch.trapezoid(raw / norm.real, qgrid[:, 0]),
                                   torch.ones_like(norm.real))
    assert s.psi_r_dia_all is s.psi_r_adi_all is None


def test_inference_tensors_reconstruct_after_explicit_update():
    s = make_solver()
    with torch.inference_mode():
        s.C_curr = s.C_curr.clone()
        initial = s.psi_r_dia.clone()
        s.C_curr *= 2
        s.update_psi_r()
        torch.testing.assert_close(s.psi_r_dia, 2*initial)


@pytest.mark.parametrize("ndof,nstates,batch", [(1, 1, 1), (2, 2, 3), (3, 3, 256)])
def test_wf_calc_nD_matches_full_complex_sum_on_independent_grid(tmp_path, ndof, nstates, batch):
    centers = torch.tensor([[-0.8, 0.2, 0.3], [0.0, -0.4, 0.6], [1.3, 0.7, -0.2]],
                           dtype=torch.double)[:, :ndof]
    widths = torch.tensor([[0.7, 1.2, 0.6], [1.3, 0.6, 0.9], [2.1, 0.9, 1.1]],
                          dtype=torch.double)[:, :ndof]
    s = make_solver(initialize=False, Nstates=nstates, qgrid=centers, alpha=widths,
                    q0=[0.0]*ndof, p0=[0.0]*ndof, k=[4.0]*ndof, mass=[1.0]*ndof,
                    potential_fn_params={"nstates": nstates})
    s.chi_overlap()
    s.build_compound_overlap()
    index = torch.arange(s.ndim, dtype=torch.double)
    c = torch.stack([torch.cos(index) + 1j*torch.sin(0.7*index),
                     torch.sin(index + 0.4) - 0.3j*torch.cos(index)])
    norms = [torch.vdot(row, s.S @ row).real.sqrt() for row in c]
    filename = tmp_path / "snapshots.pt"
    torch.save(dict(qgrid=s.qgrid, alpha=s.alpha, eigvecs=s.eigvecs,
                    C_save=list(c), time=[0., 2.5], norm=norms), filename)
    axes = [torch.linspace(-0.93, 1.17, 4+i, dtype=torch.double) for i in range(ndof)]
    result = wf_calc_nD(filename, axes, potential, {"nstates": nstates},
                        times=[2.5, 0.], batch_size=batch)
    points = torch.stack(torch.meshgrid(*axes, indexing="ij"), dim=-1).reshape(-1, ndof)
    expected = torch.stack([explicit_sum(s, points, c[i]) / norms[i] for i in [1, 0]])
    local_frames = torch.linalg.eigh(potential(points.T, {"nstates": nstates}))[1]
    expected_adi = torch.stack([
        torch.stack([local_frames[p].conj().T @ row[p] for p in range(len(points))])
        for row in expected
    ])
    shape = (2, *(len(a) for a in axes), nstates)
    torch.testing.assert_close(result["time"], torch.tensor([2.5, 0.], dtype=torch.double))
    torch.testing.assert_close(result["psi_r_dia"], expected.reshape(shape), atol=2e-12, rtol=2e-12)
    torch.testing.assert_close(result["psi_r_adi"], expected_adi.reshape(shape), atol=2e-12, rtol=2e-12)
    torch.testing.assert_close(result["rho_dia"], result["psi_r_dia"].abs().square())
    torch.testing.assert_close(result["rho_adi"], result["psi_r_adi"].abs().square())
    torch.testing.assert_close(result["rho_total"], result["rho_adi"].sum(-1))


def test_wf_calc_nD_reads_actual_saved_history_and_preserves_reference_phases(tmp_path):
    s = make_solver(prefix=str(tmp_path / "solution"), save_every_n_steps=1,
                    properties_to_save=["time", "norm", "C_save", "psi_r_dia", "psi_r_adi"])
    s.propagate()
    s.save()
    filename = tmp_path / "solution.pt"
    result = wf_calc_nD(filename, [s.qgrid[:, 0]], potential, normalize=False, batch_size=1)
    torch.testing.assert_close(result["psi_r_dia"], s.psi_r_dia_all)
    torch.testing.assert_close(result["rho_adi"], s.psi_r_adi_all.abs().square())

    saved = torch.load(filename, weights_only=True)
    phase = torch.exp(1j*torch.arange(s.ndim, dtype=torch.double).reshape(s.nstates, s.ngrids))
    saved["eigvecs"] = saved["eigvecs"] * phase.T[:, None, :]
    saved["C_save"] = [row * phase.conj().reshape(-1) for row in saved["C_save"]]
    torch.save(saved, filename)
    phased = wf_calc_nD(filename, [s.qgrid[:, 0]], potential, normalize=False, batch_size=2)
    torch.testing.assert_close(phased["psi_r_dia"], result["psi_r_dia"])
    torch.testing.assert_close(phased["rho_adi"], result["rho_adi"])


def test_wf_calc_nD_normalization_matches_analytic_gaussian(tmp_path):
    s = make_solver(qgrid=[[0.]], alpha=1., prefix=str(tmp_path / "gaussian"),
                    initial_electronic_state=[1., 2.j], properties_to_save=["time", "norm", "C_save"])
    s.C_curr *= 3.
    s.save_results(0)
    s.save()
    x = torch.linspace(-8., 8., 1001, dtype=torch.double)
    filename = tmp_path / "gaussian.pt"
    result = wf_calc_nD(filename, [x], potential)
    expected = (2/torch.pi)**0.5 * torch.exp(-2*x.square())
    torch.testing.assert_close(result["rho_total"][0], expected, atol=2e-12, rtol=2e-12)
    torch.testing.assert_close(torch.trapezoid(result["rho_total"][0], x), torch.tensor(1., dtype=torch.double))
    torch.testing.assert_close(result["rho_adi"].sum(-1), result["rho_total"])

    saved = torch.load(filename, weights_only=True)
    saved.pop("norm")
    torch.save(saved, filename)
    raw = wf_calc_nD(filename, [x], normalize=False)
    torch.testing.assert_close(raw["rho_total"][0], 9*expected, atol=2e-12, rtol=2e-12)
    assert raw["psi_r_adi"] is None and raw["rho_adi"] is None
    with pytest.raises(ValueError, match="requires a saved norm"):
        wf_calc_nD(filename, [x])


@pytest.mark.parametrize("problem,match", [
    ("times", "requested time"), ("coefficients", "time and C_save"),
    ("frames", "saved eigvecs"), ("axes", "axes must contain"),
    ("batch", "batch_size"), ("norm", "norms must be finite and positive"),
])
def test_wf_calc_nD_rejects_missing_or_inconsistent_data(tmp_path, problem, match):
    saved = dict(qgrid=torch.zeros(1, 1), alpha=1., eigvecs=torch.ones(1, 1, 1),
                 C_save=[torch.ones(1)], time=[0.], norm=[1.])
    options = dict(axes=[[0.]], times=[0.])
    if problem == "times":
        options["times"] = [0.5]
    elif problem == "coefficients":
        saved["C_save"] = []
    elif problem == "frames":
        saved["eigvecs"] = None
    elif problem == "axes":
        options["axes"] = [[0.], [0.]]
    elif problem == "batch":
        options["batch_size"] = 0
    elif problem == "norm":
        saved["norm"] = [0.]
    filename = tmp_path / "invalid.pt"
    torch.save(saved, filename)
    with pytest.raises(ValueError, match=match):
        wf_calc_nD(filename, **options)
