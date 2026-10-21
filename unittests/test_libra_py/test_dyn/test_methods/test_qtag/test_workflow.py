# Copyright (C) 2022-2026 Matthew Dutra, Alexey V. Akimov
# SPDX-License-Identifier: GPL-3.0-or-later

"""Initialization, moving/frozen basis evolution, persistence and plotting tests."""
import numpy as np
import pytest
from libra_py.dyn.methods.qtag import core, compute, initialize, save, plot
from libra_py.dyn.methods.qtag._backend import zeros, to_numpy, asarray


def harmonic(q, params):
    """Two displaced harmonic wells with constant Hermitian coupling."""
    H = zeros((q.shape[1], 2, 2), q, True)
    dh = zeros((q.shape[1], 1, 2, 2), q, True)
    hh = zeros((q.shape[1], 1, 1, 2, 2), q, True)
    H[:, 0, 0] = .5*q[0]**2
    H[:, 1, 1] = .5*(q[0]-.7)**2+.1
    H[:, 0, 1] = H[:, 1, 0] = .08
    dh[:, 0, 0, 0] = q[0]
    dh[:, 0, 1, 1] = q[0]-.7
    hh[:, 0, 0, 0, 0] = hh[:, 0, 0, 1, 1] = 1
    return core.HamiltonianData(H, dh, hh)


def setup(backend='numpy', **options):
    if backend == 'torch':
        pytest.importorskip('torch')
    pars = dict(backend=backend, wfc_q0=[-.7], wfc_p0=[.4], wfc_a0=[.8],
                wfc_s0=[.3], states=[0, 1], grid_dims=[7], rho_cut=-1,
                grid_min=[-2.2], grid_max=[1.2], alp_scl=[3.])
    pars.update(options)
    n, q, p, a, s, labels = initialize.gbfs(pars)
    q0, p0, a0, s0 = initialize.wfc_centers(pars)
    c = initialize.coeffs(q0, p0, a0, s0, [0], q, p, a, s, labels)
    return q, p, a, s, labels, c, asarray([1.], q)


def run(backend='numpy', moving=True, nsteps=5):
    args = setup(backend)
    params = dict(dt=.005, nsteps=nsteps, states=[0, 1], target_states=[0, 1],
                  q_update_method=int(moving), p_update_method=int(moving),
                  a_update_method=int(moving), qtag_pot_approx_method=1,
                  decpl_den=.01, q_sync_method=1, p_sync_method=1, a_sync_method=1)
    return compute.run_qtag(*args, params, harmonic)


@pytest.mark.parametrize('backend', ['numpy', 'torch'])
def test_initialization(backend):
    """Projection normalization, phase preservation and deterministic sampling."""
    q, p, a, s, states, c, im = setup(backend)
    overlap = compute.time_overlap(q, p, a, s, states, q, p, a, s, states)
    pops = compute.qtag_pops(states, c, overlap, [0, 1])
    np.testing.assert_allclose(to_numpy(pops), [1., 0.], atol=1e-12)
    np.testing.assert_allclose(to_numpy(s), .3)
    first = initialize.gaussian(dict(backend=backend, seed=11, rho_cut=.8,
                                     wfc_q0=[0., 1.], grid_dims=[2, 3]))
    second = initialize.gaussian(dict(backend=backend, seed=11, rho_cut=.8,
                                      wfc_q0=[0., 1.], grid_dims=[2, 3]))
    assert first[0] == 12
    np.testing.assert_allclose(to_numpy(first[1]), to_numpy(second[1]))
    assert first[1].shape == (2, 12)


@pytest.mark.parametrize('backend', ['numpy', 'torch'])
def test_frozen_norm_energy_and_times(backend):
    """A fixed Hermitian basis conserves norm and energy at every saved time."""
    result = run(backend, moving=False, nsteps=8)
    h = result.history
    np.testing.assert_allclose(to_numpy(h['norm']), 1., atol=2e-11)
    np.testing.assert_allclose(to_numpy(h['Etot']), to_numpy(h['Etot'][0]), atol=2e-11)
    np.testing.assert_allclose(to_numpy(h['time']), .005*np.arange(9))
    np.testing.assert_allclose(to_numpy(h['pops']).sum(axis=1), to_numpy(h['norm']), atol=1e-12)
    assert float(h['pops'][-1, 1]) > 0
    # Final populations must describe the saved final coefficients/basis.
    S, H = core.qtag_hamiltonian_and_overlap(result.q, result.p, result.alp, result.s,
                                            result.states, asarray([1.], result.q), harmonic,
                                            dyn_params={'qtag_pot_approx_method': 1})
    np.testing.assert_allclose(to_numpy(h['Etot'][-1]), to_numpy(compute.qtag_energy(result.coeff, H)))
    assert len(tuple(result)) == 6


def test_backend_parity_and_no_input_mutation():
    """Compare all moving-basis snapshots on NumPy and PyTorch."""
    pytest.importorskip('torch')
    numpy_result, torch_result = run(), run('torch')
    for key in numpy_result.history:
        np.testing.assert_allclose(to_numpy(torch_result.history[key]), numpy_result.history[key],
                                   rtol=2e-8, atol=3e-9, err_msg=key)
    args = setup()
    before = [x.copy() for x in args[:4]]
    compute.run_qtag(*args, {'dt': .005, 'nsteps': 2}, harmonic)
    for actual, expected in zip(args[:4], before):
        np.testing.assert_array_equal(actual, expected)


@pytest.mark.parametrize('backend', ['numpy', 'torch'])
def test_linear_fit_and_nodes(backend):
    """One Gaussian has exact affine log-gradient; nodes never divide by zero."""
    q, p, a, s, labels, c, im = setup(backend)
    # A single nonzero coefficient makes the wavefunction one Gaussian even
    # though its logarithmic derivative is sampled at all basis centers.
    c = zeros((q.shape[1],), q, True)
    c[0] = 1
    mom, real, gmom, greal = compute.qtag_momentum(q, p, a, s, c)
    np.testing.assert_allclose(to_numpy(mom), to_numpy(p[:, :1])+np.zeros(q.shape), atol=1e-10)
    np.testing.assert_allclose(to_numpy(real), to_numpy(-a[:, :1]*(q-q[:, :1])), atol=1e-10)
    np.testing.assert_allclose(to_numpy(gmom), 0., atol=1e-10)
    np.testing.assert_allclose(to_numpy(greal), -float(a[0, 0]), atol=1e-10)
    ratios, slopes, rho = core.momentum_fields(q, p, a, s, 0*c)
    np.testing.assert_array_equal(to_numpy(ratios), 0)
    np.testing.assert_array_equal(to_numpy(slopes), 0)
    np.testing.assert_array_equal(to_numpy(rho), 0)


@pytest.mark.parametrize('backend', ['numpy', 'torch'])
def test_save_load_and_reconstruction(tmp_path, backend):
    """NPZ retains phases, coefficients, shapes and all saved time samples."""
    result = run(backend)
    file = save.save_qtag_data(result, tmp_path/'trajectory.npz')
    loaded = save.load_qtag_data(file, backend)
    for key in loaded:
        np.testing.assert_array_equal(to_numpy(loaded[key]), to_numpy(result.history[key]))
    params = dict(xmin=[-3.], xmax=[3.], npoints=[31], snaps=[0, 5], which_states=[0, 1])
    original, restored = compute.wf_calc_nD(result, params), compute.wf_calc_nD(loaded, params)
    np.testing.assert_allclose(to_numpy(original['psi']), to_numpy(restored['psi']))
    assert restored['density'].shape == (2, 2, 31)
    # Direct evaluation checks that reconstruction uses the saved phase.
    c = compute.coefficient_matrix(result.history['coeffs'][0], result.states, [0, 1])
    direct = core.qtag_psi(original['points'], result.history['q'][0], result.history['p'][0],
                          result.history['a'][0], result.history['s'][0], c)
    np.testing.assert_allclose(to_numpy(original['psi'][0]), to_numpy(direct))
    assert save.save_text(result, tmp_path).exists()


def test_hdf5_and_individual_plots(tmp_path):
    """Exercise explicit output boundaries and individual headless PNG files."""
    pytest.importorskip('matplotlib')
    result = run(nsteps=1)
    paths = plot.energy_and_pops(result, tmp_path)+plot.trajectories(result, tmp_path)
    data = compute.wf_calc_nD(result, dict(xmin=[-3.], xmax=[3.], npoints=[31], snaps=[0, 1]))
    paths += plot.wf_plot(data, tmp_path)
    assert len(paths) == 10
    assert all(p.read_bytes().startswith(b'\x89PNG') for p in paths)
    # 2D plotting uses one separate file per state/snapshot.
    data2 = dict(axes=[np.arange(3), np.arange(4)], density=np.ones((1, 2, 3, 4)),
                 snaps=[0], states=[0, 1], times=[0.])
    assert len(plot.plot_wf_2D(data2, tmp_path)) == 2
    pytest.importorskip('h5py')
    filename = save.save_qtag_data(result, tmp_path/'trajectory.h5')
    loaded = save.load_qtag_data(filename)
    np.testing.assert_array_equal(loaded['coeffs'], result.history['coeffs'])


def test_bad_controls_and_regularization():
    args = setup()
    with pytest.raises(NotImplementedError, match='phase'):
        compute.run_qtag(*args, dict(dt=.01, nsteps=1, s_update_method=1), harmonic)
    with pytest.raises(ValueError, match='dt'):
        compute.run_qtag(*args, dict(dt=-.01), harmonic)
    # Exact duplicate basis functions require truncation, not a matrix inverse.
    S = np.ones((2, 2), dtype=complex)
    c = np.array([.5, .5], dtype=complex)
    evolved = core.propagate_electronic_qtag(.2, c, .3*S, S)
    np.testing.assert_allclose(evolved, np.exp(-.06j)*c, atol=1e-12)


def test_torch_end_to_end_gradient():
    """Differentiate a short moving-basis trajectory with nondegenerate metric."""
    torch = pytest.importorskip('torch')
    args = setup('torch', states=[0], grid_dims=[4])
    q, p, a, s, labels, c, im = args
    q = q.detach().clone().requires_grad_(True)
    pars = dict(dt=.002, nsteps=2, states=[0], target_states=[0],
                q_update_method=1, p_update_method=1, a_update_method=1,
                qtag_pot_approx_method=1)
    result = compute.run_qtag(q, p, a, s, labels, c, im, pars, harmonic)
    objective = result.history['Etot'][-1]
    objective.backward()
    assert torch.isfinite(q.grad).all()
    eps = 1e-6
    plus, minus = q.detach().clone(), q.detach().clone()
    plus[0, 1] += eps
    minus[0, 1] -= eps
    with torch.no_grad():
        eplus = compute.run_qtag(plus, p, a, s, labels, c, im, pars, harmonic).history['Etot'][-1]
        eminus = compute.run_qtag(minus, p, a, s, labels, c, im, pars, harmonic).history['Etot'][-1]
    np.testing.assert_allclose(to_numpy(q.grad[0, 1]), to_numpy((eplus-eminus)/(2*eps)), rtol=2e-4, atol=2e-6)


def test_without_cpp_or_sympy():
    """A fresh process can initialize and propagate with both modules blocked."""
    from pathlib import Path
    import subprocess
    import sys
    src = Path(__file__).resolve().parents[5]/'src'
    program = f'''
import sys
sys.path.insert(0, {str(src)!r})
sys.modules['liblibra_core'] = None
sys.modules['cyglibra_core'] = None
sys.modules['sympy'] = None
import numpy as np
from libra_py.dyn.methods.qtag import core, initialize, compute
n,q,p,a,s,labels = initialize.gbfs(dict(states=[0], grid_dims=[3]))
q0,p0,a0,s0 = initialize.wfc_centers(dict(states=[0]))
c = initialize.coeffs(q0,p0,a0,s0,[0],q,p,a,s,labels)
def model(q, params):
    n=q.shape[1]
    return core.HamiltonianData(np.zeros((n,1,1)), np.zeros((n,1,1,1)))
r = compute.run_qtag(q,p,a,s,labels,c,np.ones(1),dict(nsteps=1,dt=.01),model)
assert np.isfinite(r.history['Etot']).all()
'''
    completed = subprocess.run([sys.executable, '-c', program], text=True, capture_output=True)
    assert completed.returncode == 0, completed.stderr
