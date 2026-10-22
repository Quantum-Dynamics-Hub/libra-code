"""Integrated QTAG regression against the frozen test-only QTAG implementation.

Compare complete histories rather than only final populations. Extra tests
exercise ownership and topology, which array-only reference tests cannot cover.
"""
import importlib.util
import os
from pathlib import Path
import subprocess
import sys
import numpy as np
import pytest
from libra_py.dyn.methods.qtag import QTAGEngine, run_qtag
from libra_py.dyn.initialization.ensembles import make_gaussian_ensemble
from libra_py.dyn.models import Holstein5Model, FerrettiModel
from libra_py.dyn.backends.arrays import to_numpy
from libra_py.dyn.gwp.storage import read_packets
from libra_py.dyn.observables.gaussian import reconstruct_density
from qtag_reference import compute as reference, initialize as ref_init

ROOT = Path(__file__).resolve().parents[5]


def settings(backend='numpy', case='holstein', method=0):
    """Small nontrivial bases with nonzero momenta/phases, both model dimensions."""
    if backend == 'torch':
        pytest.importorskip('torch')
    if case == 'holstein':
        spec = dict(wfc_q0=[-.7], wfc_p0=[.4], wfc_a0=[.8], wfc_s0=[.3],
                    grid_dims=[7], grid_min=[-2.2], grid_max=[1.2], alp_scl=[3.], rho_cut=-1)
        model = Holstein5Model(backend=backend, params=dict(E_n=[0., .1], x_n=[0., .7],
            k_n=[1., 1.2], V=[[0., .08], [.08, 0.]], alpha=[[0., .3], [.3, 0.]],
            x_nm=[[0., .2], [.2, 0.]], ex_cpl_A=.08, ex_cpl_B=.3, ex_cpl_C=.2))
        mass = [1.]
    else:
        spec = dict(wfc_q0=[0., .7], wfc_p0=[.3, -.1], wfc_a0=[.8, 1.1],
                    wfc_s0=[.1, -.2], grid_dims=[3, 3], alp_scl=[3., 2.], rho_cut=.01)
        model = FerrettiModel(backend=backend, params=dict(Ky=.01))
        mass = [.5, .7]
    spec.update(backend=backend, states=[0, 1])
    pars = dict(dt=.005, nsteps=6, states=[0, 1], target_states=[0, 1],
                qtag_pot_approx_method=method, a_update_method=1,
                q_sync_method=1, p_sync_method=0, a_sync_method=0,
                decpl_den=.01, linfit_beta=1e-3)
    return spec, mass, model, pars


def oracle(spec, mass, model, pars, manager=None, tid=0):
    """Run the preserved array workflow from the same live packet state."""
    manager = manager or make_gaussian_ensemble(spec, mass)
    slots, packet, labels, coeff, im = read_packets(manager.storage, manager.trajectories[tid])
    return reference.run_qtag(*packet, labels, coeff, im, pars, model, model.params)


def compare(actual, expected, atol=3e-9):
    for key in ('time', 'Etot', 'dEtot', 'norm', 'pops', 'coeffs'):
        np.testing.assert_allclose(to_numpy(actual[key]), to_numpy(expected.history[key]),
                                   atol=atol, rtol=2e-8, err_msg=key)
    for key, old in [('q', 'q'), ('p', 'p'), ('gaussian_alpha', 'a'), ('gaussian_phase', 's')]:
        np.testing.assert_allclose(to_numpy(actual[key]).swapaxes(1, 2),
                                   to_numpy(expected.history[old]), atol=atol, rtol=2e-8)


@pytest.mark.parametrize('backend', ['numpy', 'torch'])
@pytest.mark.parametrize('case,method', [('holstein', k) for k in range(4)]+[('ferretti',0),('ferretti',1)])
def test_reference_history_and_density(backend, case, method):
    """All supported model/quadrature combinations match at every short step."""
    spec, mass, model, pars = settings(backend, case, method)
    expected = oracle(spec, mass, model, pars)
    manager = make_gaussian_ensemble(spec, mass)
    actual = run_qtag(manager, model, pars).histories[0]
    compare(actual, expected)
    d = len(mass)
    grid = dict(xmin=[-3.]*d, xmax=[3.]*d, npoints=[11]*d, snaps=[0, 6])
    density = reconstruct_density(actual, **grid)
    ref_density = reference.wf_calc_nD(expected, grid)
    np.testing.assert_allclose(to_numpy(density['density']), to_numpy(ref_density['density']), atol=3e-9)
    tbf = manager.active_tbfs(0)[0]
    np.testing.assert_array_equal(to_numpy(tbf.view.gaussian_alpha), to_numpy(actual['gaussian_alpha'][-1, 0]))
    assert manager.storage.timestep == 6
    assert manager.trajectories[0].norm == pytest.approx(float(actual['norm'][-1]))


@pytest.mark.parametrize('backend', ['numpy', 'torch'])
def test_independent_trajectories_pruning_and_compaction(backend):
    """Noncontiguous active slots, unequal basis counts and separate wavefunctions."""
    spec, mass, model, pars = settings(backend)
    pars['q_sync_method'] = 0  # unequal surface counts cannot synchronize by index
    other = dict(spec, grid_dims=[5], wfc_q0=[.6])
    manager = make_gaussian_ensemble([spec, other], mass)
    # Remove one packet on state 0. Both surfaces still have basis coverage.
    manager.deactivate_tbf(manager.trajectories[0].tbf_object_by_slot[2])
    expected0 = oracle(spec, mass, model, pars, manager, 0)
    expected1 = oracle(other, mass, model, pars, manager, 1)
    out = run_qtag(manager, model, pars)
    compare(out.histories[0], expected0)
    compare(out.histories[1], expected1)
    assert 2 not in to_numpy(out.histories[0]['slot_ids'][0])
    assert manager.storage.timestep == pars['nsteps']
    # Global compaction must preserve Gaussian fields and registered views.
    manager.deactivate_tbf(manager.trajectories[1].tbf_object_by_slot[2])
    before = [read_packets(manager.storage, t)[1] for t in manager.trajectories.values()]
    manager.compact_storage()
    for t, old in zip(manager.trajectories.values(), before):
        for a, b in zip(read_packets(manager.storage, t)[1], old):
            np.testing.assert_array_equal(to_numpy(a), to_numpy(b))
    frozen = to_numpy(manager.storage.q[1]).copy()
    manager.deactivate_trajectory(1, deactivate_tbfs=False)
    QTAGEngine(manager, model, pars).step()
    np.testing.assert_array_equal(to_numpy(manager.storage.q[1]), frozen)


def test_torch_gradient_and_storage_precision():
    """A scalar loss differentiates through shared storage and two moving steps."""
    torch = pytest.importorskip('torch')
    spec, mass, model, pars = settings('torch')
    # Identical bases on two surfaces give exactly degenerate overlap
    # eigenvalues; eigenvector derivatives there are not well-defined.
    spec.update(states=[0], grid_dims=[4])
    pars.update(nsteps=2, momentum_fit='direct', states=[0], target_states=[0])
    manager = make_gaussian_ensemble(spec, mass, nstates=2)
    initial = manager.storage.q.detach().clone().requires_grad_(True)
    manager.storage.q = initial
    result = run_qtag(manager, model, pars)
    loss = result.histories[0]['Etot'][-1]
    loss.backward()
    assert bool(torch.isfinite(initial.grad).all())
    assert float(abs(initial.grad).sum()) > 0
    values = []
    eps = 1e-6
    for sign in (1, -1):
        perturbed = make_gaussian_ensemble(spec, mass, nstates=2)
        perturbed.storage.q = initial.detach().clone()
        perturbed.storage.q[0, 1, 0] += sign*eps
        values.append(float(run_qtag(perturbed, model, pars).histories[0]['Etot'][-1]))
    assert float(initial.grad[0, 1, 0]) == pytest.approx((values[0]-values[1])/(2*eps), rel=2e-4, abs=2e-6)
    for backend in ('numpy', 'torch'):
        spec['backend'], spec['dtype'] = backend, 'float32'
        storage = make_gaussian_ensemble(spec, mass).storage
        assert str(storage.q.dtype).endswith('float32')
        assert str(storage.ampl_dia.dtype).endswith('complex64')


def test_no_reference_or_sympy_runtime_dependency():
    """A fresh subprocess blocks reference/C++/SymPy imports during propagation."""
    code = '''
import importlib.abc, sys
class Block(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.startswith(('qtag_reference','sympy','liblibra_core')):
            raise ImportError('blocked: '+fullname)
sys.meta_path.insert(0, Block())
from libra_py.dyn.initialization.ensembles import make_gaussian_ensemble
from libra_py.dyn.methods.qtag import run_qtag
from libra_py.dyn.models import FerrettiModel
m=make_gaussian_ensemble(dict(wfc_q0=[0.,1.],grid_dims=[2,2],alp_scl=[6.,6.]),[1.,1.])
r=run_qtag(m,FerrettiModel(),dict(dt=.001,nsteps=2))
assert r.storage.timestep==2
'''
    env = dict(os.environ, PYTHONPATH=str(ROOT/'src'))
    subprocess.run([sys.executable, '-c', code], env=env, check=True, capture_output=True, text=True)


@pytest.mark.skipif(os.environ.get('QTAG_FULL_TESTS') != '1', reason='set QTAG_FULL_TESTS=1 for full tutorial parity')
@pytest.mark.parametrize('backend', ['numpy', 'torch'])
@pytest.mark.parametrize('case', ['holstein', 'ferretti'])
def test_full_tutorial(backend, case):
    """Full histories and snapshots agree with the preserved method implementation."""
    if backend == 'torch':
        torch = pytest.importorskip('torch')
        torch.set_num_threads(1)
    file = ROOT/'examples/libra_py/dyn/methods/qtag/1_basics/run.py'
    spec = importlib.util.spec_from_file_location('integrated_qtag_example', file)
    example = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(example)
    initial, mass, model, pars, grid = example.configuration(case, backend)
    expected = oracle(initial, mass, model, pars)
    result = run_qtag(make_gaussian_ensemble(initial, mass), model, pars)
    h = result.histories[0]
    compare(h, expected, atol=2e-6)
    n = pars['nsteps']
    snaps = [1, n//4, n//2, 3*n//4, n]
    density = reconstruct_density(h, **grid, snaps=snaps)
    reference_density = reference.wf_calc_nD(expected, dict(grid, snaps=snaps))
    np.testing.assert_allclose(to_numpy(density['density']), to_numpy(reference_density['density']), atol=2e-6, rtol=2e-5)


@pytest.mark.parametrize('backend', ['numpy', 'torch'])
def test_noncontiguous_electronic_labels_and_slot_order(backend):
    """Storage electronic indices need not be contiguous or sorted in the basis."""
    spec, mass, _, pars = settings(backend)
    spec.update(states=[2, 0], initial_state=2, grid_dims=[4])
    model = Holstein5Model(backend=backend, params=dict(E_n=[0., .1, .2], x_n=[0., .3, .6],
        k_n=[1., 1.1, 1.2], V=[[0., .08, .08], [.08, 0., .08], [.08, .08, 0.]],
        alpha=[[.3]*3]*3, x_nm=[[.2]*3]*3))
    pars.update(states=[2, 0], target_states=[2, 1, 0])
    manager = make_gaussian_ensemble(spec, mass)
    manager.trajectories[0].tbf_ids = [0, 4, 1, 5, 2, 6, 3, 7]
    expected = oracle(spec, mass, model, pars, manager)
    result = run_qtag(manager, model, pars)
    compare(result.histories[0], expected)
    np.testing.assert_array_equal(to_numpy(result.histories[0]['pops'][:, 1]), 0)


def test_ensemble_model_timestep_and_invalid_inputs():
    """All trajectories see the same old/new model times during one split step."""
    specs = [dict(states=[0], grid_dims=[3], wfc_q0=[q]) for q in (0., 1.)]
    manager = make_gaussian_ensemble(specs, [1.])
    calls = []
    def model(R, P, storage, traj):
        calls.append((traj.id, storage.timestep))
        return dict(H_dia=(.5*R[:, 0]**2+storage.timestep)[:, None, None]+0j,
                    dH_dia=R[:, :, None, None]+0j)
    engine = QTAGEngine(manager, model, dict(dt=.01))
    engine.step()
    assert calls == [(0, 0), (1, 0), (0, 1), (1, 1)]
    manager.storage.iM[0, 1] = 2.
    with pytest.raises(ValueError, match='common positive inverse masses'):
        engine.step()
    assert manager.storage.timestep == 1
    with pytest.raises(ValueError, match='dt'):
        QTAGEngine(manager, model, {'dt': float('nan')})


def test_run_reuses_pair_matrices_but_step_rebuilds():
    """One model evaluation per saved time; external edits remain visible to step."""
    spec, mass, model, pars = settings()
    manager = make_gaussian_ensemble(spec, mass)
    calls = []
    def counted(**kwargs):
        calls.append(kwargs['storage'].timestep)
        return model(**kwargs)
    engine = QTAGEngine(manager, counted, pars)
    engine.run(3)
    assert calls == [0, 1, 2, 3]
    manager.storage.q[0, 0, 0] += .01
    engine.step()
    assert calls[-2:] == [3, 4]


def test_common_engine_clock_recording_and_cleanup():
    """Common initialize/run/step share one clock and obey saver cadence."""
    from libra_py.dyn.engine import DynamicsEngine, DynamicsRunResult
    spec, mass, model, pars = settings()
    pars.update(nprint=2, nsteps=3)
    manager = make_gaussian_ensemble(spec, mass)
    class Saver:
        def __init__(self): self.records, self.closed = [], False
        def save_step(self, step, data, metadata=None):
            self.records.append((step, float(data['time']), metadata['method']))
        def close(self): self.closed = True
    saver = Saver()
    with DynamicsEngine(manager, model, method='qtag', params=pars, saver=saver) as engine:
        assert isinstance(engine, DynamicsEngine)
        out = engine.run()
        assert isinstance(out, DynamicsRunResult)
        assert [s[0] for s in saver.records] == [0, 2]
        assert engine.time == pytest.approx(.015)
        engine.step(.02)
        assert engine.time == pytest.approx(.035)
        assert manager.storage.timestep == 4
        assert saver.records[-1] == (4, engine.time, 'qtag')
        restarted = engine.run(0)
        assert float(restarted.histories[0]['time'][0]) == pytest.approx(.035)
    assert saver.closed


def test_common_engine_initialize_and_mapping_savers():
    """All active trajectories save once per time; context closes mapped savers."""
    from libra_py.dyn.engine import DynamicsEngine
    spec, mass, model, pars = settings()
    manager = make_gaussian_ensemble([spec, spec], mass)
    class Saver:
        def __init__(self): self.steps, self.closed = [], False
        def save_step(self, step, data, metadata=None): self.steps.append(step)
        def close(self): self.closed = True
    savers = {0: Saver(), 1: Saver()}
    with DynamicsEngine(manager=manager, model_fn=model, method='qtag', params=pars, saver=savers) as engine:
        engine.initialize()
        engine.step()
        assert all(s.steps == [0, 1] for s in savers.values())
        assert engine.storage.timestep == 1
    assert all(s.closed for s in savers.values())
