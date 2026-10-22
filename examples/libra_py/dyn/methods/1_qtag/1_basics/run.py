# Copyright (C) 2026 Alexey V. Akimov
# SPDX-License-Identifier: GPL-3.0-or-later
"""Storage-backed QTAG tutorial; save individual PNGs and generic NPZ snapshots.

PYTHONPATH=src python examples/libra_py/dyn/methods/1_qtag/1_basics/run.py --backend numpy
Use --backend torch for tensors, --case ferretti for 2D, --steps 10 for a smoke run.
"""
import argparse
from contextlib import nullcontext
from pathlib import Path
from libra_py.dyn.initialization.ensembles import make_gaussian_ensemble
from libra_py.dyn.engine import DynamicsEngine
from libra_py.dyn.models import Holstein5Model, FerrettiModel
from libra_py.dyn.observables.gaussian import reconstruct_density
from libra_py.dyn.observables import plotting
from libra_py.dyn.savers import FaultTolerantSaver, HDF5TimeSeriesSaver
from libra_py.dyn.backends.arrays import to_numpy


def configuration(case='holstein', backend='numpy', device='cpu'):
    """Return explicit initial conditions, masses, model and QTAG parameters."""
    if case == 'holstein':
        # Two harmonic wells and a Gaussian electronic coupling; atomic units.
        model = Holstein5Model(backend=backend, params={
            'E_n': [0., 15.8114], 'x_n': [0., 1.], 'k_n': [10., 10.],
            'V': [[0., 1.], [1., 0.]],
            'alpha': [[0., 1.5811], [1.5811, 0.]],
            'x_nm': [[0., 2.], [2., 0.]],
            'ex_cpl_A': 1., 'ex_cpl_B': 1.5811, 'ex_cpl_C': 2.,
        })
        initial = dict(wfc_q0=[-2.], wfc_p0=[0.], wfc_a0=[1.], wfc_s0=[0.],
                       grid_dims=[25], alp_scl=[6.])
        inverse_masses = [1.]
        dt, nsteps, approximation = .01, 800, 3  # BAT plus exact coupling
        grid = dict(xmin=[-4.], xmax=[4.], npoints=[100])
    elif case == 'ferretti':
        model = FerrettiModel(backend=backend, params={
            'X1': 4., 'X2': 3., 'X3': 3., 'Kx': .02, 'Ky': .01,
            'Delta': .01, 'alpha': 3., 'beta': 1.5, 'gamma': .005,
        })
        initial = dict(wfc_q0=[0., 2.], wfc_p0=[0., 0.], wfc_a0=[1., 1.],
                       wfc_s0=[0., 0.], grid_dims=[5, 5], alp_scl=[6., 6.])
        inverse_masses = [1/20000., 1/6667.]
        dt, nsteps, approximation = 1., 100, 0  # BAT
        grid = dict(xmin=[-2., -4.], xmax=[6., 4.], npoints=[25, 25])
    else:
        raise ValueError('case must be holstein or ferretti')
    # wfc_a0 is the initial GWP width. Basis alpha=wfc_a0*alp_scl multiplies
    # -(r-q)^2/2, so a value 6 gives exp[-3*(r-q)^2]. rho_cut sets grid bounds.
    initial.update(states=[0, 1], initial_state=0, init_placement=0,
                   rho_cut=1e-12, backend=backend, device=device)
    dynamics = dict(dt=dt, nsteps=nsteps, states=[0, 1], target_states=[0, 1],
                    qtag_pot_approx_method=approximation,
                    q_update_method=1, p_update_method=1, a_update_method=1,
                    s_update_method=0, q_sync_method=1, p_sync_method=0,
                    a_sync_method=0, s_sync_method=0, decpl_den=.01,
                    mom_calc_type=1, momentum_fit='legacy', linfit_beta=1e-3,
                    d_weight=1, overlap_rcond=1e-12)
    return initial, inverse_masses, model, dynamics, grid


def main(backend='numpy', case='holstein', steps=None, output=None, device='cpu', save_format='npz'):
    """Initialize TBFs, evolve shared tensors, save and plot coherent observables."""
    initial, mass, model, dynamics, grid = configuration(case, backend, device)
    if steps is not None:
        dynamics['nsteps'] = steps
    output = Path(output) if output else Path(__file__).resolve().parent / 'output' / f'{case}_{backend}'
    context = nullcontext()
    if backend == 'torch':
        import torch
        context = torch.no_grad()
    with context:
        # One Trajectory owns the whole coherent wavefunction; each Gaussian
        # on each surface is a registered TBF, not an independent trajectory.
        manager = make_gaussian_ensemble(initial, mass)
        if save_format == 'hdf5':
            saver = HDF5TimeSeriesSaver(output/'trajectory.hdf', mode='w', flush_stride=100)
        elif save_format == 'npz':
            saver = FaultTolerantSaver(output/'snapshots', mode='w', compressed=True)
        else:
            raise ValueError('save_format must be npz or hdf5')
        try:
            result = DynamicsEngine(manager=manager, model_fn=model, method='qtag',
                                    params=dynamics, saver=saver).run()
        finally:
            if hasattr(saver, 'close'):
                saver.close()
        history = result.histories[0]
        plotting.energy_and_pops(history, output)
        plotting.trajectories(history, output)
        n = dynamics['nsteps']
        # Source tutorial files omit t=0. Select equal physical times, not
        # equal file row numbers: 1,200,400,600,800 in the 800-step history.
        snaps = sorted({1, *[max(1, int(n*f)) for f in (.25, .5, .75, 1)]}) if n else [0]
        density = reconstruct_density(history, **grid, snaps=snaps)
        plotting.wf_plot(density, output)
    print(f'Method: qtag; backend: {backend}; model: {case}; completed steps: {dynamics["nsteps"]}')
    print('Final populations:', to_numpy(history['pops'][-1]))
    print('Final norm:', float(history['norm'][-1]))
    print('Energy change (Ha):', float(history['dEtot'][-1]))
    print('Saved:', output)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--backend', choices=['numpy', 'torch'], default='numpy')
    parser.add_argument('--case', choices=['holstein', 'ferretti'], default='holstein')
    parser.add_argument('--steps', type=int)
    parser.add_argument('--output')
    parser.add_argument('--device', default='cpu')
    parser.add_argument('--save-format', choices=['npz', 'hdf5'], default='npz')
    main(**vars(parser.parse_args()))
