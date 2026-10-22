# QTAG with shared TBF/tensor infrastructure

This example reproduces the Holstein and Ferretti tutorial calculations using
`dyn.initialization`, `dyn.models`, `dyn.engine.DynamicsEngine`, `dyn.observables`,
and `dyn.savers`. The engine delegates numerical steps to `dyn.methods.qtag`. The old implementation
is retained only in the test suite for numerical comparisons.

Run from the repository root:

```bash
PYTHONPATH=src python examples/libra_py/dyn/methods/1_qtag/1_basics/run.py --backend numpy
PYTHONPATH=src python examples/libra_py/dyn/methods/1_qtag/1_basics/run.py --backend torch
PYTHONPATH=src python examples/libra_py/dyn/methods/1_qtag/1_basics/run.py --case ferretti
```

`--steps 10` selects a short run. `--backend torch --device cuda` selects GPU
storage when available. `--output PATH` controls output (default: `output/<case>_<backend>` beside the script); rerunning replaces
that example's snapshot series and plots. `OPENBLAS_NUM_THREADS=1` and
`OMP_NUM_THREADS=1` can reduce overhead for the small dense matrices.

`configuration()` shows all initial wavepacket, mass, potential, update,
synchronization and fit parameters explicitly. `main()` creates a manager,
registers the basis as TBFs, runs QTAG, saves generic NPZ records and generates
individual PNGs. There is no separate example-local model implementation:
analytic derivatives/Hessians come from the shared model classes.

A single `Trajectory` owns a whole coherent Gaussian expansion. Its TBFs are
stored in `(trajectory,slot,dof)` tensors, with coefficients in `ampl_dia`
on each packet's labeled state. `result.histories[0]` provides the first
trajectory's copied analysis history in `(time,slot,dof)` order. To evolve
several independent wavefunctions, pass a list of initial-condition dictionaries
to `make_gaussian_ensemble`; see `src/libra_py/dyn/methods/qtag/README.md`.

The default Holstein run uses 25 packets per surface and 800 steps of 0.01 a.u.;
Ferretti uses 5x5 packets per surface and 100 steps of 1 a.u. Both retain the
source tutorial's iterative momentum fit with tolerance 0.001. Holstein density
snapshots use times 0.01, 2, 4, 6 and 8 a.u. The saved history includes t=0,
so its indices are 1, 200, 400, 600 and 800 (the tutorial's text row indices
are one smaller).

Outputs include:

- `snapshots/step_*.npz` and `manifest.jsonl`, using `FaultTolerantSaver`.
- Energy, energy change, norm, and state-population PNGs.
- Position, momentum, width and phase PNGs for each DOF.
- Individual 1D density snapshots, or a 2D image per snapshot/state.

Reload the snapshots with existing saver tools:

```python
from libra_py.dyn.savers import load_saved_steps, stack_saved_observables
from libra_py.dyn.observables.gaussian import reconstruct_density

history = stack_saved_observables(load_saved_steps("examples/libra_py/dyn/methods/1_qtag/1_basics/output/holstein_numpy/snapshots"))
density = reconstruct_density(history, xmin=[-4.], xmax=[4.], npoints=[100],
                              snaps=[1, 200, 400, 600, 800])
```

The full Holstein result has populations approximately `[0.738272,0.256036]`,
norm `0.994307`, and energy change `-0.538847` Ha. Full trajectories and density
snapshots on both backends are tested against the unchanged reference method.
The finite-basis norm/energy losses are reported without renormalization.

For faster output, add `--save-format hdf5` (requires h5py). This uses the
existing buffered `HDF5TimeSeriesSaver`, writing `trajectory.hdf` with a flush
every 100 records and on close. Load it with `load_hdf5_timeseries`. The default
NPZ mode retains per-step compression and durable filesystem synchronization;
that stronger per-record durability can dominate runtime on slow filesystems.
