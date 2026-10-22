# Tully-1 EHRENFEST with DynamicsEngine

Run from the repository root (or set PYTHONPATH to the absolute src directory):

```bash
PYTHONPATH=src python examples/libra_py/dyn/methods/3_ehrenfest/1_basics/run.py
```

Use `--steps 10 --ntraj 4` for a smoke run, `--save-stride 5` for the
output interval, `--output PATH` for a different destination, and
`--save-format hdf5` for buffered HDF5 (requires h5py). Default output is
`output/tully1_numpy` beside the script. Rerunning replaces saved records.

`configuration()` defines integration/output parameters, `make_system()`
constructs a seeded incoming Gaussian ensemble, and `main()` runs
`DynamicsEngine(method="ehrenfest")`, saves snapshots, and makes plots.
Both Tully examples use identical initial conditions, 12 trajectories,
1200 steps of 0.5 a.u., and four electronic substeps.

Outputs: `snapshots/step_*.npz` plus its manifest (or `trajectory.hdf`),
`energy.png`, `energy_change.png`, `norm.png`, `populations.png`,
`position_dof0.png`, and `momentum_dof0.png`. FSSH additionally plots
active states and surface-hopping populations. All PNGs are headless.

Electronic populations are adiabatic ensemble averages; norm is their sum.
Energy change is mean total energy minus its initial value, whereas the saved
`dEtot_ave` is the ensemble energy spread. Nuclear plots show every trajectory.
Ehrenfest active labels are bookkeeping, not physical hopping populations.
QTAG instead reports coherent diabatic populations and the metric norm.
