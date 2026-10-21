# QTAG basics: NumPy and PyTorch scripts

This is a pure-Python adaptation of
`Tutorials_Libra/6_dynamics/5_qtag/1_basics/tutorial.ipynb`.
It follows the notebook's model definition, basis placement, initial-state
projection, QTAG propagation, density reconstruction and visualization.
There are no notebook cells or C++ matrix objects.

- `run.py` constructs parameters explicitly and executes the workflow.
- `models.py` supplies the Holstein/Ferretti Hamiltonians and analytic
  derivatives, using `libra_py.dyn.models` for values/first derivatives.

## Run from the repository root

```bash
# Complete 1D Holstein example, 800 steps of dt=0.01 a.u.
PYTHONPATH=src python examples/libra_py/dyn/methods/qtag/1_basics/run.py --backend numpy --output out_numpy
PYTHONPATH=src python examples/libra_py/dyn/methods/qtag/1_basics/run.py --backend torch --output out_torch

# Complete 2D Ferretti example, 100 steps of dt=1 a.u.
PYTHONPATH=src python examples/libra_py/dyn/methods/qtag/1_basics/run.py --case ferretti --backend numpy --output out_2d

# Short smoke run or a different approximation
PYTHONPATH=src python examples/libra_py/dyn/methods/qtag/1_basics/run.py --steps 10 --output out_smoke
PYTHONPATH=src python examples/libra_py/dyn/methods/qtag/1_basics/run.py --method 2 --output out_lha_exact
```

With Libra installed (including an editable install), you can instead run
`python run.py` from this directory. `--backend torch --device cuda` selects a
supported GPU. Matplotlib is needed for plots; Torch is optional for NumPy.
SymPy and the C++ extension are not required. HDF5 is not used by the example.

For small dense matrices, setting `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1`
can avoid excessive BLAS thread overhead. In a restricted environment, set
`MPLCONFIGDIR` to a writable cache directory if Matplotlib cannot use its
normal configuration directory.

## Correspondence to the notebook

The Holstein example has two harmonic diabatic wells, a Gaussian coupling,
and 25 packets per state. The notebook's custom and built-in 1D choices are
the same model, so they are represented by one callback here. Method 3 is BAT
on diagonal elements plus exact Gaussian coupling, as in the notebook.

The Ferretti example has two DOFs with masses 20000 and 6667, two states,
and a 5x5 grid per state. It uses BAT (method 0). `Ky=0.01` follows the
notebook's executable cell; its introductory prose instead says `0.10`.
Exact Gaussian methods 2/3 are rejected for this case because its coupling
is not the sum of Gaussians assumed by those methods. LHA (1) is available.

Both examples use the notebook's q/p/width updates, fixed phases, population
threshold, density-weighted momentum fit, and synchronization choices. A
spectral overlap cutoff and direct linear fit replace fragile explicit matrix
inversions and the iterative fit. Initial projected coefficients are normalized.

The input comments distinguish initial GWP widths (`wfc_a0`) from QTAG basis
widths (`alp = wfc_a0*alp_scl`), explain atomic units and show which parameters
control the basis, potential, and propagation. PyTorch runs use `no_grad()` to
avoid retaining the graph for hundreds of snapshots; the library's core and
short trajectories can also be differentiated.

## Files produced

The output directory contains:

- `trajectory.npz`: complete history, complex coefficients and state labels.
- `observables.txt`: time, unnormalized energy expectation, norm and populations.
- `energy.png`, `energy_change.png`, `norm.png`, `populations.png`.
- Separate `position_dof*.png`, `momentum_dof*.png`, `width_dof*.png`,
  `phase_dof*.png` for each DOF.
- `density_step*.png` snapshots in 1D, or separate snapshots per state in 2D.
- `potential_dof*.png`: diabatic/adiabatic potential cuts for each DOF, with
  other coordinates held at zero.

The default full runs create 14 PNGs in 1D and 24 in 2D. No plot calls `show()`.
Snapshot 0 is the actual initial state; final index is `nsteps`. Five evenly
spaced saved snapshots are selected (duplicates removed for very short runs).
Explicit output paths are overwritten on a repeated run.

## Numerical validation and interpretation

Both full-length cases were run on NumPy and CPU PyTorch in float64 without
SymPy. For the parameters above:

| Case | Final populations | Final norm | Energy change (Ha) |
| --- | --- | --- | --- |
| Holstein, 800 steps | about `[0.83881, 0.15306]` | about `0.991874` | about `-1.5697` |
| Ferretti, 100 steps | about `[0.99999879, 4.61e-11]` | about `0.99999879` | about `-3.32e-6` |

Over the recorded trajectories, the largest NumPy/PyTorch population difference
was about `3.4e-6` in the long Holstein run and `3e-14` in Ferretti. Differences
can accumulate in ill-conditioned moving Gaussian bases. These values are
validation observations, not an exact reference trajectory for every BLAS/device.

The Holstein norm loss and roughly 7.5% drift of its initial `C^H H C` energy
show that the notebook parameter set is not a converged accuracy benchmark.
The code reports these quantities rather than renormalizing or hiding the drift.
Refine the time step/basis and compare approximation choices before interpreting
long-time populations quantitatively. Tests of a fixed basis conserve norm and
energy; moving-basis projection adds a distinct source of numerical error.

Useful follow-up exercises mirror the notebook: vary the duration, number of
packets, approximation, frozen/moving widths, and synchronization threshold.
Use the individual norm/energy plots to assess each change.
