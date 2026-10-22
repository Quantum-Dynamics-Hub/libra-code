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
spectral overlap cutoff replaces explicit matrix inversions. The momentum fit
uses `momentum_fit="legacy", linfit_beta=1e-3`, reproducing the notebook
solver with NumPy/PyTorch operations. Initial projected coefficients are normalized.
The earlier direct-fit implementation changed the basis motion and did not
reproduce the Holstein snapshots; using a more accurate fit is a different
numerical experiment, not a drop-in replacement for this tutorial.

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
The archive includes the actual initial state at index 0. For Holstein, plots
select indices `[1,200,400,600,800]`, at times `[0.01,2,4,6,8]` a.u. These match
the notebook file rows `[0,199,399,599,799]`, which omit the initial state.
Shorter runs use the first completed step and quarter-duration checkpoints
(duplicates removed); a zero-step run plots the initial state.
Explicit output paths are overwritten on a repeated run.

## Numerical validation and interpretation

The corrected Holstein case was run on NumPy and CPU PyTorch; Ferretti was
run on NumPy, all in float64 without SymPy. For the parameters above:

| Case | Final populations | Final norm | Energy change (Ha) |
| --- | --- | --- | --- |
| Holstein, 800 steps | about `[0.738272, 0.256036]` | about `0.994307` | about `-0.538847` |
| Ferretti, 100 steps | about `[0.99999948, 4.72e-11]` | about `0.99999948` | about `-2.08e-6` |

The five Holstein density snapshots were compared point by point on the
100-point spatial grid with a fresh run of the original Python/C++ tutorial
workflow. Maximum absolute density errors were below `1.5e-8` for both NumPy
and CPU PyTorch. The unit suite includes an optional full-run regression
against these independently generated reference densities. These are floating
point comparisons, not a promise of bitwise equality across BLAS/devices.

The finite moving basis still loses some norm and energy. These quantities are
reported without renormalization; reproducing the tutorial does not establish
convergence with respect to timestep or basis size.

Useful follow-up exercises mirror the notebook: vary the duration, number of
packets, approximation, frozen/moving widths, and synchronization threshold.
Use the individual norm/energy plots to assess each change.
