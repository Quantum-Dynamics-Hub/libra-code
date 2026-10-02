# Manuscript LDR and DVR examples

Each model has its own `settings.py`, `run.py`, and `plot.py`. From that model's
folder, run:

```bash
python run.py
python plot.py
```

No environment variables, installation of the source checkout, or command-line
parameters are required. Python must have PyTorch, NumPy, and Matplotlib
installed. The compiled Libra extension is not required. Scripts also work when
invoked by absolute path from another directory. They locate the repository's
Python sources themselves and always write beside their own `settings.py`.

| Folder | Calculation | Output directory |
| --- | --- | --- |
| `tully1` | Tully simple avoided crossing, 1300 a.u. | `tully1/output/` |
| `tully2` | Tully dual avoided crossing, 1500 a.u. | `tully2/output/` |
| `flv_2d` | FLV weak and strong coupling, 3000 a.u. each | `flv_2d/output/weak/` and `flv_2d/output/strong/` |

`run.py` runs both methods sequentially. `plot.py` creates PNG and PDF figures:

- `comparison`: LDR solid curves and DVR open circles for populations (plus
  phase-aligned coherences in 1D), kinetic/potential/total energies, and mean
  positions. FLV has weak and strong coupling in adjacent columns, like the
  manuscript's Figure 4; its combined figure is in `flv_2d/output/`.
- `density_snapshots`: total nuclear probability densities at selected times.
  FLV shows LDR and DVR in separate rows at 0, 1000, and 2000 a.u., with a
  common color scale for each time, in each coupling's output folder.

Each method saves a compact `.pt` file containing observables and selected
densities, plus a JSON copy of its settings. Large LDR overlap, Hamiltonian,
and propagator matrices are not written to disk. Output directories are created
automatically and are ignored by Git. Repeating a run replaces that method's
results in its output directory. Plotting rejects mismatched LDR/DVR settings.

## Parameters and physical conventions

All quantities are in atomic units. Edit the local `settings.py` to change a
model's parameters or `OUTPUT_DIR`. Both methods use dt=1 and record observables
every 10 steps, including the initial and final times.

Tully 1 and 2 use 841 regular centers on [-20, 22], spacing 0.05, Gaussian
exponent 200, mass 2000, q0=-6, p0=30, and probability standard deviation
1/sqrt(2). The initial electronic reference is the ground state at q0.

FLV uses the manuscript's Kx=0.02, Ky=0.10, X1=4, X2=X3=3, Delta=0.01,
a=3, b=1.5, and coupling gamma=0.01 or 0.08. Masses are (20000, 6667),
q0=(2, 0), p0=(0, 0), and probability standard deviations are (0.150, 0.197).
The initial electronic reference is the upper state at q0. The LDR grid has
7381 centers on [1, 7] x [-1.5, 1.5], with spacing (0.05, 0.05) and exponents
(200, 200). The DVR box is [0, 10] x [-2, 2], spacing (0.05, 0.05), matching
the larger reference box in the archived FLV scripts.

The full FLV LDR calculation is a large dense calculation: its compound matrix
dimension is 14762, and one complex128 matrix alone occupies about 3.49 GB.
Several such matrices and diagonalization workspaces are held simultaneously.
Run the full 2D defaults on a machine with sufficient memory for tens of GB of
working storage. The full grids are not silently reduced by the scripts.

The shared implementation in `_common.py` calls `ldr_torch.ldr_solver` for
matrix construction, corrected initial projection, and propagation. The DVR
reference uses `exact_torch.exact_tdse_solver_multistate` potential/kinetic
operators and electronic transformations with a small example adapter for
float64/complex128 grids, `torch.fft.fftfreq` ordering (including odd grids),
and recording only the required observables/snapshots. DVR is periodic SOFT:
half potential step, full Fourier kinetic step, half potential step.
The listed endpoints are sampled; the FFT period is `number_of_points * spacing`.

Both methods start from the **same fixed-reference product state**
`chi0(q) * phi_istate(q0)`, as defined in the manuscript's Eq. 22. LDR supplies
overlaps with all electronic states and solves `S C0 = b`; DVR evaluates this
same product on its grid. This corrects the archived scripts' discrepancy
between fixed-reference and coordinate-dependent adiabatic initial states.
Consequently these are corrected comparisons, not a promise to reproduce the
older plotted curves exactly.

LDR uses the symmetrized endpoint potential and Lowdin density-matrix estimator;
DVR populations/coherences are integrated in the adiabatic representation.
Their finite-grid difference is part of the comparison. Adjacent eigenvector
phases are aligned consistently in 1D to compare coherences. The 2D figures
show populations, not gauge-dependent coherences. Nuclear density snapshots
are reconstructed from the full LDR wavefunction, including electronic overlaps.

These examples cover the manuscript's regular-grid model benchmarks. They do
not perform a quasi-regular/thermal-grid convergence study.

## Optional execution checks and separate runs

In any of the three model folders:

```bash
python run.py --quick
python plot.py --quick
```

Quick mode runs 20 steps on smaller grids and writes to **`output/quick/`**,
leaving the default results separate. It tests execution, not convergence or
scattering accuracy. Tully quick runs retain spacing 0.05 to resolve momentum
30, but use a smaller [-12, 4] box; FLV quick runs use LDR spacing 0.3 and DVR
spacing 0.1. Full and quick outputs both record their actual settings.

To rerun just one method or select a device:

```bash
python run.py --method dvr
python run.py --method ldr --device cuda
```

The default is CPU with four PyTorch threads; `--threads` can change that.
