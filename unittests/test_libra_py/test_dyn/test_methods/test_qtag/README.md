# QTAG tests

From the repository root:

```bash
cmake --build _build --target libra_core -j 4
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=src:_build/src python -m pytest -q -rs unittests/test_libra_py/test_dyn/test_methods/test_qtag
```

Or run `python -m pytest -q -rs` from this directory. `conftest.py` adds the
repository source directory. C++ parity checks require the rebuilt extension
on the Python path (e.g. `PYTHONPATH=../../../../../_build/src` from here).
They skip when the extension is not available. Torch cases skip if Torch is
absent; plotting/HDF5 checks require Matplotlib/h5py respectively. SymPy is
not required, including for the scalar-autograd derivative checks.

`test_core.py` covers:

- All exposed C++ QTAG operations, including all four potential methods and
  generalized coefficient propagation from the electronic module.
- Unequal widths, multiple spatial dimensions, complex coefficients and phases.
- Variable endpoint Hessians and mixed-coordinate derivatives to exercise the
  C++ LHA corrections; permutation of state labels to test assembly reuse.
- Independent 2D numerical quadrature of local quadratic Taylor polynomials
  and exact Gaussian couplings, rather than relying on cross-language agreement.
- Log-gradient and second-derivative finite differences.
- Norm conservation/time reversal in a nonorthogonal fixed basis.
- NumPy/PyTorch parameter derivatives compared with scalar finite differences.

`test_workflow.py` covers:

- Normalized projection and reproducible initialization, including nonzero phases.
- Fixed-basis energy/norm conservation and correct initial/final snapshot times.
- Agreement of complete moving-basis NumPy/PyTorch histories and no input mutation.
- Linear momentum fitting, zero-density nodes and metric rank deficiency.
- Short end-to-end Torch gradients (nondegenerate metric, no explicit gradient
  seeds that would import SymPy).
- NPZ/HDF5 round trips, wavefunction reconstruction and individual headless PNGs.
- A fresh process with C++ and SymPy imports blocked, proving the workflow's
  numerical NumPy path is independent of both.

Direct formula comparisons use relative/absolute tolerances around `1e-10` /
`1e-11`; numerical derivatives use looser tolerances for finite-difference
error. End-to-end parity tests use short trajectories to avoid treating
accumulated differences of ill-conditioned long adaptive runs as exact equality.
Full 800-step Holstein and 100-step Ferretti CPU runs are documented in the
example README; the full Holstein regression is opt-in as described below.

The momentum-fit tests compare NumPy and Torch against C++ `solve_linsys`
at three stopping tolerances. They check both momentum and real logarithmic
gradient, including their slopes, and ensure that an exact fit is distinguishable.

To additionally run the full Holstein example regression on both backends:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 QTAG_FULL_TESTS=1 python -m pytest -q test_tutorial.py
```

This takes longer, especially for Torch's iterative fit. It checks the actual
snapshot selection and both state densities at all five tutorial times against
the original Python/C++ workflow. See `data/README.md` for reference provenance.
All outputs go to pytest's temporary directory; no source tutorial is required.


# Algorithm tests

`test_python_npi.py` tests NPI against its C++ implementation.
`test_engine_qtag.py` tests the integrated TBF/storage QTAG driver against the test-only
`qtag_reference` implementation, without requiring C++ or SymPy.

The QTAG tests now live in `../test_methods/test_qtag`. Run from that directory:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q test_engine_qtag.py
```

Short tests compare all saved energies, populations, coefficients, coordinates,
momenta, widths, phases and reconstructed densities on NumPy and Torch. They
cover both 1D Holstein and 2D Ferretti, BAT/LHA/exact-coupling variants,
noncontiguous slots/labels, independent trajectories, inactive trajectories,
pruning/compaction, timestep alignment, input validation, and a scalar Torch
gradient checked against finite differences. The gradient test uses a
nondegenerate metric; repeated eigenvalues need not have finite eigenvector
derivatives. A subprocess blocks the reference package, SymPy and C++ imports.

Full 800-step Holstein and 100-step Ferretti comparisons are opt-in because
the iterative Torch momentum fit is slower:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 QTAG_FULL_TESTS=1 \
  python -m pytest -q test_engine_qtag.py
```

These four tests compare entire histories and all five density snapshots with
the preserved method implementation, on both backends. Reference densities
are computed at runtime, so no new binary trajectory fixture is needed.

Tests for the shared modules live beside their corresponding source packages:

- `test_core/test_gaussian_storage.py`: allocation, spawning, expansion, views,
  parent copying and compaction of optional Gaussian fields.
- `test_initialization/test_gaussian_initialization.py`: deterministic grids,
  seeded sampling, phases, noncontiguous labels and projection coefficients.
- `test_gwp/test_packets.py`: unequal-width Gaussian kernels vs the reference.
- `test_models/test_gaussian_hessians.py`: finite differences of gradients,
  including mixed derivatives and the Ferretti Y=0 point.
- `test_observables/test_gaussian_observables.py`: coherent cross terms, generic
  saver round trips, density reconstruction, PNGs and disk-only histories.

All tests use float64 except explicit allocation precision checks. Torch tests
skip if unavailable. C++ is only needed for the pre-existing NPI tests.

`test_core.py` and `test_workflow.py` exercise the frozen test-only oracle;
`test_engine_qtag.py` compares the production common engine with that oracle.
`test_tutorial.py` exercises the migrated production example against the C++
density fixture. Production code never imports `qtag_reference`.
