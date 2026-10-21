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
example README; they are integration validation, not part of the unit suite.
