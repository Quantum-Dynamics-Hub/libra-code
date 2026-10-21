# Gaussian wavepacket integral tests

`test_integrals.py` tests `libra_py.dyn.gwp.integrals`. It combines agreement
with the Python-exposed C++ implementation, independent numerical integration,
physical identities, and differentiation checks. C++ agreement alone would
not detect a formula error shared by both implementations.

## Running the suite

From the repository root, with NumPy and pytest installed:

```bash
PYTHONPATH=src:_build/src python -m pytest -q -rs unittests/test_libra_py/test_dyn/test_gwp/test_integrals.py
```

For C++ comparisons, first build the corrected extension:

```bash
cmake --build _build --target gwp libra_core -j 4
```

Use your build directory in `PYTHONPATH` if it differs from `_build`. The tests
import `liblibra_core` on Linux and `cyglibra_core` on Cygwin. An older extension
may still import successfully but fail the regression comparisons, including
the scalar tuple-returning product-decomposition overload.

PyTorch tests require PyTorch; these tests do not require SymPy.
NumPy tests run without PyTorch. Missing C++/PyTorch modules cause their tests
to skip; `-rs` shows the reasons.

To run only the independent numerical integration checks:

```bash
PYTHONPATH=src python -m pytest -v unittests/test_libra_py/test_dyn/test_gwp/test_integrals.py -k quadrature
```

You can also run directly from this directory, with the package installed or
its source directory on your Python path:

```bash
python -m pytest -q -rs test_integrals.py
```

For an uninstalled source checkout, use paths relative to this directory:

```bash
PYTHONPATH=../../../../src:../../../../_build/src python -m pytest -q -rs test_integrals.py
```

## What each group establishes

| Test | Reference and purpose |
| --- | --- |
| `test_cpp_integrals` | Four operators, both backends, scalar/vector/column/legacy overloads. Multiple DOFs and hbar != 1 exercise the corrected legacy kinetic factors. |
| `test_cpp_mass_and_matrices` | C++ rectangular matrices (two versus four packets), electronic-state masks, and per-DOF inverse masses. |
| `test_cpp_values` | C++ density, wavefunction, legacy values and available spatial gradients. Exercises the corrected exponent initialization. |
| `test_cpp_product_and_triple` | All vector/scalar decomposition outputs and triple integrals, including the pi normalization factor and scalar Python binding. |
| `test_quadrature` | Direct numerical integrals of the defining wavefunctions and derivatives, pointwise product identity, and conjugation of both first packets in the triple integral. Independent of C++ and completed-square formulas. |
| `test_symmetries_and_diagonal` | Unit self-overlap, expected diagonal position/coupling/kinetic values, Hermitian and anti-Hermitian exchange identities, and conversion between general and legacy conventions. |
| `test_autograd` | PyTorch Jacobians versus finite differences for all pair parameters and product-decomposition outputs. |
| `test_matrix_autograd` | Gradients through rectangular pairwise broadcasting, including inverse-mass gradients for kinetic matrices. |
| `test_value_gradient` | Spatial derivatives versus centered finite differences of wavefunction and density values. |
| `test_invalid_inputs` | Python exceptions for invalid widths, dimensions, inverse masses, zero hbar, and matrix input shapes. Invalid C++ inputs are not called because some terminate the process. |
| `test_torch_dtype_device` | float32/float64 inputs retain device and produce complex64/complex128 outputs. CPU always runs when PyTorch is present; CUDA runs only when available. |

## Numerical choices and scope

The `packets` helper uses fixed random seeds and positive widths in `[0.3, 2.0)`.
Inputs are in atomic units unless a legacy overload explicitly supplies hbar.
Most comparisons use float64/complex128 with relative tolerance `2e-12` and
absolute tolerance `2e-13`; the absolute tolerance matters near zero.

Quadrature uses 400 Gauss-Legendre points on a finite interval extending
`12/sqrt(min(alp1, alp2))` beyond the outermost centers. For these parameters
the omitted Gaussian tails are negligible. This is a one-dimensional reference;
multidimensional behavior is covered by C++ comparisons and algebraic identities.
It is not an accuracy guarantee for arbitrary extreme input values.

Spatial finite differences use a step of `1e-5` and looser tolerances (`1e-8`
relative, `1e-10` absolute) for truncation and roundoff error. PyTorch `gradcheck`
has been replaced by the local `check_jacobian` helper, which uses double
precision, a centered-difference step of `1e-6`, relative tolerance `1e-3`,
and absolute tolerance `1e-5`. It compares the full Jacobian for all real and
imaginary output components using scalar autograd calls with implicit unit
seeds. This avoids the symbolic-shape path that imports SymPy in PyTorch
`gradcheck`. `test_jacobian_checker` checks the helper itself, including that
it rejects a deliberately detached term. Test-only
conversion helpers detach tensors for value assertions; gradient checks call
the implementation directly and do not detach its outputs.

With C++ and PyTorch available the suite collects 66 cases. A passing CPU run
does not imply CUDA execution; that depends on `torch.cuda.is_available()`.
