# citools Examples

These examples show how `libra_py.citools` represents many-electron
electronic states and computes overlaps between them.

The examples use the determinant convention used by the Python CI tools:

- `+i` is an alpha electron in spatial orbital `i`
- `-i` is a beta electron in spatial orbital `i`
- orbital indices start from 1
- a closed-shell two-orbital determinant is written as `(1, -1, 2, -2)`

In nonadiabatic dynamics, the central quantity is often a time-overlap
matrix such as `<Phi_I(t) | Phi_J(t+dt)>`.  At the one-electron level this
comes from molecular-orbital overlaps. At the many-electron level the
overlap is a determinant of the occupied-orbital overlap submatrix. For
restricted calculations, alpha and beta spin functions are orthogonal, so
the spin-orbital overlap matrix has zero alpha/beta off-diagonal blocks.

## Examples

1. `01_slater_determinant_overlaps.py`

   Builds determinant overlaps directly. It illustrates why a doubled
   spin-orbital matrix has the block form
   `[[S_alpha, 0], [0, S_beta]]` and compares the factorized alpha/beta
   determinant formula with a full spin-orbital determinant.

2. `02_spin_adapted_csfs.py`

   Constructs singlet and triplet configuration state functions (CSFs)
   from two open-shell determinants. The example shows how the spin-adapted
   singlet is the antisymmetric spin combination while the triplet
   `M_s = 0` component is symmetric.

3. `03_ci_time_overlap.py`

   Computes a CI-state time-overlap matrix from a small model set of
   one-electron time-overlaps and CI coefficients. This mirrors the flow
   used by package interfaces: build Slater-determinant overlaps, transform
   them to spin-adapted CSFs, then contract with CI amplitudes.

Run from the repository root with:

```bash
PYTHONPATH=src python examples/citools/01_slater_determinant_overlaps.py
PYTHONPATH=src python examples/citools/02_spin_adapted_csfs.py
PYTHONPATH=src python examples/citools/03_ci_time_overlap.py
```

