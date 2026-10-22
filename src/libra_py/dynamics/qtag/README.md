# QTAG implementations

This directory contains the original QTAG workflow based on C++ `MATRIX`,
`CMATRIX`, `nHamiltonian`, and legacy text/HDF5 savers.

For the NumPy/PyTorch implementation, use:

```python
from libra_py.dyn.engine import DynamicsEngine
from libra_py.dyn.initialization import make_gaussian_ensemble
```

The new implementation and migration details are documented in
[`dyn/methods/qtag/README.md`](../../dyn/methods/qtag/README.md). It uses
array/tensor model callbacks with `DynamicsEngine(method="qtag")`, returns
shared storage and per-trajectory histories, and saves
plots to individual images. The legacy import path and callable signatures in
this directory remain unchanged for existing C++ workflows.
