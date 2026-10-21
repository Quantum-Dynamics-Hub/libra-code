# QTAG implementations

This directory contains the original QTAG workflow based on C++ `MATRIX`,
`CMATRIX`, `nHamiltonian`, and legacy text/HDF5 savers.

For the NumPy/PyTorch implementation, use:

```python
from libra_py.dyn.methods.qtag import core, compute, initialize, save, plot
```

The new implementation and migration details are documented in
[`dyn/methods/qtag/README.md`](../../dyn/methods/qtag/README.md). It uses
array/tensor model callbacks, returns a final-state/history object, and saves
plots to individual images. The legacy import path and callable signatures in
this directory remain unchanged for existing C++ workflows.
