# Hop-proposal examples

These examples reproduce and extend the numerical cases used by the hopping
unit tests. Run them directly from this directory; repository-root invocation
is intentionally not supported.

```bash
cd examples/libra_py/dyn/hopping
python example_direct_probabilities.py
python example_lz_zn.py
python example_interface_and_sampling.py
```

The examples cover:

- FSSH with and without uphill Boltzmann scaling;
- GFSH, original GFSH, FSSH2, FSSH3, MSSH, and MASH;
- diabatic and adiabatic Landau–Zener probabilities;
- multidimensional Zhu–Nakamura probabilities;
- the common `tsh_method` interface, batched trajectories, and stochastic
  state selection.

Option 5 (DISH) remains an explicit placeholder in the common interface
because the C++ implementation handles it through a separate decoherence-event
scheduler.
