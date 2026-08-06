# Hop-proposal examples

These examples reproduce and extend the numerical cases used by the hopping
unit tests. Run them directly from this directory; repository-root invocation
is intentionally not supported.

```bash
cd examples/libra_py/dyn/hopping
python example_direct_probabilities.py
python example_lz_zn.py
python example_interface_and_sampling.py
python example_acceptance_probabilities.py
python example_acceptance_and_rescaling.py
```

The examples cover:

- FSSH with and without uphill Boltzmann scaling;
- GFSH, original GFSH, FSSH2, FSSH3, MSSH, and MASH;
- diabatic and adiabatic Landau–Zener probabilities;
- multidimensional Zhu–Nakamura probabilities;
- the common `tsh_method` interface, batched trajectories, and stochastic
  state selection.
- quantum, classical, and harmonic-oscillator thermal probabilities;
- energy-, derivative-coupling-, and force-based hop acceptance;
- uniform, derivative-coupling, and frustrated-hop momentum rescaling.

Option 5 (DISH) remains an explicit placeholder in the common interface
because the C++ implementation handles it through a separate decoherence-event
scheduler.

The active Libra environment must already provide `libra_py`. The scripts do
not alter `sys.path` or locate source files themselves.
