# Decoherence functions and times

Run this example from this directory:

```bash
cd examples/libra_py/dyn/decoherence
python example_decoherence.py
```

It demonstrates EDC, the dephasing-informed correction, Schwartz-1,
Schwartz-2, Gu--Franco, conversion of rates to DISH coherence intervals, SDM,
projection, collapse, and deterministic DISH event detection. It prints the
small numerical arrays directly and creates no output directories.

The active Libra environment must provide `libra_py`; the script deliberately
does not modify `sys.path`.
