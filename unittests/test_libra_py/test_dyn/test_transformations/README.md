# Transformation unit tests

Run these tests directly from this directory using an environment in which
`libra_py` is installed:

```bash
cd unittests/test_libra_py/test_dyn/test_transformations
pytest
```

The test files do not modify `sys.path` and do not reconstruct or access the
repository root. `test_three_state_tracking.py` covers
cyclic three-state identity changes, complex phases, deterministic and
stochastic tracking, batched projectors, force costs, LD/SVD updates, and
phase-corrected state mappings.

`test_reordering_regression.py` adapts the legacy `test_reordering.py` data:
eight difficult mixed-state overlap matrices and every exact permutation for
two, three, and four states. The legacy test data are copied into this local
suite, so the tests do not import or access files outside this directory.
