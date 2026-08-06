# Decoherence

This package is the NumPy translation of three C++ sources:

- `dyn_decoherence_time.cpp` -> `times.py`
- `dyn_decoherence_methods.cpp` -> `methods.py`
- `dyn_methods_dish.cpp` -> `dish.py`

All Python arrays are trajectory-major. Amplitudes have shape
`(ntraj, nstates)`, rates `(ntraj, nstates, nstates)`, and state forces
`(ntraj, nstates, ndof)`. A single-trajectory vector or matrix is accepted by
the elementary functions where it is unambiguous.

The `DynamicsEngine` uses the original integer choices for
`decoherence_times_type`: `-1` off, `0` supplied rates, `1` EDC, `2`
Schwartz-1, `3` Schwartz-2, `4` Schwartz-1 interaction-width, and `5`
Gu--Franco. It likewise preserves `decoherence_algo=0` for SDM, `1` for
instantaneous decoherence, and `7` for revised DISH. Old DISH remains
`tsh_method=5` and, consistently with C++, must be paired with
`decoherence_algo=-1`.

AFSSH, BCSH, and the independent-trajectory exact-factorization primitives
need auxiliary nuclear variables. Their reusable mathematical operations live
in `methods.py`; engine choices requiring orchestration not yet represented by
the Python storage workflow fail explicitly rather than silently doing
nothing.
