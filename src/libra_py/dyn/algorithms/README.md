# Numerical algorithms

NPI lives in `npi.py`. Dynamics method implementations live in `../methods`;
see [QTAG](../methods/qtag/README.md) for the common-engine specialization.

`quantum_potential.py` owns Gaussian quantum-momentum/logarithmic-gradient
estimation and direct/legacy linear fits. It uses shared Gaussian packet kernels
and preserves NumPy/PyTorch types. It does not advance the basis or project
coefficients; those operations remain in `propagation/gaussian.py`.
