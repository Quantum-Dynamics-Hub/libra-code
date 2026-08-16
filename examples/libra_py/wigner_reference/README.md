# Water Wigner-sampling reference

This directory contains a compact water-molecule fixture for the Wigner
sampling example and unit tests.

- `water_equilibrium.xyz`: equilibrium geometry in angstrom.
- `mode_1.xyz` through `mode_3.xyz`: mass-orthogonal Cartesian displacement
  vectors for the bend, symmetric stretch, and asymmetric stretch.
- `water_hessian_au.txt`: a 9 x 9 Cartesian Hessian in hartree/bohr^2.

The geometry and frequencies (1594.75, 3657.05, and 3755.93 cm^-1) are
representative gas-phase water values. The displacement vectors are idealized
water symmetry coordinates. The Hessian was reconstructed from those
mass-orthonormal modes using O and H masses of 15.999 and 1.00784 Da. It is a
realistic harmonic reference fixture, not the verbatim output of a particular
electronic-structure program.

