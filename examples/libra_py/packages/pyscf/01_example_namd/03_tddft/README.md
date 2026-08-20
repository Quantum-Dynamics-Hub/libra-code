# Mixed-spin TDDFT NAMD example

This example applies fewest-switches surface hopping to Al3 with the same
eight-state ordering as `01_casscf`: the UKS/PBE0 reference and first
linear-response TDDFT root of the doublet, each expanded over two Ms
projections, followed by the quartet UKS/PBE0 reference expanded over four
projections. The manifolds are independent because spin-orbit coupling is not
included.

Run the dynamics and export its coordinates from this directory:

```bash
LIBRA_NTRAJ=1 LIBRA_NSTEPS=5 python run.py
python make_trajectory.py
```

The adapter supplies energies, ground- and excited-state analytic gradients,
and factorized determinant overlaps between consecutive unrestricted TDDFT
snapshots. Analytic derivative-coupling vectors are disabled; the dynamics
obtains time-dependent couplings from overlaps and uses gradient differences
for momentum rescaling (`momenta_rescaling_algo = 211`). Output is written
under `FSSH2_/` and XYZ trajectories under `trajectories/`.
