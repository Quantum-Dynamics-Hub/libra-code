# Mixed-spin CISD NAMD example

This example applies fewest-switches surface hopping to Al3 with the same
eight-state ordering as `01_casscf`: two doublet roots expanded over their two
Ms projections, followed by one quartet root expanded over four projections.
The manifolds are independent because spin-orbit coupling is not included.

The PySCF backend uses canonical UHF references and UCISD, which are required
for analytic open-shell CISD gradients. To keep the example practical and
avoid zero occupied/virtual spin spaces in PySCF's UCISD gradient code, the
calculation correlates five frontier orbitals (orbitals 17--21 in the STO-3G
ordering) and freezes all remaining occupied and virtual orbitals.

Run the dynamics and export its coordinates from this directory:

```bash
LIBRA_NTRAJ=1 LIBRA_NSTEPS=5 python run.py
python make_trajectory.py
```

The adapter supplies energies, all state gradients, and consecutive-geometry
UCISD wavefunction overlaps. Analytic derivative-coupling vectors are disabled;
the dynamics obtains time-dependent couplings from overlaps and uses gradient
differences for momentum rescaling (`momenta_rescaling_algo = 211`). Output is
written under `FSSH2_/` and XYZ trajectories under `trajectories/`.
