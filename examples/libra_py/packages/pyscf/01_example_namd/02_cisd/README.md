# Mixed-spin CISD NAMD example

This example applies fewest-switches surface hopping to Al3 with the same
eight-state ordering as `01_casscf`: two doublet roots expanded over their two
Ms projections, followed by one quartet root expanded over four projections.
The manifolds are independent because spin-orbit coupling is not included.

The PySCF backend uses canonical UHF references and UCISD, which are required
for analytic open-shell CISD gradients. The 18 doubly occupied core orbitals
are frozen, while the full virtual space is retained. Freezing virtual orbitals
must be avoided because PySCF's analytic UCISD gradient is then inconsistent
with the corresponding frozen-space energy.

The active state and its coupled partner belong to the doublet block, for
which both analytic gradients are evaluated. The quartet is uncoupled without
spin-orbit coupling, and its unused gradient is not requested because PySCF's
UCISD gradient cannot handle the resulting empty correlated beta-occupied
space. The second doublet root uses an enlarged Davidson space and cycle limit;
the adapter raises an error rather than returning an unconverged root.

Run the dynamics and export its coordinates from this directory:

```bash
LIBRA_NTRAJ=1 LIBRA_NSTEPS=5 python run.py
python make_trajectory.py
```

The nuclear timestep is 0.1 fs. The adapter supplies energies, the required
doublet gradients, and consecutive-geometry UCISD wavefunction overlaps.
Analytic derivative-coupling vectors are disabled;
the dynamics obtains time-dependent couplings from overlaps and uses gradient
differences for momentum rescaling (`momenta_rescaling_algo = 211`). Output is
written under `FSSH2_/` and XYZ trajectories under `trajectories/`.
