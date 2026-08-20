# PySCF nonadiabatic molecular dynamics examples

This directory demonstrates the same Libra FSSH workflow with three PySCF
electronic-structure strategies:

- `01_casscf`: mixed doublet/quartet state-averaged CASSCF for Al3;
- `02_cisd`: UHF/UCISD doublet and quartet blocks for Al3;
- `03_tddft`: UKS/PBE0 TDDFT doublet and quartet blocks for Al3.

Each subdirectory is self-contained and contains `run.py`, a local
`recipes/fssh2.py`, `make_trajectory.py`, and a method-specific README. Run the
commands from the selected subdirectory:

```bash
python run.py
python make_trajectory.py
```

For a quick end-to-end check, reduce the ensemble and trajectory length:

```bash
LIBRA_NTRAJ=1 LIBRA_NSTEPS=1 python run.py
python make_trajectory.py
```

The dynamics use electronic time overlaps to update the nonadiabatic coupling
matrix. Analytic PySCF derivative-coupling vectors are disabled, and frustrated
hop momentum treatment uses state-gradient differences. Generated HDF5 data,
plots, and XYZ trajectories are ignored by Git.
