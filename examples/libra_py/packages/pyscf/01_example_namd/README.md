# PySCF nonadiabatic molecular dynamics

This example connects the PySCF CASSCF strategy to Libra through
`libra_py.packages.pyscf.methods.pyscf_compute_adi`.

## Mixed doublet/quartet example

`run.py` follows the structure and local FSSH2 recipe of the corresponding
OpenMolcas example. It runs three Al3 trajectories for five steps by default,
using separate spin-constrained PySCF CASSCF calculations for two doublet roots
and one quartet root. Expansion over Ms produces the same eight-state ordering
as the OpenMolcas example.

Run the dynamics first, then convert the stored coordinates to XYZ:

```bash
python run.py
python make_trajectory.py
```

`run.py` writes the simulation data and plots under `FSSH2_/`. In particular,
`make_trajectory.py` requires `FSSH2_/mem_data.hdf`; it reads the `q/data`
dataset, converts coordinates from Bohr to Angstrom, and writes one multi-frame
XYZ file per trajectory under `trajectories/`:

```text
trajectories/traj_0.xyz
trajectories/traj_1.xyz
trajectories/traj_2.xyz
```

The number of files follows `LIBRA_NTRAJ`, and the number of frames in each
file follows `LIBRA_NSTEPS`. Run both commands from this example directory. If
the molecule is changed in `run.py`, update `ATOM_LABELS` in
`make_trajectory.py` to match.

For a faster one-trajectory, one-step smoke test:

```bash
LIBRA_NTRAJ=1 LIBRA_NSTEPS=1 python run.py
python make_trajectory.py
```

The calculation uses time overlaps for nonadiabatic propagation and gradients
for momentum rescaling. Explicit PySCF NAC response vectors are disabled
because spin-constrained state-averaged CASSCF NAC response is not currently
supported. Without spin-orbit coupling, matrix blocks between the doublet and
quartet manifolds are zero, so this example does not model intersystem
crossing. Generated results and plots are written under `FSSH2_/`.

## Earlier adapter example

From this directory, run either configuration:

```bash
python run_namd.py fssh_nacv
python run_namd.py nbra
```

`fssh_nacv` runs non-NBRA FSSH using PySCF derivative-coupling vectors.
`nbra` runs the single-trajectory NBRA configuration using time overlaps.
The example requires an installed Libra Python package and PySCF.
