# OpenMolcas nonadiabatic molecular dynamics

This example runs a short FSSH2 nonadiabatic molecular-dynamics calculation
for an Al3 cluster using the Libra OpenMolcas interface. It then converts the
stored nuclear coordinates into one XYZ file per trajectory.

## What the example does

`run.py` initializes three slightly displaced Al3 trajectories on adiabatic
state 2 and propagates them for five steps with a 0.5 fs nuclear time step. The
electronic structure consists of two doublet roots and one quartet root. After
expansion over their spin projections, the Libra Hamiltonian contains eight
states: four doublet components and four quartet components.

The dynamics settings are supplied by the local `recipes/fssh2.py` module.
OpenMolcas calculations are performed in separate directories for every
trajectory and spin manifold. `run.py` stores the dynamics data and diagnostic
plots under `FSSH2_/`.

`make_trajectory.py` reads `FSSH2_/mem_data.hdf`, converts coordinates from
Bohr to Angstrom, and writes:

- `trajectories/traj_0.xyz`
- `trajectories/traj_1.xyz`
- `trajectories/traj_2.xyz`

## Requirements

- A working Libra installation available through `PYTHONPATH`.
- OpenMolcas with its runtime environment configured.
- NumPy, HDF5/h5py, and Matplotlib.

Before running, update `model_params_molcas["exe"]` in `run.py` if the
OpenMolcas `pymolcas` executable is installed at a different location.

## Running

Run both commands from this directory and in this order:

```bash
python run.py
python make_trajectory.py
```

The second command requires the `FSSH2_/mem_data.hdf` file produced by the
first. If the molecule or number of atoms is changed in `run.py`, update the
`labels` list in `make_trajectory.py` accordingly.

The initial nuclear perturbations are sampled without a fixed random seed, so
separate runs will not produce identical trajectories.

## Generated data

The calculation creates `FSSH2_/`, `trajectories/`, an `al3/` scratch
directory, and `wd_molcas_*` OpenMolcas working directories. These files can be
large and are excluded from version control by the accompanying `.gitignore`.
