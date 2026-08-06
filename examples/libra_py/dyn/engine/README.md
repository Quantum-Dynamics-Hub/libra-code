# DynamicsEngine examples

Each numbered example is self-contained and must be run from its own
directory. The scripts rely on the active Libra environment and do not modify
`sys.path`. Computed snapshots and plots remain nested inside the selected
example directory.

## 01_tully1_fssh_ehrenfest

This example compares twelve-trajectory FSSH and Ehrenfest propagation through
the Tully model 1 avoided crossing. FSSH uses derivative-coupling acceptance
option 20 and energy-conserving derivative-coupling rescaling option 200. Both
methods use the symmetric two-point vibronic-Hamiltonian integrator.

```bash
cd examples/libra_py/dyn/engine/01_tully1_fssh_ehrenfest
python compute.py
python plot.py
```

Snapshots are written to `output/tsh` and `output/ehrenfest`. The plotting
script creates energy, population, nuclear-motion, and trajectory-resolved
active-state PNG files under `plots/`.

## 02_tully1_tsh_decoherence

This example compares FSSH (`tsh_method=0`), GFSH (`tsh_method=1`), and FSSH2
(`tsh_method=7`). Each hopping method is run coherently, with instantaneous
ID-A, and with EDC-driven SDM, giving nine calculations initialized from the
same incoming ensemble.

```bash
cd examples/libra_py/dyn/engine/02_tully1_tsh_decoherence
python compute.py
python plot.py
```

For a quick smoke calculation:

```bash
python compute.py --ntraj 2 --nsteps 4 --save-stride 1
python plot.py
```

Snapshots are placed under `output/tsh_decoherence_comparison/`. The plotting
script creates:

- `plots/tsh_decoherence_comparison/surface_hopping_population.png`
- `plots/tsh_decoherence_comparison/electronic_population.png`
- `plots/tsh_decoherence_comparison/energy_drift.png`
- `plots/tsh_decoherence_comparison/active_states.png`

`energy_drift.png` plots `Etot_ave(t) - Etot_ave(0)`. The similarly named
saved quantity `dEtot_ave` is the ensemble standard deviation and is not an
energy-conservation diagnostic.
