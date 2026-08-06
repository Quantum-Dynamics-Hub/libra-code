# Unified dynamics engine example

Run twelve trajectories through the Tully model 1 avoided crossing using both
FSSH and Ehrenfest dynamics, and save step-wise NPZ snapshots:

```bash
cd examples/libra_py/dyn/engine
python example_dynamics_engine.py
```

Results are written under `output/tsh` and `output/ehrenfest` using
`FaultTolerantSaver`. Each step is a self-contained NPZ file and each output
directory has a JSON-lines manifest. This reproducible example opens the saver
with `mode="w"`, replacing snapshots from its previous run; use the default
`mode="a"` for restart/append workflows.

The FSSH calculation pairs derivative-coupling acceptance (option 20) with
derivative-coupling momentum rescaling (option 200), so accepted hops conserve
the trajectory's kinetic plus active-surface energy. Both calculations use the
symmetric two-point vibronic-Hamiltonian electronic integrator (option 4).

Create plots from those saved results:

```bash
python plot_dynamics.py
```

The plotting script writes:

- `plots/energies.png`
- `plots/populations.png`
- `plots/nuclear_motion.png`
- `plots/active_state.png` (a trajectory-by-trajectory active-state map)

The active Libra environment must already provide `libra_py`; neither script
modifies `sys.path`.
