## Integrated QTAG

QTAG is a specialization of `dyn.engine.DynamicsEngine`. The previous standalone
implementation has moved to `unittests/test_libra_py/test_dyn/qtag_reference` as
an isolated test oracle. No production code imports it.

| Responsibility | Shared location |
| --- | --- |
| Tensor allocation, views, lifecycle | `core/storage.py`, `views.py`, `manager.py` |
| Array placement and dtype/device helpers | `backends/arrays.py`, `torch_backend.py` |
| Gaussian evaluation, log-gradients, overlap/kinetic | `gwp/packets.py`, existing `gwp/integrals.py` |
| BAT/LHA moments, exact Gaussian coupling | `gwp/potentials.py` |
| Storage slices and coherent coefficient writes | `gwp/storage.py` |
| Local model evaluation and pair-matrix assembly | `hamiltonians/engine.py`, `hamiltonians/gaussian.py` |
| Analytic potentials and Hessians | `models/holstein.py`, `models/ferretti.py` |
| Sampling and projection into a Gaussian basis | `initialization/ensembles.py` |
| Positive-subspace metric solve | `utils/linalg.py` |
| Generalized electronic evolution | `propagation/nonorthogonal.py`, shared `electronic.hermitian_propagator` |
| Quantum momentum and logarithmic-gradient fitting | `algorithms/quantum_potential.py` |
| Basis updates and coefficient projection | `propagation/gaussian.py`, shared `nuclear.drift_coordinates` |
| Coherent populations, energy and density reconstruction | `observables/gaussian.py` |
| Individual PNG plots | `observables/plotting.py` |
| Disk snapshots | Existing `savers.FaultTolerantSaver` and HDF5 savers |
| Symmetric splitting and ensemble coordination | `engine.py`, `methods/qtag/method.py` |

### State ownership and conventions

A `Trajectory` represents one coherent wavefunction. Each of its registered
TBFs is one Gaussian on one diabatic electronic state. Independent wavefunctions
occupy different leading trajectory rows in `TensorStorage`.

- `q`, `p`, `iM`, `gaussian_alpha`, `gaussian_phase`: `(ntraj,capacity,ndof)`.
- `act_states_dia`: `(ntraj,capacity)`, one electronic label per packet.
- `ampl_dia`: `(ntraj,capacity,nstates)`; only the labeled entry in each row is
  populated, holding the **coherent expansion coefficient**. It is not normalized
  separately per packet. `dm_dia` holds the corresponding local outer product.
- Local model fields `ham_dia`, `d1ham_dia`, `d2ham_dia` use the existing storage
  definitions and are written by `HamiltonianEngine`.
- Pair overlap and Hamiltonian matrices are transient, sized to the active
  packets. They are not confused with the per-TBF electronic matrices.

Gaussian fields are allocated explicitly with `allocate_gaussian_vars()`.
They follow storage expansion, parent copying, compaction and `TensorView`
access. The initializer handles allocation and registration automatically.
No numerical state is stored on TBF topology objects or in a second live array
container. History snapshots are copies for analysis.

For each DOF, the normalized packet is
`(alpha/pi)**0.25 * exp[-alpha*(r-q)**2/2 + i*p*(r-q) + i*phase]`.
Thus `gaussian_alpha` is twice the width passed to `gwp.integrals`.
Initial `wfc_a0` is a GWP width; the projected initial wavefunction uses
`exp[-wfc_a0*(r-q0)**2]`. Coordinates, momenta, masses and times use atomic units.
The array kernels use `(ndof,npackets)`, matching the integral API; transposition
occurs at the storage boundary.

Coherent observables include inter-packet overlaps: norm is `C†SC`, energy is
`C†HC`, and populations are state-restricted quadratic forms. Do not substitute
the independent-path population/energy estimators or normalize each TBF's
coefficient. The saved energy is unnormalized, matching the reference method.

### Usage

```python
from libra_py.dyn.initialization import make_gaussian_ensemble
from libra_py.dyn.engine import DynamicsEngine
from libra_py.dyn.models import FerrettiModel
from libra_py.dyn.savers import FaultTolerantSaver

initial = dict(
    backend="numpy",  # "torch" selects tensor storage; device/dtype are optional
    states=[0, 1], initial_state=0,
    wfc_q0=[0., 2.], wfc_p0=[0., 0.],
    wfc_a0=[1., 1.], wfc_s0=[0., 0.],
    grid_dims=[5, 5], alp_scl=[6., 6.], rho_cut=1e-12,
)
manager = make_gaussian_ensemble(initial, [1/20000., 1/6667.])
model = FerrettiModel(params={"Ky": .01}, backend="numpy")
engine = DynamicsEngine(manager=manager, model_fn=model, method="qtag", params=dict(
    dt=1., qtag_pot_approx_method=0, a_update_method=1,
    q_sync_method=1, p_sync_method=0, a_sync_method=0,
    linfit_beta=1e-3, decpl_den=.01,
))
result = engine.run(100, savers=FaultTolerantSaver("output/snapshots"))
history = result.histories[0]
# Actual current state remains accessible through manager.storage / TBF views.
```

For multiple wavefunctions, supply a list of initial-condition dictionaries.
They may have different packet counts but must share dimension, backend, dtype,
device and electronic storage dimension. `savers` then maps trajectory ids to
separate saver instances. The global timestep advances once per ensemble step;
models see the old timestep for the first half-step and the new timestep for
the final half-step. The model callback follows the standard
`model(R, P, storage, traj)` contract, with `R.shape=(nactive,ndof)`.

`methods.qtag.run_qtag(manager, model, params)` is a convenience wrapper using
`params['nsteps']`. `engine.step()` allows external lifecycle changes between
steps. Inactive trajectories and pruned slots are skipped. Automatic spawning
criteria and physical redistribution/projection of coefficients at a manual
spawn are not supplied. In particular, copying a parent's coefficient to a
new packet alone does not preserve the represented wavefunction.

### Numerical controls and limitations

Potential methods are 0 BAT, 1 LHA, 2 LHA with exact Gaussian coupling, and
3 BAT with exact coupling. LHA requires `d2H_dia`; no finite-difference Hessian
is silently inserted. Methods 2/3 require `exact_coupling=(A,B,C)` in dynamics
parameters or `ex_cpl_A/B/C` in the model parameters. They represent the same
`sum_d A exp[-B*(r_d-C)**2]` coupling for every off-diagonal pair; do not select
this option for a model with a different coupling, such as Ferretti.

The momentum/basis options follow the preserved QTAG implementation:
`mom_calc_type`, `d_weight`, `momentum_fit`, `linfit_beta`, `fit_rcond`,
`node_tol`, `q/p/a/s_update_method`, `q/p/a/s_sync_method`, and `decpl_den`.
`momentum_fit='direct'` is the default unless `linfit_beta` is supplied;
`momentum_fit='legacy', linfit_beta=1e-3` reproduces the tutorials. The finite
iterative stopping tolerance affects basis motion and must not be silently
replaced with an exact solve. Synchronization requires equal packet counts on
the participating surfaces; disable synchronization for unequal counts.
Phase evolution is unsupported (use `s_update_method=0`); stored initial phases
are preserved. Nonpositive widths raise instead of being clipped.

`overlap_rcond=1e-12` discards small metric eigenvalues. Projection can lose
norm; it is never hidden by renormalizing during propagation. A coherent
trajectory requires common positive inverse masses for all its packets.

Torch arrays retain device, precision and autograd graphs. Gradients are valid
only at fixed spectral rank and fixed discrete synchronization/fit choices;
eigenvector derivatives at repeated eigenvalues can be undefined. Use
`torch.no_grad()` for full trajectories when derivatives are not needed.
The new backend facade supports this path; it does not convert unrelated
NumPy-only dynamics methods into Torch implementations.

Histories include the initial state and have `(time,npackets,ndof)` nuclear
arrays. `states`, `slot_ids` and `target_states` are saved per record. Dense
history stacking assumes fixed topology during `run`; `step` plus individual
saved records supports external topology changes. Set `keep_history=False`
for disk-only recording. Savers detach tensors only at the disk boundary.

### Validation

Tests are under matching `unittests/test_libra_py/test_dyn` subdirectories.
Short tests compare every stored step and reconstructed density against the
test-only `qtag_reference`, cover storage lifecycle and independent trajectories,
finite-difference Hessians/gradients, and generic saving/plotting. A subprocess
blocks reference, C++ and SymPy imports to verify runtime independence.

Full tutorial parity is opt-in:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 QTAG_FULL_TESTS=1 PYTHONPATH=src \
  python -m pytest -q unittests/test_libra_py/test_dyn/test_methods/test_qtag/test_engine_qtag.py
```

Both the 800-step Holstein and 100-step Ferretti runs have been compared with
the preserved implementation on NumPy and CPU Torch, including all five density
snapshots. The runnable scripts live in `examples/libra_py/dyn/methods/1_qtag/1_basics`.

### Performance

An uninterrupted `run()` reuses final S/H matrices for the next first half-step,
as the reference method does. Model callbacks are assumed to be pure functions
of the current physical state/time during that run; saver callbacks must not
mutate live dynamics state. Standalone `step()` always rebuilds matrices, so
external edits between calls remain supported. History recording uses the
already-owned snapshot arrays instead of copying them a second time.

### Common engine interface

`DynamicsEngine(manager, model, method="qtag", params=...)` and the explicit
keyword form above are equivalent. `methods.qtag.QTAGEngine` is only a convenience
constructor returning the common engine. `algorithms.qtag` has been removed.

The engine owns initialization, active-trajectory selection, time/timestep,
recording, `nprint` saver cadence, saver mappings, and context-manager cleanup.
`initialize()` prepares and saves the initial record. `step(dt=None)` permits
external storage/topology edits and rebuilds caches; `run(nsteps=None, dt=None,
savers=None, keep_history=True)` reuses pair matrices internally. The timestep
can change between calls without resetting physical time. Energy changes refer
to the first recorded state of that engine. A new run includes its starting
snapshot, even when continuing. Dense histories assume fixed topology per run.

The protocol in `methods/base.py` exposes prepare, advance, finish and snapshots.
The engine increments its clock between advance and finish so every active
trajectory sees the same old/new model timestep. QTAG owns the splitting order,
metric cutoffs, temporary pair matrices and coherent estimators. Existing
adiabatic/Ehrenfest/TSH constructors and return values are unchanged.
