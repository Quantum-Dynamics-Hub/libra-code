# Python dynamics architecture

This package implements dynamics over shared trajectory-basis-function (TBF)
and tensor-storage infrastructure. Full simulations run through
`dyn.engine.DynamicsEngine`; methods specialize its numerical operations instead
of maintaining independent simulation engines. These principles apply to new
methods and to refactoring existing implementations.

## State and topology

- A TBF represents one nuclear wavepacket; its object stores identity, topology,
  metadata, and a view into shared storage.
- A `Trajectory` conceptually represents one electron–nuclear wavefunction;
  a batch represents independent wavefunctions. Preserve the existing
  independent-path initializer's `(storage, trajectory)` interface, which groups
  independent paths into slots for batched FSSH/Ehrenfest calculations. Do not
  interpret those slots as a coherent expansion.
- `TensorStorage` owns live numerical state. `TensorView` provides writable
  access; `TrajectoryManager` coordinates registration and topology changes.
  Do not keep a second live state container in a method or on TBF objects.
- Respect documented array axes: storage generally starts with
  `(ntraj, capacity, ...)`. Use actual active slots and electronic labels;
  do not assume contiguous slots, sorted labels, equal basis sizes, or fixed
  numbers of states/DOFs. Adapt kernel layouts explicitly at module boundaries.
- Optional arrays are allocated explicitly and must survive expansion, parent
  copying, compaction, and view remapping. Owned analysis snapshots and transient
  numerical workspaces are legitimate copies, distinct from live state.

See [core](core/README.md) for storage and lifecycle details.

## One engine, method-specific numerical operations

`DynamicsEngine` owns simulation orchestration: initialization, step/run
lifecycle, time and timestep updates, trajectory selection where supported,
history collection, saver cadence, and resource cleanup. Use its context manager
when attaching savers. Avoid a second method-specific run loop.

Method-specific physics belongs under `methods/<method>/`. QTAG implements the
[method protocol](methods/base.py) in `methods/qtag/method.py`:

- `prepare`: rebuild numerical caches for the selected trajectories.
- `advance`: perform the part of the step evaluated at the old time.
- `finish`: complete the step after the engine advances the ensemble clock.
- `snapshots`: provide owned observable records keyed by trajectory ID.

These hooks must not increment the clock, save files, or own duplicate live
state. Reusable numerical kernels belong in shared domain modules, even when
first introduced for one method. Method modules compose those kernels and retain
only method-specific choices and splitting rules.

Current selectors include `qtag`, `ehrenfest`, `tsh`, and `adiabatic`.
FSSH is a sub-method of `tsh`, selected by `tsh_method=0`, not a separate
engine selector. Other hopping variants use their corresponding algorithm numbers.
QTAG accepts a manager and returns `DynamicsRunResult` with per-trajectory
histories. Existing independent-path methods accept `(trajectory, storage,
model_fn)` and return lists of `StepResult`. Preserve these supported interfaces;
engine integration does not imply identical initialization or result types.

## Module responsibilities and imports

| Location | Responsibility |
| --- | --- |
| `core` | Storage, views, identity, topology, and lifecycle; no physics |
| `backends` | Array creation, dtype/device handling, backend operations |
| `initialization/ensembles.py` | Shared sampling and independent/coherent ensemble construction |
| `gwp` | Gaussian integrals, packet evaluation, potential moments, storage adapters |
| `models` | Model potentials and analytic derivatives/Hessians |
| `hamiltonians` | Model evaluation, representation handling, Hamiltonian assembly |
| `algorithms` | Reusable numerical algorithms, including NPI and quantum-momentum/gradient estimation in `quantum_potential.py` |
| `propagation` | Reusable nuclear/electronic and nonorthogonal propagation kernels |
| `hopping`, `decoherence`, `spawning` | Domain-specific algorithms and criteria |
| `observables` | Physical estimators, reconstruction, and reusable plotting |
| `savers` | Serialization, disk formats, buffering, and loading |
| `utils` | Small shared utilities such as metric linear algebra |
| `methods/<method>` | Method-specific numerical composition |
| `engine.py` | Shared simulation orchestration |

Here **A → B means A imports/depends on B**. Higher-level orchestration and
physics may import lower-level interfaces, for example
`engine → methods → propagation/gwp/core`. `core` may use backend support but
must not import the engine, methods, models, or physics algorithms. Shared kernels
must not depend on a particular method or its run loop. Avoid circular imports;
there is no requirement to force every domain module into one linear chain.

Keep implementations in named modules, not `__init__.py`; package initializers
provide documentation and intentional public exports. Put reusable models in
`models`, not example-local copies. Reuse existing helpers before introducing a
parallel implementation. In particular, extend `initialization/ensembles.py`
instead of restoring the removed `initialization/gaussian.py` module.

## Numerical semantics, backends, and performance

Coherent Gaussian coefficients and independent-path electronic amplitudes have
different meanings. QTAG uses overlap-aware norm `C†SC`, energy `C†HC`, and
state-restricted populations; its saved energy is unnormalized. Never normalize
coefficients independently per Gaussian or substitute independent-path
estimators. Label representation and normalization in outputs: QTAG examples
report coherent diabatic populations, while Tully examples report adiabatic
ensemble populations and, for hopping, separate active-state fractions.

Document units, shapes, and width conventions. QTAG storage uses
`exp[-alpha*(r-q)^2/2]`, while the integral/initial-wavepacket convention uses
`exp[-a*(r-q)^2]`, so `alpha = 2*a` for the same packet. See the
[QTAG module map and conventions](methods/qtag/README.md).

Use shared backend helpers to preserve NumPy/PyTorch dtype, device, and autograd
where supported. Do not silently detach tensors or round-trip through NumPy in
numerical kernels. Conversion to CPU/NumPy belongs at explicit plotting or disk
I/O boundaries. Document backend limitations rather than claiming support that
has not been implemented. Pure-Python dynamics must not require C++ bindings or
SymPy at runtime.

Preserve vectorized operations and reuse valid overlap/Hamiltonian evaluations
between uninterrupted steps. Rebuild caches when external state or topology
changes; do not trade correctness for stale cached data. Avoid redundant history
copies and model evaluations. Use existing buffered savers when appropriate;
keep I/O cost separate from propagation cost when measuring performance.

## Examples and verification

Runnable method examples live under
[`examples/libra_py/dyn/methods`](../../../examples/libra_py/dyn/methods/README.md):
`1_qtag`, `2_tsh`, and `3_ehrenfest`. Prefer a self-contained `run.py` with explicit,
commented parameter construction, shared initialization, a `DynamicsEngine` call,
shared savers, and headless individual PNG outputs. Keep script structure, CLI
options, summary printouts, and comparable properties consistent across methods,
while explaining physical differences. The decoherence example combines
computation and plotting in one script with optional stage selectors.

Tests mirror source locations under `unittests/test_libra_py/test_dyn` and should
run from the repository root and their own directories. Cover numerical formulas,
backend agreement, and affected lifecycle/API contracts. Use independent
references, finite differences, or conservation checks as appropriate. The former
standalone QTAG implementation lives only in `test_dyn/qtag_reference` as a
regression oracle; production code and examples must not import it. Keep that
oracle independent rather than changing it to reproduce a new implementation.
Optional C++ comparisons and long tutorial regressions should remain explicit.

When moving examples or APIs, update documentation links and regression-test paths
in the same change. Run checks proportionate to the change; short example runs
should verify saved records and plots, not just successful imports.

## Guidance for future agents

Read this document before changing the dynamics architecture. Scoped `AGENTS.md`
files in this source package and the corresponding example/test trees point to
these principles. Preserve the common engine and shared numerical modules when
adding methods; do not recreate the standalone QTAG architecture or duplicate
initializers. Update this document when an intentional architectural change
supersedes a principle, and explain the compatibility and numerical consequences.
