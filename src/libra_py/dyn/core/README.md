# dyn.core

The `core` package is the stable foundation of the Python dynamics prototype.
It defines data layout and topology only. Physics algorithms, propagation,
force evaluation, decoherence models, and spawning criteria should live in
higher-level modules.


## Design Goals

- Keep heavy numerical data in shared tensor buffers.
- Keep Python objects lightweight: identity, topology, metadata, and views.
- Use names that match the existing C++ headers where practical, especially
  `src/dyn/dyn_variables.h` and `src/nhamiltonian/nHamiltonian.h`.
- Make tensor indexing explicit and batch-friendly.
- Avoid data copies unless a manager operation intentionally copies one TBF
  slot into another.
- Keep method-specific memory allocation explicit, because some arrays are very
  large.


## Core Concepts

### TensorStorage

`TensorStorage` owns the numerical arrays. Its leading dimensions are:

```text
(ntraj, ntbf_capacity, ...)
```

where `ntraj` indexes independent trajectories and `ntbf_capacity` indexes
trajectory basis-function slots. The logical number of allocated TBF slots is
tracked by `ntbf`; the boolean mask `alive[ntraj, ntbf_capacity]` records which
slots are active for each trajectory.

The storage class follows C++ naming:

```text
q, p, iM, f
ampl_adi, ampl_dia
dm_adi, dm_dia
act_states, act_states_dia
ham_adi, ham_dia
hvib_adi, hvib_dia
nac_adi, nac_dia
basis_transform
time_overlap_adi, time_overlap_dia
```

For now, `ndia == nadi == nstates`; separate adiabatic and diabatic fields are
still kept so both representations can coexist.


### Allocation Policy

Constructor allocation is intentionally limited to common/core data:

- `allocate_nuclear_vars()`
- `allocate_electronic_vars()`
- `allocate_hamiltonian_vars()`

Method-specific arrays are allocated only by explicit calls:

```python
storage.allocate_afssh()
storage.allocate_fssh2()
storage.allocate_hamiltonian_derivatives(der_lvl=1)
storage.allocate_simple_decoherence()
```

This mirrors the C++ `allocate_*` style and avoids accidental allocation of
large tensors such as `d2ham_*`, which scale as `ndof * ndof`.

Status flags such as `electronic_vars_status`, `afssh_vars_status`, and
`dc1_adi_mem_status` are kept to make allocation state explicit.


### TensorView

`TensorView` is a writable slice into `TensorStorage` for one
`(trajectory, TBF slot)` pair. It does not own data.

Example:

```python
view.q[:] = [0.0, 0.1, 0.2]
view.ampl_adi[0] = 1.0
view.ham_adi[:] = H
```

These writes directly update the shared storage arrays.

Generic access is also available:

```python
view.get("q")
view.set("act_states", 0)
view.maybe("dR")
```

Optional fields raise a clear `AttributeError` if accessed before their
allocation function is called. `maybe()` returns `None` instead.


### TrajectoryBasisFunction

`TrajectoryBasisFunction` is topology and metadata for one nuclear wavepacket.
It stores:

- object/registry id
- owning trajectory id
- `TensorView`
- parent/spawn metadata
- lifecycle flag

Important distinction:

```text
tbf.id     -> object id used by TrajectoryManager
tbf.tbf_id -> TensorStorage slot index
```

This distinction matters because the storage slot is used for tensor indexing,
while the object id is used for registries and metadata.


### Trajectory

`Trajectory` represents one electron-nuclear wavefunction. It keeps a coherent
set of TBF storage slots:

```python
trajectory.tbf_ids
```

These are storage slot ids, so they can be used directly:

```python
storage.q[trajectory.id, trajectory.tbf_ids]
```

It also tracks TBF object ids separately:

```python
trajectory.tbf_object_ids
trajectory.tbf_object_by_slot
```

Wavefunction-level metadata such as `norm` and `active` live on the trajectory.


### TrajectoryManager

`TrajectoryManager` owns the object registries and connects topology to
storage. It handles:

- trajectory creation
- TBF object creation
- spawning a new storage slot
- copying data from parent TBFs
- activating/deactivating TBFs and trajectories
- removing TBFs
- compacting storage and remapping views
- summaries of topology state

The manager is the correct place to coordinate storage-slot ids and object ids.
Direct physics logic should not be added here.


## Lifecycle

A typical workflow is:

```python
storage = TensorStorage(backend=np, ntraj=10, ndof=6, nstates=2)
manager = TrajectoryManager(storage)

traj = manager.create_trajectory()
tbf0 = manager.create_tbf(traj, idx=0)

tbf0.view.q[:] = q0
tbf0.view.p[:] = p0
tbf0.view.ampl_adi[0] = 1.0

tbf1 = manager.spawn_tbf(traj, parent=tbf0, spawn_time=time)
manager.deactivate_tbf(tbf1)
manager.compact_storage()
```


## Module Boundaries

Allowed in `core`:

- tensor allocation
- tensor slicing
- object identity
- topology bookkeeping
- lifecycle flags
- registry operations

Not allowed in `core`:

- force evaluation
- electronic propagation
- nuclear propagation
- hopping decisions
- decoherence formulas
- spawning criteria
- model Hamiltonian construction

Those belong in `propagation`, `hamiltonians`, `spawning`, `decoherence`, or
other higher-level modules.


## Files

`storage.py`
: Tensor allocation and allocation-status flags.

`views.py`
: Writable per-TBF views into storage.

`tbf.py`
: TBF identity, metadata, lifecycle, and view ownership.

`trajectory.py`
: Wavefunction-level topology and TBF-slot bookkeeping.

`manager.py`
: Registry and lifecycle coordination across storage, trajectories, and TBFs.


## Examples

The files below are runnable demonstrations:

```bash
python test2.py  # TensorStorage
python test3.py  # TensorView
python test4.py  # TrajectoryBasisFunction
python test5.py  # Trajectory
python test6.py  # TrajectoryManager
```

Run them from `src/libra_py/dyn/core`.
