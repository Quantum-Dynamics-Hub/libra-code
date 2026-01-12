# 1. Design principles (explicit)

- TrajectoryBasisFunction (TBF) = nuclear wavepacket

- Trajectory = one electron–nuclear wavefunction

- Batch = ensemble of independent wavefunctions

- Heavy numerics live in shared tensor buffers

- Objects store topology and views, not data copies

- Tensor ranks are not hard-coded



# 2. High-level principles for module splitting

Data structures ≠ algorithms

Topology ≠ numerics

Interfaces ≠ implementations

Stable core ≠ experimental physics

No circular imports


# 3. Import direction rule (very important)

Allowed direction:
core → propagation → spawning → experiments

Forbidden:
propagation → core
spawning → core

Core must never depend on physics.

