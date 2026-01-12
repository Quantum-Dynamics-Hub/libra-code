**High-level principles for module splitting**

Data structures ≠ algorithms

Topology ≠ numerics

Interfaces ≠ implementations

Stable core ≠ experimental physics

No circular imports


**Import direction rule (very important)**

Allowed direction:
core → propagation → spawning → experiments

Forbidden:
propagation → core
spawning → core

Core must never depend on physics.

