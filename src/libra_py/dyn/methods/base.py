"""Contract for coherent/ensemble method specializations of DynamicsEngine."""
from typing import Protocol


class DynamicsMethodImplementation(Protocol):
    """Numerical hooks; orchestration and lifecycle remain in dyn.engine.

    prepare rebuilds caches after external edits. advance uses the old time;
    finish uses the new time after the engine increments the ensemble clock.
    snapshots returns trajectory-id -> owned snapshot mapping. Implementations
    must not advance the clock, save files, or keep a second live state container.
    """
    def prepare(self, trajectories): ...
    def advance(self, trajectories, dt): ...
    def finish(self, trajectories, dt): ...
    def snapshots(self, trajectories, time): ...
