"""
Trajectory basis-function topology.

TrajectoryBasisFunction is intentionally light: it owns identity and lifecycle
metadata, while all numerical state lives in TensorStorage and is accessed
through TensorView.
"""

from .views import TensorView


class TrajectoryBasisFunction:
    """
    Trajectory Basis Function (TBF):
    a single nuclear wavepacket with electronic amplitudes.

    Parameters
    ----------
    id : int
        Object/registry id assigned by TrajectoryManager.
    view : TensorView
        Writable view into the shared TensorStorage arrays.
    trajectory_id : int
        Trajectory object id this TBF belongs to.
    parent_id : int, optional
        Object id of the parent TBF, if this TBF was spawned.
    spawn_time : float, optional
        Simulation time at which this TBF was spawned.
    """

    def __init__(
        self,
        id: int,
        view: TensorView,
        trajectory_id: int,
        parent_id=None,
        spawn_time=None,
        metadata=None,
    ):
        self.id = id
        self.view = view
        self.trajectory_id = trajectory_id

        # Bookkeeping
        self.alive = True
        self.parent_id = parent_id
        self.spawn_time = spawn_time
        self.metadata = dict(metadata or {})

    @property
    def object_id(self):
        """Object/registry id assigned by TrajectoryManager."""

        return self.id

    @property
    def tbf_id(self):
        """TensorStorage TBF slot index addressed by this object."""

        return self.view.tbf_id

    @property
    def storage_tbf_id(self):
        """Alias for the TensorStorage TBF slot index."""

        return self.tbf_id

    @property
    def traj_id(self):
        """TensorStorage trajectory index addressed by this object."""

        return self.view.traj_id

    @property
    def storage(self):
        """Backing TensorStorage object."""

        return self.view.storage

    def is_alive_in_storage(self):
        """Return True if this TBF slot is active in TensorStorage."""

        return bool(self.storage.alive[self.traj_id, self.tbf_id])

    def activate(self, update_storage=True):
        """
        Mark this TBF active.

        Parameters
        ----------
        update_storage : bool, optional
            If True, also updates TensorStorage.alive.
        """

        self.alive = True
        if update_storage:
            self.storage.alive[self.traj_id, self.tbf_id] = True

    def deactivate(self, update_storage=True):
        """
        Mark this TBF inactive.

        Parameters
        ----------
        update_storage : bool, optional
            If True, also updates TensorStorage.alive.
        """

        self.alive = False
        if update_storage:
            self.storage.alive[self.traj_id, self.tbf_id] = False

    def set_parent(self, parent_id):
        """Record the object id of the parent TBF."""

        self.parent_id = parent_id

    def set_spawn_time(self, time):
        """Record the simulation time at which this TBF was spawned."""

        self.spawn_time = time

    def to_dict(self):
        """Return a compact metadata summary."""

        return {
            "id": self.id,
            "trajectory_id": self.trajectory_id,
            "traj_id": self.traj_id,
            "tbf_id": self.tbf_id,
            "alive": self.alive,
            "alive_in_storage": self.is_alive_in_storage(),
            "parent_id": self.parent_id,
            "spawn_time": self.spawn_time,
            "metadata": dict(self.metadata),
        }

    def __repr__(self):
        return (
            "TrajectoryBasisFunction("
            f"id={self.id}, trajectory_id={self.trajectory_id}, "
            f"traj_id={self.traj_id}, tbf_id={self.tbf_id}, "
            f"alive={self.alive})"
        )
