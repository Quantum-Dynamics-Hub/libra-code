"""
Trajectory topology.

Trajectory represents one electron-nuclear wavefunction and keeps track of the
TensorStorage TBF slots that currently belong to it.
"""

from .tbf import TrajectoryBasisFunction


class Trajectory:
    """
    Electron–nuclear wavefunction represented as
    a coherent set of coupled TBFs.
    """

    def __init__(self, id: int):
        self.id = id

        # TensorStorage TBF slot ids. These are the indices used to slice
        # storage arrays, e.g. storage.q[trajectory.id, trajectory.tbf_ids].
        self.tbf_ids: list[int] = []

        # TBF object/registry ids, useful for manager bookkeeping.
        self.tbf_object_ids: list[int] = []

        # Map TensorStorage TBF slot id -> TBF object id.
        self.tbf_object_by_slot: dict[int, int] = {}

        # Wavefunction-level properties.
        self.norm = 1.0
        self.active = True
        self.metadata = {}

    def add_tbf(self, tbf: TrajectoryBasisFunction):
        """
        Attach a TBF object to this trajectory.

        The trajectory stores the TBF storage slot in tbf_ids because the
        propagation and Hamiltonian layers use these values for tensor
        indexing.
        """

        if tbf.traj_id != self.id:
            raise ValueError(
                f"TBF view points to trajectory {tbf.traj_id}, "
                f"but this trajectory id is {self.id}"
            )

        slot = tbf.tbf_id
        if slot not in self.tbf_ids:
            self.tbf_ids.append(slot)

        if tbf.id not in self.tbf_object_ids:
            self.tbf_object_ids.append(tbf.id)

        self.tbf_object_by_slot[slot] = tbf.id
        tbf.trajectory_id = self.id

    def remove_tbf(self, tbf_id: int):
        """
        Remove a TBF storage slot from this trajectory.

        Parameters
        ----------
        tbf_id : int
            TensorStorage TBF slot id, not the TBF object id.
        """

        self.tbf_ids.remove(tbf_id)
        object_id = self.tbf_object_by_slot.pop(tbf_id, None)
        if object_id in self.tbf_object_ids:
            self.tbf_object_ids.remove(object_id)

    def has_tbf(self, tbf_id: int):
        """Return True if this trajectory contains a TBF storage slot."""

        return tbf_id in self.tbf_ids

    def active_tbf_ids(self, storage=None):
        """
        Return active TBF storage slots.

        If storage is supplied, TensorStorage.alive is used as the source of
        truth; otherwise the current bookkeeping list is returned.
        """

        if storage is None:
            return list(self.tbf_ids)

        return [
            tbf_id
            for tbf_id in self.tbf_ids
            if bool(storage.alive[self.id, tbf_id])
        ]

    def deactivate(self):
        """Mark the trajectory inactive."""

        self.active = False

    def activate(self):
        """Mark the trajectory active."""

        self.active = True

    def set_norm(self, value):
        """Set the trajectory wavefunction norm."""

        self.norm = float(value)

    def __len__(self):
        return len(self.tbf_ids)

    def __iter__(self):
        return iter(self.tbf_ids)

    def to_dict(self):
        """Return a compact metadata summary."""

        return {
            "id": self.id,
            "tbf_ids": list(self.tbf_ids),
            "tbf_object_ids": list(self.tbf_object_ids),
            "norm": self.norm,
            "active": self.active,
            "metadata": dict(self.metadata),
        }

    def __repr__(self):
        return (
            "Trajectory("
            f"id={self.id}, tbf_ids={self.tbf_ids}, "
            f"norm={self.norm}, active={self.active})"
        )
