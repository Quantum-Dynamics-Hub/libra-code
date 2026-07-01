"""
Trajectory/TBF registry.

TrajectoryManager owns the topology objects and the shared TensorStorage. It is
the place where object ids, storage slots, and lifecycle bookkeeping meet.
"""

from .storage import TensorStorage
from .tbf import TrajectoryBasisFunction
from .trajectory import Trajectory
from .views import TensorView


class TrajectoryManager:
    """
    Owns all TBFs, Trajectories, and shared tensor storage.
    """

    def __init__(self, storage: TensorStorage):
        self.storage = storage

        self.tbfs: dict[int, TrajectoryBasisFunction] = {}
        self.trajectories: dict[int, Trajectory] = {}

        self._next_tbf_id = 0
        self._next_traj_id = 0

    def create_trajectory(self, metadata=None) -> Trajectory:
        """
        Create and register a new trajectory.

        Parameters
        ----------
        metadata : dict, optional
            User metadata to attach to the Trajectory.
        """

        traj = Trajectory(self._next_traj_id)
        traj.metadata.update(metadata or {})
        self.trajectories[traj.id] = traj
        self._next_traj_id += 1
        return traj

    def get_trajectory(self, trajectory):
        """
        Return a Trajectory from either an object or an id.
        """

        if isinstance(trajectory, Trajectory):
            return trajectory
        return self.trajectories[trajectory]

    def get_tbf(self, tbf):
        """
        Return a TrajectoryBasisFunction from either an object or an object id.
        """

        if isinstance(tbf, TrajectoryBasisFunction):
            return tbf
        return self.tbfs[tbf]

    def create_tbf(
        self,
        trajectory: Trajectory,
        idx: int,
        parent_id=None,
        spawn_time=None,
        metadata=None,
        activate_storage=True,
    ) -> TrajectoryBasisFunction:
        """
        Create and register a TBF object for an existing storage slot.

        Parameters
        ----------
        trajectory : Trajectory or int
            Owning trajectory or trajectory id.
        idx : int
            TensorStorage TBF slot id.
        parent_id : int, optional
            Parent TBF object id.
        spawn_time : float, optional
            Simulation time at which this TBF was created.
        metadata : dict, optional
            User metadata to attach to the TBF object.
        activate_storage : bool, optional
            Whether to mark storage.alive[trajectory.id, idx] as True.
        """

        trajectory = self.get_trajectory(trajectory)

        if idx < 0 or idx >= self.storage.ntbf_capacity:
            raise IndexError(
                f"TBF storage slot {idx} is outside capacity "
                f"{self.storage.ntbf_capacity}"
            )

        view = TensorView(self.storage, trajectory.id, idx)
        tbf = TrajectoryBasisFunction(
            self._next_tbf_id,
            view,
            trajectory.id,
            parent_id=parent_id,
            spawn_time=spawn_time,
            metadata=metadata,
        )

        self.tbfs[tbf.id] = tbf
        trajectory.add_tbf(tbf)
        if activate_storage:
            tbf.activate(update_storage=True)

        self._next_tbf_id += 1
        return tbf

    def spawn_tbf(
        self,
        trajectory,
        parent=None,
        spawn_time=None,
        copy_from_parent=True,
        metadata=None,
    ) -> TrajectoryBasisFunction:
        """
        Allocate a new storage slot and register a TBF object for it.

        Parameters
        ----------
        trajectory : Trajectory or int
            Owning trajectory or trajectory id.
        parent : TrajectoryBasisFunction or int, optional
            Parent TBF object or object id.
        spawn_time : float, optional
            Simulation time at which this TBF is spawned.
        copy_from_parent : bool, optional
            If True and parent is supplied, copy core storage fields from the
            parent slot into the new slot.
        metadata : dict, optional
            User metadata for the spawned TBF.
        """

        trajectory = self.get_trajectory(trajectory)
        parent_tbf = self.get_tbf(parent) if parent is not None else None

        slot = self.storage.spawn(trajectory.id)
        tbf = self.create_tbf(
            trajectory,
            slot,
            parent_id=parent_tbf.id if parent_tbf is not None else None,
            spawn_time=spawn_time,
            metadata=metadata,
            activate_storage=False,
        )

        if parent_tbf is not None and copy_from_parent:
            self.copy_tbf_data(parent_tbf, tbf)

        return tbf

    def copy_tbf_data(self, source, target, fields=None):
        """
        Copy selected tensor slices from one TBF object to another.

        Parameters
        ----------
        source, target : TrajectoryBasisFunction or int
            Source and target TBF objects or object ids.
        fields : sequence of str, optional
            TensorStorage field names to copy. Defaults to core nuclear and
            electronic state fields.
        """

        source = self.get_tbf(source)
        target = self.get_tbf(target)
        fields = fields or (
            "q",
            "p",
            "iM",
            "f",
            "ampl_adi",
            "ampl_dia",
            "q_mm",
            "p_mm",
            "proj_adi",
            "dm_adi",
            "dm_dia",
            "act_states",
            "act_states_dia",
        )

        for name in fields:
            value = source.view.maybe(name)
            if value is not None:
                target.view.set(name, value)

    def deactivate_tbf(self, tbf, remove_from_trajectory=False):
        """
        Deactivate a TBF object and its storage slot.

        Parameters
        ----------
        tbf : TrajectoryBasisFunction or int
            TBF object or object id.
        remove_from_trajectory : bool, optional
            If True, also removes the storage slot from trajectory bookkeeping.
        """

        tbf = self.get_tbf(tbf)
        tbf.deactivate(update_storage=True)

        if remove_from_trajectory:
            traj = self.trajectories[tbf.trajectory_id]
            if traj.has_tbf(tbf.tbf_id):
                traj.remove_tbf(tbf.tbf_id)

        return tbf

    def activate_tbf(self, tbf):
        """
        Activate a TBF object and its storage slot.
        """

        tbf = self.get_tbf(tbf)
        tbf.activate(update_storage=True)
        traj = self.trajectories[tbf.trajectory_id]
        if not traj.has_tbf(tbf.tbf_id):
            traj.add_tbf(tbf)
        return tbf

    def remove_tbf(self, tbf):
        """
        Unregister a TBF object and deactivate its storage slot.
        """

        tbf = self.get_tbf(tbf)
        self.deactivate_tbf(tbf, remove_from_trajectory=True)
        return self.tbfs.pop(tbf.id)

    def active_tbf_ids(self, trajectory):
        """
        Return active storage TBF slots for a trajectory.
        """

        trajectory = self.get_trajectory(trajectory)
        return trajectory.active_tbf_ids(self.storage)

    def active_tbfs(self, trajectory):
        """
        Return active TBF objects for a trajectory.
        """

        trajectory = self.get_trajectory(trajectory)
        return [
            self.tbfs[trajectory.tbf_object_by_slot[slot]]
            for slot in trajectory.active_tbf_ids(self.storage)
        ]

    def deactivate_trajectory(self, trajectory, deactivate_tbfs=True):
        """
        Mark a trajectory inactive, optionally deactivating all its TBFs.
        """

        trajectory = self.get_trajectory(trajectory)
        trajectory.deactivate()

        if deactivate_tbfs:
            for object_id in list(trajectory.tbf_object_ids):
                self.deactivate_tbf(object_id)

        return trajectory

    def activate_trajectory(self, trajectory, activate_tbfs=False):
        """
        Mark a trajectory active, optionally reactivating all its TBFs.
        """

        trajectory = self.get_trajectory(trajectory)
        trajectory.activate()

        if activate_tbfs:
            for object_id in list(trajectory.tbf_object_ids):
                self.activate_tbf(object_id)

        return trajectory

    def compact_storage(self):
        """
        Compact TensorStorage and refresh TBF slot ids.

        This operation rebuilds trajectory/TBF slot bookkeeping after dropping
        globally inactive storage slots. TBF object ids are preserved.
        """

        old_slots = [
            slot
            for slot, keep in enumerate(self.storage.alive.any(axis=0))
            if bool(keep)
        ]
        slot_map = {old: new for new, old in enumerate(old_slots)}

        self.storage.compact()

        for traj in self.trajectories.values():
            new_tbf_ids = []
            new_object_by_slot = {}

            for old_slot in traj.tbf_ids:
                if old_slot not in slot_map:
                    continue

                new_slot = slot_map[old_slot]
                object_id = traj.tbf_object_by_slot[old_slot]
                tbf = self.tbfs[object_id]
                tbf.view = TensorView(self.storage, traj.id, new_slot)
                new_tbf_ids.append(new_slot)
                new_object_by_slot[new_slot] = object_id

            traj.tbf_ids = new_tbf_ids
            traj.tbf_object_by_slot = new_object_by_slot

        return slot_map

    def summary(self):
        """
        Return a compact manager summary.
        """

        return {
            "n_trajectories": len(self.trajectories),
            "n_tbfs": len(self.tbfs),
            "trajectories": {
                key: traj.to_dict()
                for key, traj in self.trajectories.items()
            },
            "tbfs": {
                key: tbf.to_dict()
                for key, tbf in self.tbfs.items()
            },
        }
