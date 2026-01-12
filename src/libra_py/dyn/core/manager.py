"""
 central authority
 easy spawning and deletion
 compatible with batching

"""

from .storage import TensorStorage
from .tbf import TrajectoryBasisFunction
from .trajectory import Trajectory

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

    def create_trajectory(self) -> Trajectory:
        traj = Trajectory(self._next_traj_id)
        self.trajectories[traj.id] = traj
        self._next_traj_id += 1
        return traj

    def create_tbf(self, trajectory: Trajectory, idx: int) -> TrajectoryBasisFunction:
        view = TensorView(self.storage, idx)
        tbf = TrajectoryBasisFunction(self._next_tbf_id, view, trajectory.id)

        self.tbfs[tbf.id] = tbf
        trajectory.add_tbf(tbf)

        self._next_tbf_id += 1
        return tbf

