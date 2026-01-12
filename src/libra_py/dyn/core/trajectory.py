"""
Key insight:
Propagation, decoherence, pruning operate here, not on individual TBFs.

"""

from .tbf import TrajectoryBasisFunction

class Trajectory:
    """
    Electron–nuclear wavefunction represented as
    a coherent set of coupled TBFs.
    """

    def __init__(self, id: int):
        self.id = id
        self.tbf_ids: list[int] = []

        # Wavefunction-level properties
        self.norm = 1.0
        self.active = True

    def add_tbf(self, tbf: TrajectoryBasisFunction):
        self.tbf_ids.append(tbf.id)
        tbf.trajectory_id = self.id

    def remove_tbf(self, tbf_id: int):
        self.tbf_ids.remove(tbf_id)

