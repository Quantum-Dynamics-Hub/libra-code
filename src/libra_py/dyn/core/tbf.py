"""
TBF responsibilities:
- holds physical state via TensorView
- knows which Trajectory it belongs to
- can spawn / die

"""

from .views import TensorView

class TrajectoryBasisFunction:
    """
    Trajectory Basis Function (TBF):
    a single nuclear wavepacket with electronic amplitudes.
    """

    def __init__(self, id: int, view: TensorView, trajectory_id: int):
        self.id = id
        self.view = view
        self.trajectory_id = trajectory_id

        # Bookkeeping
        self.alive = True
        self.parent_id = None
        self.spawn_time = None

