import numpy as np

from libra_py.dyn.core.storage import TensorStorage
from libra_py.dyn.core.tbf import TrajectoryBasisFunction
from libra_py.dyn.core.views import TensorView


def test_tbf_identity_lifecycle_and_metadata():
    storage = TensorStorage(
        backend=np,
        ntraj=2,
        ndof=2,
        nstates=2,
        ntbf_initial=1,
    )
    view = TensorView(storage, traj_id=1, tbf_id=0)
    tbf = TrajectoryBasisFunction(
        id=7,
        view=view,
        trajectory_id=1,
        metadata={"label": "initial"},
    )

    assert tbf.object_id == 7
    assert tbf.tbf_id == 0
    assert tbf.traj_id == 1
    assert tbf.is_alive_in_storage()

    tbf.view.q[:] = [0.0, 0.5]
    assert np.allclose(storage.q[1, 0], [0.0, 0.5])

    tbf.deactivate()
    assert not tbf.alive
    assert not bool(storage.alive[1, 0])

    tbf.activate()
    tbf.set_parent(3)
    tbf.set_spawn_time(12.5)

    summary = tbf.to_dict()
    assert summary["id"] == 7
    assert summary["tbf_id"] == 0
    assert summary["parent_id"] == 3
    assert summary["spawn_time"] == 12.5
    assert summary["metadata"] == {"label": "initial"}
