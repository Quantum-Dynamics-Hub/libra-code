from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[4]
SRC = ROOT / "src"

if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from libra_py.dyn.backends import backend
from libra_py.dyn.core.storage import TensorStorage
from libra_py.dyn.core.trajectory import Trajectory
from libra_py.dyn.propagation import (
    drift,
    kick,
    state_specific_forces,
    tdse_step,
    update_density,
    velocity,
)


def _storage(ndof=2):
    storage = TensorStorage(
        backend=backend,
        ntraj=1,
        ndof=ndof,
        nstates=2,
        ntbf_initial=1,
        ntbf_capacity=2,
    )
    storage.allocate_hamiltonian_derivatives()
    return storage


def _trajectory():
    traj = Trajectory(0)
    traj.tbf_ids = [0]
    return traj


def test_nuclear_drift_kick_and_velocity():
    storage = _storage()
    traj = _trajectory()
    storage.q[0, 0] = [0.0, 1.0]
    storage.p[0, 0] = [2.0, -1.0]
    storage.iM[0, 0] = [0.5, 0.25]
    storage.f[0, 0] = [-0.2, 0.4]

    np.testing.assert_allclose(velocity(storage, traj), [[1.0, -0.25]])
    drift(storage, traj, 0.1)
    kick(storage, traj, 0.1)

    np.testing.assert_allclose(storage.q[0, 0], [0.1, 0.975])
    np.testing.assert_allclose(storage.p[0, 0], [1.98, -0.96])


def test_tdse_step_updates_amplitudes_and_density():
    storage = _storage(ndof=1)
    traj = _trajectory()
    storage.ampl_adi[0, 0] = [1.0 + 0.0j, 0.0 + 0.0j]
    storage.hvib_adi[0, 0] = np.diag([0.0, 0.2])

    tdse_step(traj, storage, dt=0.5, rep="adiabatic", hamiltonian_type="vibronic")
    update_density(storage, traj, rep="adiabatic")

    np.testing.assert_allclose(storage.ampl_adi[0, 0], [1.0 + 0.0j, 0.0 + 0.0j])
    np.testing.assert_allclose(storage.dm_adi[0, 0, 0, 0], 1.0 + 0.0j)


def test_state_specific_forces_preserve_multi_dof_shape():
    storage = _storage(ndof=2)
    traj = _trajectory()
    storage.act_states[0, 0] = 1
    storage.d1ham_adi[0, 0] = np.array(
        [
            [[0.1, 0.0], [0.0, 0.3]],
            [[-0.2, 0.0], [0.0, 0.4]],
        ],
        dtype=complex,
    )

    forces = state_specific_forces(storage, traj)

    assert forces.shape == (1, 2)
    np.testing.assert_allclose(forces, [[-0.3, -0.4]])
