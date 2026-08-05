from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[4]
SRC = ROOT / "src"

if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from libra_py.dyn.backends import backend
from libra_py.dyn.control_params import DynControlParams
from libra_py.dyn.core.storage import TensorStorage
from libra_py.dyn.core.trajectory import Trajectory
from libra_py.dyn.engine import DynamicsEngine


def _storage():
    storage = TensorStorage(
        backend=backend,
        ntraj=1,
        ndof=1,
        nstates=2,
        ntbf_initial=1,
        ntbf_capacity=2,
    )
    storage.allocate_hamiltonian_derivatives()
    storage.iM[0, 0] = [1.0]
    storage.p[0, 0] = [1.0]
    storage.ampl_adi[0, 0] = [1.0 + 0.0j, 0.0 + 0.0j]
    storage.act_states[0, 0] = 1
    return storage


def _trajectory():
    traj = Trajectory(0)
    traj.tbf_ids = [0]
    return traj


def _adiabatic_model(R, P, storage, traj):
    return {
        "H_adi": np.array([np.diag([0.0, 0.2])], dtype=complex),
        "NAC_adi": np.zeros((1, 2, 2), dtype=complex),
        "dH_adi": np.array([[[[0.1, 0.0], [0.0, 0.3]]]], dtype=complex),
    }


def test_adiabatic_engine_uses_active_state_force_without_tdse():
    storage = _storage()
    traj = _trajectory()

    engine = DynamicsEngine(
        traj,
        storage,
        _adiabatic_model,
        method="adiabatic",
        rep="adiabatic",
    )
    result = engine.step(0.5)

    assert result.method == "adiabatic"
    assert result.timestep == 1
    np.testing.assert_allclose(result.forces, [[-0.3]])
    np.testing.assert_allclose(storage.p[0, 0], [0.85])
    np.testing.assert_allclose(storage.ampl_adi[0, 0], [1.0 + 0.0j, 0.0 + 0.0j])


def test_ehrenfest_engine_propagates_tdse_and_mean_field_force():
    storage = _storage()
    traj = _trajectory()
    storage.act_states[0, 0] = 0

    params = DynControlParams(force_method=2)
    engine = DynamicsEngine(
        traj,
        storage,
        _adiabatic_model,
        params=params,
        method="ehrenfest",
        rep="adiabatic",
    )
    result = engine.step(0.1)

    assert result.method == "ehrenfest"
    np.testing.assert_allclose(result.forces, [[-0.1]])
    np.testing.assert_allclose(abs(result.amplitudes[0, 0]), 1.0)
    np.testing.assert_allclose(storage.dm_adi[0, 0, 0, 0], 1.0 + 0.0j)


def test_dynamics_engine_tsh_runs_tdse_but_does_not_hop():
    storage = _storage()
    traj = _trajectory()
    storage.act_states[0, 0] = 1
    storage.ampl_adi[0, 0] = [0.0 + 0.0j, 1.0 + 0.0j]

    engine = DynamicsEngine(
        traj,
        storage,
        _adiabatic_model,
        method="tsh",
        rep="adiabatic",
    )
    result = engine.step(0.1)

    assert result.method == "tsh"
    np.testing.assert_array_equal(storage.act_states[0, 0:1], [1])
    np.testing.assert_allclose(abs(result.amplitudes[0, 1]), 1.0)
