import numpy as np

from libra_py.dyn.backends import backend
from libra_py.dyn.core.storage import TensorStorage
from libra_py.dyn.core.trajectory import Trajectory
from libra_py.dyn.hamiltonians import HamiltonianEngine


def assert_engine_builds_adiabatic(model, q_values):
    q_values = np.asarray(q_values, dtype=float)
    if q_values.ndim == 1:
        q_values = q_values.reshape(1, -1)
    ndof, npoints = q_values.shape
    nstates = model.evaluate(q_values)["H_dia"].shape[-1]
    storage = TensorStorage(
        backend=backend,
        ntraj=1,
        ndof=ndof,
        nstates=nstates,
        ntbf_initial=npoints,
        ntbf_capacity=npoints,
    )
    storage.allocate_hamiltonian_derivatives(der_lvl=1)
    storage.q[0, :, :] = q_values.T
    storage.p[0, :, :] = 1.0
    storage.iM[0, :, :] = 0.5
    traj = Trajectory(0)
    traj.tbf_ids = list(range(npoints))

    HamiltonianEngine(backend).evaluate(traj, storage, model, rep="adiabatic")

    expected_e = np.linalg.eigvalsh(storage.ham_dia[0, traj.tbf_ids].real)
    assert np.allclose(
        np.diagonal(storage.ham_adi[0, traj.tbf_ids], axis1=-2, axis2=-1).real,
        expected_e,
    )
    assert np.allclose(
        storage.hvib_adi[0, traj.tbf_ids],
        storage.ham_adi[0, traj.tbf_ids] - 1j * storage.nac_adi[0, traj.tbf_ids],
    )
