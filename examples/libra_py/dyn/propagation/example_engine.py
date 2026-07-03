"""
Educational example: DynamicsEngine and TSHEngine.

The TSH path currently runs TD-SE amplitudes and state-specific forces without
hopping decisions or momentum rescaling.

Run from the repository root:

    PYTHONPATH=src python examples/libra_py/dyn/propagation/example_engine.py
"""

import numpy as np

from libra_py.dyn.backends import backend
from libra_py.dyn.control_params import DynControlParams
from libra_py.dyn.core.storage import TensorStorage
from libra_py.dyn.core.trajectory import Trajectory
from libra_py.dyn.engine import DynamicsEngine, TSHEngine


def section(title):
    print(f"\n--- {title} ---")


def adiabatic_model(R, P, storage, traj):
    return {
        "H_adi": np.array([np.diag([0.0, 0.2])], dtype=complex),
        "NAC_adi": np.array([[[0.0, 0.03], [-0.03, 0.0]]], dtype=complex),
        "dH_adi": np.array([[[[0.1, 0.0], [0.0, 0.3]]]], dtype=complex),
        "DC1_adi": np.zeros((1, 1, 2, 2), dtype=complex),
    }


def make_system(active_state=0, amplitudes=None):
    storage = TensorStorage(
        backend=backend,
        ntraj=1,
        ndof=1,
        nstates=2,
        ntbf_initial=1,
    )
    storage.allocate_hamiltonian_derivatives()
    storage.iM[0, 0] = [1.0]
    storage.p[0, 0] = [1.0]
    storage.act_states[0, 0] = active_state
    storage.ampl_adi[0, 0] = amplitudes or [1.0 + 0.0j, 0.0 + 0.0j]

    traj = Trajectory(0)
    traj.tbf_ids = [0]
    return storage, traj


def main():
    section("Adiabatic dynamics")
    storage, traj = make_system(active_state=1)
    engine = DynamicsEngine(
        traj,
        storage,
        adiabatic_model,
        method="adiabatic",
        rep="adiabatic",
    )
    result = engine.step(0.1)
    print("time:", result.time)
    print("active state:", result.active_states)
    print("force:", result.forces)

    section("Ehrenfest dynamics")
    storage, traj = make_system(amplitudes=[0.6 + 0.0j, 0.8 + 0.0j])
    engine = DynamicsEngine(
        traj,
        storage,
        adiabatic_model,
        params=DynControlParams(force_method=2),
        method="ehrenfest",
        rep="adiabatic",
    )
    result = engine.step(0.1)
    print("force:", result.forces)
    print("amplitudes:", result.amplitudes)

    section("TSH preparation")
    storage, traj = make_system(active_state=0)
    engine = TSHEngine(traj, storage, adiabatic_model, rep="adiabatic")
    result = engine.step(0.1)
    print("active state:", result.active_states)
    print("amplitudes:", result.amplitudes)


if __name__ == "__main__":
    main()
