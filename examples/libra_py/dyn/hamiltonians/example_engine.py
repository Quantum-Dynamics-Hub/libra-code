"""
Educational example: HamiltonianEngine.

HamiltonianEngine writes model results directly into TensorStorage and can
return active Hamiltonian slices for propagation setup.

Run from the repository root:

    PYTHONPATH=src python examples/libra_py/dyn/hamiltonians/example_engine.py
"""

import numpy as np

from libra_py.dyn.backends import backend
from libra_py.dyn.core.storage import TensorStorage
from libra_py.dyn.core.trajectory import Trajectory
from libra_py.dyn.hamiltonians import HamiltonianEngine


def section(title):
    print(f"\n--- {title} ---")


def adiabatic_model(R, P, storage, traj):
    return {
        "H_adi": np.array([np.diag([0.0, 0.2])], dtype=complex),
        "NAC_adi": np.array([[[0.0, 0.03], [-0.03, 0.0]]], dtype=complex),
    }


def main():
    storage = TensorStorage(
        backend=backend,
        ntraj=1,
        ndof=1,
        nstates=2,
        ntbf_initial=1,
    )
    traj = Trajectory(0)
    traj.tbf_ids = [0]

    engine = HamiltonianEngine(backend)

    section("Evaluate model")
    engine.evaluate(traj, storage, adiabatic_model, rep="adiabatic")
    print("ham_adi:\n", storage.ham_adi[0, 0])
    print("hvib_adi:\n", storage.hvib_adi[0, 0])

    section("Active matrix")
    active_hvib = engine.active_matrix(
        storage,
        traj,
        rep="adiabatic",
        kind="vibronic",
    )
    print("active Hvib batch shape:", active_hvib.shape)


if __name__ == "__main__":
    main()
