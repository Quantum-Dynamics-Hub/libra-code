"""
Educational example: TensorView.

TensorView is a writable view into TensorStorage for one pair:

    (trajectory id, TBF storage slot)

It does not copy data. Assigning through a view changes the shared storage.

Run from the repository root:

    PYTHONPATH=src python examples/libra_py/dyn/core/example_tensor_view.py
"""

import numpy as np

from libra_py.dyn.core.storage import TensorStorage
from libra_py.dyn.core.views import TensorView


def section(title):
    print(f"\n--- {title} ---")


def main():
    storage = TensorStorage(
        backend=np,
        ntraj=2,
        ndof=3,
        nstates=2,
        ntbf_initial=1,
    )
    view = TensorView(storage, traj_id=1, tbf_id=0)

    section("Writable slices")

    view.q[:] = [0.1, 0.2, 0.3]
    view.p[:] = [1.0, 0.0, -1.0]
    view.ampl_adi[:] = [1.0 + 0.0j, 0.0 + 0.0j]
    view.ham_adi[:] = np.diag([0.0, 0.5])

    print("view q:", view.q)
    print("same storage data:", storage.q[1, 0])
    print("view ham_adi:\n", view.ham_adi)

    section("Generic access")

    # set/get are useful when field names are chosen dynamically.
    view.set("act_states", 0)
    print("active state:", view.get("act_states"))

    # Group helpers return dictionaries of writable slices.
    nuclear = view.nuclear()
    electronic = view.electronic()
    hamiltonian = view.hamiltonian()

    nuclear["q"][0] += 1.0
    electronic["q_mm"][:] = np.sqrt(2.0) * np.real(view.ampl_adi)
    hamiltonian["time_overlap_adi"][:] = np.eye(storage.nstates)

    print("q after grouped edit:", view.q)
    print("q_mm:", view.q_mm)
    print("time_overlap_adi:\n", view.time_overlap_adi)

    section("Optional fields")

    # maybe() returns None for unallocated optional fields.
    print("dR before allocation:", view.maybe("dR"))

    storage.allocate_afssh()
    view.dR[0, 0, 1] = 0.123
    print("dR after allocation:\n", view.dR)
    print("allocated optional fields:", sorted(view.optional().keys()))


if __name__ == "__main__":
    main()
