"""
Demonstrate TensorView usage.

Run from this directory:

    python test3.py
"""

import numpy as np
from pathlib import Path
import sys

SRC_DIR = Path(__file__).resolve().parents[3]
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from libra_py.dyn.core import storage as st
from libra_py.dyn.core.views import TensorView


def show(name, value):
    print(f"{name:34s}: {value}")


def main():
    # TensorStorage owns the actual arrays. TensorView only points to one
    # (trajectory, TBF) slot inside those arrays.
    storage = st.TensorStorage(
        backend=np,
        ntraj=2,
        ndof=3,
        nstates=2,
        ntbf_initial=1,
        ntbf_capacity=3,
    )

    view = TensorView(storage, traj_id=1, tbf_id=0)

    show("view indices", (view.traj_id, view.tbf_id))
    show("storage q shape", storage.q.shape)
    show("view q shape", view.q.shape)

    # Properties such as view.q and view.ampl_adi are writable slices.
    # Updating the view updates the shared TensorStorage array.
    view.q[:] = [0.1, 0.2, 0.3]
    view.p[:] = [1.0, 0.0, -1.0]
    view.iM[:] = [1.0, 0.5, 0.25]
    view.ampl_adi[:] = [1.0 + 0.0j, 0.0 + 0.0j]
    view.ampl_dia[:] = [0.2 + 0.0j, 0.8 + 0.0j]
    view.set("act_states", 0)

    show("view q", view.q.tolist())
    show("storage q[1,0]", storage.q[1, 0].tolist())
    show("view ampl_adi", view.ampl_adi.tolist())
    show("view ampl_dia", view.ampl_dia.tolist())

    # The same applies to matrix-valued electronic/Hamiltonian data.
    view.dm_adi[:] = np.outer(view.ampl_adi, view.ampl_adi.conjugate())
    view.ham_adi[:] = np.diag([0.0, 0.25])
    view.nac_adi[0, 1] = 0.01
    view.hvib_adi[:] = view.ham_adi - 1j * view.nac_adi
    view.time_overlap_adi[:] = np.eye(storage.nstates)

    show("dm_adi", view.dm_adi.tolist())
    show("hvib_adi", view.hvib_adi.tolist())

    # Generic access lets code work with field names dynamically.
    view.set("f", [-0.5, 0.0, 0.5])
    force = view.get("f")
    show("force via get/set", force.tolist())

    # Grouped helpers return dictionaries of writable slices.
    nuclear = view.nuclear()
    electronic = view.electronic()
    hamiltonian = view.hamiltonian()

    nuclear["q"][0] += 1.0
    electronic["q_mm"][:] = np.sqrt(2.0) * np.real(view.ampl_adi)
    hamiltonian["cum_phase_corr"][:] = np.eye(storage.nstates)

    show("q after nuclear dict edit", view.q.tolist())
    show("q_mm from electronic dict", view.q_mm.tolist())
    show("cum_phase_corr", view.cum_phase_corr.tolist())

    # Optional fields are not allocated until the matching storage allocator
    # is called. maybe() is useful for probing availability.
    show("dR before allocate_afssh", view.maybe("dR"))

    try:
        _ = view.dR
    except AttributeError as err:
        show("direct dR access", str(err))

    storage.allocate_afssh()
    view.dR[0, 0, 1] = 0.123
    show("dR after allocate_afssh", view.dR.tolist())

    storage.allocate_hamiltonian_derivatives(der_lvl=1)
    view.dc1_adi[0, 0, 1] = 0.04
    view.d1ham_adi[2, 1, 1] = -0.7
    show("dc1_adi shape", view.dc1_adi.shape)
    show("d1ham_adi[2]", view.d1ham_adi[2].tolist())
    show("d2ham_adi before der_lvl=2", view.maybe("d2ham_adi"))

    storage.allocate_hamiltonian_derivatives(der_lvl=2)
    view.d2ham_adi[0, 1, 0, 0] = 2.5
    show("d2ham_adi shape", view.d2ham_adi.shape)

    storage.allocate_simple_decoherence()
    view.coherence_clocks[0, 1] = 12.0
    view.gaps_curr[0, 1] = 0.33
    show("coherence clocks", view.coherence_clocks.tolist())
    show("current gaps", view.gaps_curr.tolist())

    # optional() collects only the optional fields that have been allocated.
    optional_fields = view.optional()
    show("allocated optional keys", sorted(optional_fields.keys()))

    assert storage.q[1, 0, 0] == view.q[0]
    assert storage.dR[1, 0, 0, 0, 1] == 0.123
    assert storage.coherence_clocks[1, 0, 0, 1] == 12.0

    print("TensorView demonstration completed successfully.")


if __name__ == "__main__":
    main()
