# *********************************************************************************
# * Copyright (C) 2026  Jieyang Gu and Alexey V. Akimov
# *
# * This file is distributed under the terms of the GNU General Public License
# * as published by the Free Software Foundation, either version 3 of
# * the License, or (at your option) any later version.
# * See the file LICENSE in the root directory of this distribution
# * or <http://www.gnu.org/licenses/>.
# *
# *********************************************************************************/
"""
.. module:: debug
   :platform: Unix, Windows
   :synopsis: Per-step console dump of the FSSH electronic state and the
       nHamiltonian fields, for interactively debugging a NA-MD run.

.. moduleauthor::
       Jieyang Gu

"""

import numpy as np
from liblibra_core import intList


def _cmatrix_to_str(m, fmt="{:+.4e}"):
    rows = []
    for i in range(m.num_of_rows):
        row = "  ".join(
            f"{fmt.format(m.get(i, j).real)}{fmt.format(m.get(i, j).imag)}j"
            for j in range(m.num_of_cols)
        )
        rows.append(f"    [{row}]")
    return "\n".join(rows)


def print_step_debug(i, dyn_var, ham, dyn_params, traj=None):
    """Print the C vector and the nHam fields relevant to FSSH for one step.

    Call this right after ``compute_dynamics(...)`` returns inside the main
    loop of ``libra_py.dynamics.tsh.compute.run_dynamics`` -- at that point
    ``dyn_var`` and ``ham`` both hold the state produced by *this* step, since
    ``compute_dynamics`` mutates them in place rather than returning a copy.

    Parameters
    ----------
    i : int
        Current step index (as used in ``for i in range(nsteps)``).

    dyn_var : dyn_variables
        Holds the electronic amplitudes and the active-state indices.

    ham : nHamiltonian
        The top-level Hamiltonian object. Per-trajectory quantities are
        fetched through the ``[0, itraj]`` id path -- the same convention
        used for ``full_id`` in the ``compute_model`` callback -- because
        ``ham.children`` is not exposed to Python.

    dyn_params : dict
        Read only for ``rep_tdse`` (which amplitude representation is
        dynamically consistent) and ``quantum_dofs`` (which dof's
        ``dc1_adi`` to print).

    traj : int, list of int, or None, optional
        Which trajectories to print. ``None`` (default) prints all of them.
        In NBRA runs only one Hamiltonian child exists regardless of
        ``ntraj``, so its per-trajectory index is clamped to 0.
    """
    ntraj = dyn_var.ntraj
    nadi = dyn_var.nadi
    ndia = dyn_var.ndia
    rep_tdse = dyn_params.get("rep_tdse", 1)
    quantum_dofs = dyn_params.get("quantum_dofs") or [0]

    if traj is None:
        traj_list = list(range(ntraj))
    elif isinstance(traj, int):
        traj_list = [traj]
    else:
        traj_list = list(traj)

    Cadi = dyn_var.get_ampl_adi()
    Cdia = dyn_var.get_ampl_dia()

    print(f"\n########## NAMD step {i} ##########")

    # NBRA runs allocate a single Hamiltonian child shared by every trajectory
    # (generic_recipe.py: ham.add_new_children(ndia, nadi, ndof, 1) when
    # isNBRA == 1), so the per-trajectory id must be clamped to 0 in that case.
    #
    # Read dyn_params["isNBRA"], NOT dyn_params["is_nbra"]: generic_recipe's
    # own comn.check_input call (compute.py:1163) unconditionally seeds
    # dyn_params["is_nbra"] = 0 into the SAME dict before run_dynamics ever
    # runs, so that lowercase key reads back as 0 regardless of what isNBRA
    # was set to. ham's actual child count was decided from isNBRA, before
    # run_dynamics was even called (compute.py:1201) -- so isNBRA is the only
    # value that reflects how many children ham really has.
    is_nbra = dyn_params.get("isNBRA", 0)

    for itraj in traj_list:
        # The C++ get_*_adi(id_) overloads take vector<int> by reference; a
        # plain Python list does not auto-convert (Boost.Python raises
        # ArgumentError), so the id path has to be built as an intList.
        ham_id = intList()
        ham_id.append(0)
        ham_id.append(0 if is_nbra == 1 else itraj)

        print(f"--- traj {itraj}  (active state = {dyn_var.act_states[itraj]}) ---")

        print("  C_adi (adiabatic amplitudes, this trajectory's column):")
        for st in range(nadi):
            c = Cadi.get(st, itraj)
            marker = "  <- active" if st == dyn_var.act_states[itraj] else ""
            print(f"    state {st}: {c.real:+.6e} {c.imag:+.6e}j  "
                  f"|C|^2={abs(c)**2:.6f}{marker}")

        print("  C_dia (diabatic amplitudes, this trajectory's column):")
        for st in range(ndia):
            c = Cdia.get(st, itraj)
            print(f"    state {st}: {c.real:+.6e} {c.imag:+.6e}j  "
                  f"|C|^2={abs(c)**2:.6f}")

        print(f"  rep_tdse = {rep_tdse}  "
              f"({'adiabatic' if rep_tdse == 1 else 'diabatic'} amplitudes are the dynamically consistent ones)")

        ham_adi = ham.get_ham_adi(ham_id)
        hvib_adi = ham.get_hvib_adi(ham_id)
        time_overlap_adi = ham.get_time_overlap_adi(ham_id)
        basis_transform = ham.get_basis_transform(ham_id)

        print("  ham_adi (adiabatic energies, diagonal):")
        print("    " + "  ".join(f"{ham_adi.get(k, k).real:+.6e}" for k in range(nadi)))

        print("  hvib_adi = Ham - i*hbar*NAC:")
        print(_cmatrix_to_str(hvib_adi))

        print("  nac_adi (derivative coupling <i|d/dt|j>):")
        try:
            # nHamiltonian.get_nac_adi(id_) is unreachable from Python: the
            # boost.python export table registers get_ham_adi's v2 overload
            # twice (src/nhamiltonian/libnhamiltonian.cpp:482) instead of
            # get_nac_adi's, so only the zero-arg, level-0 get_nac_adi() is
            # bound. That overload does not carry per-trajectory data, so it
            # raises Boost.Python.ArgumentError here. hvib_adi above still
            # carries the NAC information (Hvib = Ham - i*hbar*NAC), so
            # nothing is actually lost -- this is cosmetic until the binding
            # table is fixed and rebuilt.
            print(_cmatrix_to_str(ham.get_nac_adi(ham_id)))
        except Exception as e:
            print(f"    <unavailable: {type(e).__name__} -- "
                  f"get_nac_adi(id_) is not exposed to Python, see "
                  f"libnhamiltonian.cpp:482; recover it from hvib_adi instead>")

        for dof in quantum_dofs:
            dc1 = ham.get_dc1_adi(dof, ham_id)
            print(f"  dc1_adi[dof={dof}] (derivative coupling vector component):")
            print(_cmatrix_to_str(dc1))

        print("  time_overlap_adi <psi(t)|psi(t+dt)>:")
        print(_cmatrix_to_str(time_overlap_adi))

        print("  basis_transform (dia->adi eigenvectors):")
        print(_cmatrix_to_str(basis_transform))
