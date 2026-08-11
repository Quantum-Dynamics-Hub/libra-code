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
.. module:: methods
   :platform: Unix, Windows
   :synopsis: This module implements an adapter between the ABC interface and the compute_adi 

.. moduleauthor::
       Alexey V. Akimov, Jieyang Gu

"""

import numpy as np

from liblibra_core import CMATRIX, CMATRIXList, Cpp2Py

from libra_py import units
from libra_py.packages.pyscf.interfaces import (
    ES_Request,
    ES_Result,
    ES_Strategy,
    MolecularGeometry,
)


def _matrix2nparray(matrix, dtype=float):
    return np.array(
        [
            [matrix.get(i, j) for j in range(matrix.num_of_cols)]
            for i in range(matrix.num_of_rows)
        ],
        dtype=dtype,
    )


def _q_to_geometry(q, itraj, atom_labels):
    coords = q.col(itraj)

    coordinates = (
        _matrix2nparray(coords, float).reshape(-1, 3) / units.Angst
    )

    return MolecularGeometry(
        atom_labels=tuple(atom_labels),
        coords_bohr=coordinates,
    )


class tmp:
    pass

# =============================================================================
# Input conversion helpers
# =============================================================================

def _get_trajectory_index(full_id):
    if isinstance(full_id, (list, tuple)):
        return int(full_id[-1])
    Id = Cpp2Py(full_id)
    return int(Id[-1])


def _params_to_request(params, itraj):
    nstates = int(params.get("nstates", 2))

    return ES_Request(
        n_singlets=nstates,
        n_triplets=int(params.get("n_triplets", params.get("n_triplet", 0))),
        H_soc=params.get("H_soc", False),
        gradient_state=_resolve_gradient_state(params, itraj),
        hessian_state=params.get("hessian_state"),
        nacv=params.get("nacv", False),
        time_overlap=params.get("time_overlap", True),
    )


def _resolve_gradient_state(params, itraj=None):
    gradient_state = params.get(
        "gradient",
        params.get("gradient_mode", params.get("gradient_state", "all")),
    )

    if gradient_state is None or gradient_state is False:
        return None

    if gradient_state is True:
        return "all"

    if isinstance(gradient_state, str):
        gradient_state = gradient_state.lower()

        if gradient_state == "none":
            return None

        if gradient_state == "all":
            return "all"

        if gradient_state != "active":
            raise ValueError(
                "Gradient mode must be one of 'active', 'all', or 'none'."
            )

        if itraj is None:
            raise ValueError(
                "Active-state gradient requests require a trajectory index."
            )

        act_states = params.get("act_state", params.get("active_states"))

        if isinstance(act_states, dict):
            act_state = act_states.get(itraj)
        elif isinstance(act_states, (list, tuple)):
            act_state = act_states[itraj]
        else:
            act_state = act_states

        if act_state is None:
            raise KeyError(
                f"Missing active state for trajectory {itraj} in params['act_state']."
            )

        return int(act_state)

    return int(gradient_state)


def _get_strategy(params, itraj):
    strategies = params.setdefault("es_strategies", {})

    if itraj in strategies:
        return strategies[itraj]

    strategy_spec = params.get(
        "strategy_factory",
        params.get("es_strategy"),
    )

    if strategy_spec is None:
        raise KeyError(
            "Missing strategy specification: expected 'strategy_factory' or 'es_strategy'."
        )

    if isinstance(strategy_spec, dict):
        strategy = strategy_spec[itraj]
    elif callable(strategy_spec):
        strategy = strategy_spec()
    elif hasattr(strategy_spec, "copy"):
        strategy = strategy_spec.copy()
    elif hasattr(strategy_spec, "clone"):
        strategy = strategy_spec.clone()
    else:
        strategy = strategy_spec

    strategies[itraj] = strategy
    return strategy


# =============================================================================
# ES calculation helper
# =============================================================================

def _compute_es_result(
    strategy: ES_Strategy,
    geometry: MolecularGeometry,
    request: ES_Request,
):
    return strategy.compute_result(geometry, request)


# =============================================================================
# ES_Result -> Libra type conversion helpers
# =============================================================================

def _numpy_to_cmatrix(array):
    array = np.asarray(array)

    nrows, ncols = array.shape
    matrix = CMATRIX(nrows, ncols)

    for i in range(nrows):
        for j in range(ncols):
            matrix.set(
                i,
                j,
                complex(array[i, j]),
            )

    return matrix


def _energies_to_ham_adi(energies):
    return _numpy_to_cmatrix(
        np.diag(energies)
    )


def _identity_cmatrix(size):
    matrix = CMATRIX(size, size)
    for i in range(size):
        matrix.set(i, i, 1.0 + 0.0j)
    return matrix


def _zero_cmatrix_list(nstates, natoms):
    matrices = CMATRIXList()
    for _ in range(3 * natoms):
        matrices.append(CMATRIX(nstates, nstates))
    return matrices


def _gradients_to_d1ham_adi(
    gradients,
    nstates,
    natoms,
):
    d1ham_adi = CMATRIXList()

    for dof in range(3 * natoms):
        d1ham_adi.append(
            CMATRIX(nstates, nstates)
        )

    for state, gradient in enumerate(gradients):

        if gradient is None:
            continue

        for atom in range(natoms):

            for xyz in range(3):

                dof = 3 * atom + xyz

                d1ham_adi[dof].set(
                    state,
                    state,
                    complex(
                        gradient[atom, xyz]
                    ),
                )

    return d1ham_adi


def _nacv_to_dc1adi(nac_vectors, nstates, natoms):
    """Convert nac_vectors (n_total,n_total,natoms,3) -> CMATRIXList of derivative couplings per DOF."""
    dc1adi = CMATRIXList()
    for _ in range(3 * natoms):
        dc1adi.append(CMATRIX(nstates, nstates))

    # nac_vectors[i,j,atom,xyz] in Bohr^-1; place into per-dof CMATRIXList
    for i in range(nstates):
        for j in range(nstates):
            nv = nac_vectors[i, j]
            if nv is None:
                continue
            for atom in range(natoms):
                for xyz in range(3):
                    dof = 3 * atom + xyz
                    dc1adi[dof].set(i, j, complex(nv[atom, xyz]))

    return dc1adi


def _time_overlap_to_cmatrix(time_overlap):
    return _numpy_to_cmatrix(
        time_overlap
    )


def _time_overlap_to_hvib(
    energies,
    time_overlap,
    dt,
):
    nstates = len(energies)

    hvib = _energies_to_ham_adi(
        energies
    )

    for i in range(nstates):

        for j in range(i + 1, nstates):

            dij = (
                time_overlap[i, j]
                - time_overlap[j, i]
            ) / (2.0 * dt)

            hvib.set(
                i,
                j,
                -1j * dij,
            )

            hvib.set(
                j,
                i,
                +1j * dij,
            )

    return hvib


# =============================================================================
# Build Libra callback result
# =============================================================================

def _es_result_to_libra(
    result: ES_Result,
    request: ES_Request,
    natoms: int,
    dt: float,
) -> tmp:
    nstates = request.n_total

    obj = tmp()

    # Energies
    obj.ham_adi = _energies_to_ham_adi(
        result.H_el
    )

    # Identity adiabatic transformation
    obj.basis_transform = _identity_cmatrix(nstates)
    obj.ovlp_adi = _identity_cmatrix(nstates)
    obj.dc1_adi = _zero_cmatrix_list(nstates, natoms)
    obj.d1ham_adi = _zero_cmatrix_list(nstates, natoms)

    # Gradients
    if result.gradients is not None:

        obj.d1ham_adi = _gradients_to_d1ham_adi(
            result.gradients,
            nstates,
            natoms,
        )

    # NAC vectors (nonadiabatic coupling vectors)
    if result.nac_vectors is not None:
        obj.dc1_adi = _nacv_to_dc1adi(
            result.nac_vectors,
            nstates,
            natoms,
        )

    # Time overlaps
    if result.time_overlap is not None:

        obj.time_overlap_adi = (
            _time_overlap_to_cmatrix(
                result.time_overlap
            )
        )

        obj.hvib_adi = (
            _time_overlap_to_hvib(
                result.H_el,
                result.time_overlap,
                dt,
            )
        )

    else:

        obj.time_overlap_adi = _identity_cmatrix(nstates)

        obj.hvib_adi = _energies_to_ham_adi(
            result.H_el
        )

    return obj


# =============================================================================
# Main Libra callback
# =============================================================================

def pyscf_compute_adi(
    q,
    params,
    full_id,
):
    # 1. Libra input -> generic ES input
    itraj = _get_trajectory_index(full_id)
    geometry = _q_to_geometry(
        q,
        itraj,
        params["atom_labels"],
    )
    request = _params_to_request(params, itraj)

    # 2. Get the stateful ES strategy for this trajectory.
    current = _get_strategy(params, itraj)

    # 3. Run ES calculation
    result = _compute_es_result(
        current,
        geometry,
        request,
    )

    # 4. Generic ES result -> Libra result
    obj = _es_result_to_libra(
        result,
        request,
        natoms=len(params["atom_labels"]),
        dt=float(params.get("dt", 41.0)),
    )

    return obj


#
# Backward-compatible names and convenience wrapper
#


def compute_model(q, params, full_id):
    """Libra compute_model callback (backwards-compatible signature).

    Accepts `params["gradient"]` with values "active", "all", or "none".
    Accepts `params["nacv"]` boolean to request NAC vectors (nonadiabatic coupling vectors).
    """
    return pyscf_compute_adi(q, params, full_id)



# Backward-compatible name used by older PySCF examples.
strategy_compute_adi = pyscf_compute_adi


def main() -> None:
    """Mini demo: two-trajectory SA-3 CASSCF STO-3G HeH+ FSSH run (warmup + short run)."""
    from liblibra_core import Random, dyn_variables, nHamiltonian, update_Hamiltonian_variables
    from libra_py.dynamics.tsh.compute import run_dynamics
    from libra_py.dynamics.tsh.recipes.fssh_h_plus import load as load_fssh_recipe
    from libra_py.packages.pyscf.implementations.casscf import CASSCF

    # Parameters (kept minimal and consistent with examples)
    NTRAJ = 2
    NSTATES = 3
    NSTEPS = 2
    DT_FS = 0.5
    ATOM_LABELS = ["He", "H"]

    def build_strategy_factory(nstates):
        def factory():
            return CASSCF(norbcas=2, nelecas=2, nroots=nstates, basis="sto-3g", charge=1, unit="Bohr")

        return factory

    dt_au = DT_FS * units.fs2au

    compute_model_fn = pyscf_compute_adi
    model_params = {
        "model": 0,
        "model0": 0,
        "nstates": NSTATES,
        "atom_labels": ATOM_LABELS,
        "strategy_factory": build_strategy_factory(NSTATES),
        "gradient": "active",   # "active"|"all"|"none"
        "nacv": True,
        "time_overlap": True,
        "dt": dt_au,
        "act_state": {itraj: 1 for itraj in range(NTRAJ)},
    }

    dyn_params = {}
    load_fssh_recipe(dyn_params)
    dyn_params.update({"ntraj": NTRAJ, "nsteps": NSTEPS, "dt": dt_au, "prefix": "hehp_sa3_casscf_sto3g"})

    # Minimal initial nuclear/electronic structures (mirrors examples)
    init_nucl = {"init_type": 0, "ndof": 6, "q": [0.0]*6, "p": [0.0]*6, "mass": [4.002602*units.amu]*3 + [1.007825*units.amu]*3}
    init_elec = {"init_type": 3, "ndia": NSTATES, "nadi": NSTATES, "nstates": NSTATES, "istates": [0.0]*NSTATES, "rep": 1, "ntraj": NTRAJ}
    init_elec["istates"][1] = 1.0

    rnd = Random()
    dyn_var = dyn_variables(NSTATES, NSTATES, 6, NTRAJ)
    dyn_var.init_nuclear_dyn_var(init_nucl, rnd)
    dyn_var.init_amplitudes(init_elec, rnd)
    dyn_var.init_density_matrix(init_elec)
    dyn_var.init_auxiliary_variables(init_elec, rnd)

    ham = nHamiltonian(NSTATES, NSTATES, 6)
    ham.add_new_children(NSTATES, NSTATES, 6, NTRAJ)
    ham.init_all(2, 1)

    warmup_params = dict(model_params)
    warmup_params["timestep"] = 0
    update_Hamiltonian_variables(dyn_params, dyn_var, ham, ham, compute_model_fn, warmup_params, 0)
    update_Hamiltonian_variables(dyn_params, dyn_var, ham, ham, compute_model_fn, warmup_params, 1)

    dyn_var.update_basis_transform(ham)
    dyn_var.update_amplitudes({"rep_tdse": init_elec["rep"]}, ham)
    dyn_var.update_density_matrix(dyn_params, ham, 1)
    dyn_var.init_active_states(init_elec, rnd)

    run_dynamics(dyn_var, dyn_params, ham, compute_model_fn, model_params, rnd)

    print("Mini HeH+ SA-3 CASSCF STO-3G FSSH run completed.")


if __name__ == "__main__":
    main()
