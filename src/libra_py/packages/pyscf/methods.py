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
from libra_py.packages.pyscf.interfaces import (
    ES_Request,
    ES_Result,
    ES_Strategy,
    MolecularGeometry,
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


def _params_to_request(params):
    nstates = params.get("nstates", 2)

    return ES_Request(
        n_singlets=nstates,
        n_triplet=0,
        H_soc=params.get("H_soc", False),
        gradient_state=params.get("gradient_state", "all"),
        hessian_state=params.get("hessian_state"),
        nacv=params.get("nacv", False),
        time_overlap=params.get("time_overlap", True),
    )


# =============================================================================
# ES calculation helper
# =============================================================================

def _compute_es_result(
    strategy: ES_Strategy,
    geometry: MolecularGeometry,
    request: ES_Request,
    previous: ES_Strategy | None,
):
    result = ES_Result()

    strategy.compute_result(
        geometry,
        request,
        result,
        previous=previous,
    )

    return result


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
    obj.basis_transform = CMATRIX(
        nstates,
        nstates,
    )

    for i in range(nstates):
        obj.basis_transform.set(
            i,
            i,
            1.0 + 0.0j,
        )

    # Gradients
    if result.gradients is not None:

        obj.d1ham_adi = _gradients_to_d1ham_adi(
            result.gradients,
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

        obj.time_overlap_adi = CMATRIX(
            nstates,
            nstates,
        )

        obj.hvib_adi = _energies_to_ham_adi(
            result.H_el
        )

    return obj


# =============================================================================
# Main Libra callback
# =============================================================================

def strategy_compute_adi(
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
    request = _params_to_request(params)

    # 2. Get current and previous ES snapshots
    previous = params.setdefault(
        "es_previous",
        {},
    ).get(itraj)

    strategy_spec = params.get(
        "strategy_factory",
        params.get("es_strategy"),
    )

    if strategy_spec is None:
        raise KeyError("Missing strategy specification: expected 'strategy_factory' or 'es_strategy'.")

    if callable(strategy_spec):
        current = strategy_spec()
    elif hasattr(strategy_spec, "copy"):
        current = strategy_spec.copy()
    elif hasattr(strategy_spec, "clone"):
        current = strategy_spec.clone()
    else:
        current = strategy_spec

    # 3. Run ES calculation
    result = _compute_es_result(
        current,
        geometry,
        request,
        previous,
    )

    # 4. Generic ES result -> Libra result
    obj = _es_result_to_libra(
        result,
        request,
        natoms=len(params["atom_labels"]),
        dt=float(params.get("dt", 41.0)),
    )

    # 5. Current becomes previous
    params["es_previous"][itraj] = current

    return obj
