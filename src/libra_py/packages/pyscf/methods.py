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
from types import SimpleNamespace

from liblibra_core import CMATRIX, CMATRIXList, Cpp2Py

from libra_py import data_conv
from libra_py.packages.pyscf.interfaces import (
    ES_Request,
    ES_Result,
    ES_Strategy,
    MolecularGeometry,
)

def _q_to_geometry(q, itraj, atom_labels):
    """Libra stores nuclear coordinates in Bohr; keep them as ``coords_bohr``."""
    coords = q.col(itraj)
    coordinates = data_conv.MATRIX2nparray(coords, np.float64).reshape(-1, 3)
    return MolecularGeometry( atom_labels=tuple(atom_labels), coords_bohr=coordinates )

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

def _resolve_gradient_state(params, itraj):
    """Translate ``params["gradient_state"]`` into what ``ES_Request`` expects.

    Accepted values:
        None    -> no gradient
        int     -> gradient for that one state
        "all"   -> gradients for every state
        "active"-> gradient for the trajectory's current active state, read
                   from ``params["act_state"][itraj]`` (written every step by
                   the dynamics driver, e.g. libra_py.dynamics.tsh.compute).

    ``act_state`` is absent on the warm-up calls that precede the propagation
    loop (libra_py.dynamics.tsh.compute calls the model at lines 1018-1019 but
    only fills act_state at line 1082), so a missing entry falls back to the
    ground state rather than raising -- the same convention the mopac, dftbplus,
    and cp2k interfaces use.
    """
    gradient_state = params.get("gradient_state", "all")

    if isinstance(gradient_state, str):
        gradient_state = gradient_state.lower()

    if gradient_state in (None, "all"):
        return gradient_state

    if gradient_state == "active":
        return params.setdefault("act_state", {}).get(itraj, 0)

    return int(gradient_state)

def _params_to_request(params, itraj): #possible fields of params: "nstates", "H_soc",
                                       #gradient_state", "hes", "nacv", "time_overlap"
    nstates = params.get("nstates", 2)

    return ES_Request(
        n_singlets=nstates,
        n_triplets=0,
        H_soc=params.get("H_soc", False),
        gradient_state=_resolve_gradient_state(params, itraj),
        hessian_state=params.get("hessian_state"),
        nacv=params.get("nacv", False),
        time_overlap=params.get("time_overlap", True),
    )

# =============================================================================
# ES_Result -> Libra type conversion helpers
# =============================================================================

def _energies_to_ham_adi(energies):
    return data_conv.nparray2CMATRIX( np.diag(energies) )

def _gradients_to_d1ham_adi(gradients, nstates, natoms):
    """Pack ES_Result.gradients into Libra's ``d1ham_adi``.

    ``gradients`` is one ``(natoms, 3)`` block per electronic state, or None for
    states whose gradient was not requested. Libra wants one
    CMATRIX(nstates, nstates) per nuclear dof, with ``dof = 3 * atom + xyz``,
    carrying dE_state/dq on its diagonal.
    """
    per_dof = np.zeros((3 * natoms, nstates, nstates), dtype=np.complex128)

    for state, gradient in enumerate(gradients):

        if gradient is None:
            continue

        per_dof[:, state, state] = np.asarray(gradient).reshape(-1)

    d1ham_adi = CMATRIXList()

    for dof in range(3 * natoms):
        d1ham_adi.append( data_conv.nparray2CMATRIX(per_dof[dof]) )

    return d1ham_adi

def _time_overlap_to_cmatrix(time_overlap):
    return data_conv.nparray2CMATRIX( np.asarray(time_overlap) )

def _nac_vectors_to_dc1_adi(nac_vectors, nstates, natoms ):
    """Pack ES_Result.nac_vectors into Libra's ``dc1_adi``.

    ``nac_vectors`` has shape ``(nstates, nstates, natoms, 3)``; Libra wants one
    CMATRIX(nstates, nstates) per nuclear dof, with ``dof = 3 * atom + xyz``.
    """
    # (nstates, nstates, natoms, 3) -> (3 * natoms, nstates, nstates)
    per_dof = np.asarray(nac_vectors).reshape(nstates, nstates, 3 * natoms)
    per_dof = np.moveaxis(per_dof, -1, 0)

    dc1_adi = CMATRIXList()

    for dof in range(3 * natoms):
        dc1_adi.append( data_conv.nparray2CMATRIX(per_dof[dof]) )

    return dc1_adi

def _time_overlap_to_hvib(  energies, time_overlap, dt):

    nstates = len(energies)
    hvib = _energies_to_ham_adi( energies  )

    for i in range(nstates):
        for j in range(i + 1, nstates):
            dij = ( time_overlap[i, j] - time_overlap[j, i] ) / (2.0 * dt)
            hvib.set( i, j, -1j * dij )
            hvib.set( j, i, +1j * dij )

    return hvib

# =============================================================================
# Build Libra callback result
# =============================================================================

def _es_result_to_libra(result: ES_Result, request: ES_Request, natoms: int, dt: float) -> tmp:

    nstates = request.n_total

    obj = tmp()

    # Energies
    obj.ham_adi = _energies_to_ham_adi( result.H_el )

    # Identity adiabatic transformation
    obj.basis_transform = CMATRIX( nstates, nstates )

    for i in range(nstates):
        obj.basis_transform.set( i, i, 1.0 + 0.0j )

    # Gradients
    if result.gradients is not None:
        obj.d1ham_adi = _gradients_to_d1ham_adi( result.gradients, nstates, natoms )

    # Derivative couplings
    if result.nac_vectors is not None:
        obj.dc1_adi = _nac_vectors_to_dc1_adi( result.nac_vectors, nstates, natoms )

    # Time overlaps
    if result.time_overlap is not None:
        obj.time_overlap_adi = _time_overlap_to_cmatrix( result.time_overlap )
        obj.hvib_adi = _time_overlap_to_hvib( result.H_el, result.time_overlap, dt )

    elif request.time_overlap:
        # Overlaps were asked for, but this is the first geometry of the
        # trajectory: there is no previous state to overlap against
        # (ES_Strategy.compute_result only computes one when
        # get_previous_state() is not None).
        #
        # S(t0, t0) is the identity -- no time has elapsed, so every state
        # overlaps perfectly with itself and no coupling has accumulated yet.
        # A zero matrix would instead assert that consecutive electronic states
        # are mutually orthogonal, which is both physically false and rejected
        # outright by nac_npi: validate_npi_input requires |S^T S - I| < 1e-6
        # (src/calculators/NPI.cpp:161), so nac_update_method=2 with
        # nac_algo=NPI would raise on the very first step.
        obj.time_overlap_adi = CMATRIX( nstates, nstates )

        for i in range(nstates):
            obj.time_overlap_adi.set( i, i, 1.0 + 0.0j )

        obj.hvib_adi = _energies_to_ham_adi( result.H_el )

    # Otherwise time-overlaps were never requested, so no attribute is set at
    # all -- the same convention used above for d1ham_adi and dc1_adi.
    # nHamiltonian copies only the attributes that are present (the hasattr
    # gate in nHamiltonian_compute_adiabatic.cpp:544-680), so leaving it off
    # means "this callback has nothing to say about time-overlaps" rather than
    # overwriting whatever the Hamiltonian already holds.

    return obj


def _copy_strategy(strategy_spec):
    """Create an independent electronic-structure strategy from a specification."""
    if callable(strategy_spec):
        return strategy_spec()
    if hasattr(strategy_spec, "copy"):
        return strategy_spec.copy()
    if hasattr(strategy_spec, "clone"):
        return strategy_spec.clone()
    return strategy_spec


def _manifold_ms_labels(nroots, multiplicity, include_ms=True):
    """Return ``(root, 2*Ms)`` labels for one spin-free manifold."""
    if multiplicity < 1:
        raise ValueError("Spin multiplicities must be positive integers")
    if include_ms:
        return [
            (root, multiplicity - 1 - 2 * projection)
            for root in range(nroots)
            for projection in range(multiplicity)
        ]
    return [(root, multiplicity - 1) for root in range(nroots)]


def _compute_spin_manifold(geometry, params, itraj, manifold, index):
    """Compute one PySCF spin manifold and expand it over requested Ms values."""
    multiplicity = int(manifold["spin"])
    nroots = int(manifold.get("nroots", 1))
    include_ms = bool(manifold.get("include_ms_projections", True))
    labels = _manifold_ms_labels(nroots, multiplicity, include_ms)

    strategy_spec = manifold.get("strategy_factory", manifold.get("es_strategy"))
    if strategy_spec is None:
        raise KeyError(
            f"spin_manifolds[{index}] requires 'strategy_factory' or 'es_strategy'"
        )

    all_strategies = params.setdefault("es_spin_strategies", {})
    trajectory_strategies = all_strategies.setdefault(itraj, {})
    if index not in trajectory_strategies:
        trajectory_strategies[index] = _copy_strategy(strategy_spec)
    strategy = trajectory_strategies[index]

    if int(strategy.nroots) != nroots:
        raise ValueError(
            f"spin_manifolds[{index}] requests {nroots} roots, but its strategy "
            f"provides {strategy.nroots}"
        )
    strategy_multiplicity = getattr(strategy, "spin_multiplicity", multiplicity)
    if int(strategy_multiplicity) != multiplicity:
        raise ValueError(
            f"spin_manifolds[{index}] has multiplicity {multiplicity}, but its "
            f"strategy provides {strategy_multiplicity}"
        )

    request = ES_Request(
        n_singlets=nroots,
        n_triplets=0,
        H_soc=False,
        gradient_state=manifold.get(
            "gradient_state", params.get("gradient_state", "all")
        ),
        hessian_state=manifold.get(
            "hessian_state", params.get("hessian_state")
        ),
        nacv=bool(manifold.get("nacv", params.get("nacv", False))),
        time_overlap=bool(
            manifold.get("time_overlap", params.get("time_overlap", True))
        ),
    )
    result = strategy.compute_result(geometry, request)

    root_of_state = [root for root, _ in labels]
    ms_of_state = [twice_ms for _, twice_ms in labels]
    nstates = len(labels)
    energies = np.asarray(result.H_el, dtype=float)
    if energies.shape != (nroots,):
        raise ValueError(
            f"spin_manifolds[{index}] returned energies with shape "
            f"{energies.shape}; expected ({nroots},)"
        )

    time_overlap = np.eye(nstates)
    if result.time_overlap is not None:
        root_overlap = np.asarray(result.time_overlap)
        time_overlap = np.zeros((nstates, nstates), dtype=root_overlap.dtype)
        for i in range(nstates):
            for j in range(nstates):
                if ms_of_state[i] == ms_of_state[j]:
                    time_overlap[i, j] = root_overlap[root_of_state[i], root_of_state[j]]

    return {
        "multiplicity": multiplicity,
        "labels": labels,
        "root_of_state": root_of_state,
        "energies": energies[root_of_state],
        "gradients": result.gradients,
        "nac_vectors": result.nac_vectors,
        "time_overlap": time_overlap,
    }


def _pyscf_compute_adi_spin_manifolds(q, params, full_id):
    """Build a block-diagonal Libra Hamiltonian from PySCF spin manifolds.

    Each entry of ``params['spin_manifolds']`` represents an independent
    spin-free PySCF calculation and defines ``spin`` (the multiplicity),
    ``nroots``, and either ``es_strategy`` or ``strategy_factory``. Every root
    is expanded over its ``2*S+1`` spin projections by default. Without SOC,
    matrix elements between different multiplicities and different Ms values
    are zero.
    """
    itraj = _get_trajectory_index(full_id)
    geometry = _q_to_geometry(q, itraj, params["atom_labels"])
    manifolds = params["spin_manifolds"]
    if not isinstance(manifolds, (list, tuple)) or not manifolds:
        raise ValueError("spin_manifolds must be a non-empty list of dictionaries")
    if any(not isinstance(item, dict) or "spin" not in item for item in manifolds):
        raise ValueError("Each spin_manifolds entry must be a dictionary with 'spin'")
    multiplicities = [int(item["spin"]) for item in manifolds]
    if len(set(multiplicities)) != len(multiplicities):
        raise ValueError("spin_manifolds must not repeat a spin multiplicity")

    results = [
        _compute_spin_manifold(geometry, params, itraj, manifold, index)
        for index, manifold in enumerate(manifolds)
    ]
    natoms = len(params["atom_labels"])
    ndof = 3 * natoms
    nstates = sum(len(result["labels"]) for result in results)
    dt = float(params.get("dt", 41.0))
    energy_zero = float(params.get("energy_zero", 0.0))

    obj = SimpleNamespace()
    obj.spin_labels = [
        (result["multiplicity"], root, twice_ms)
        for result in results for root, twice_ms in result["labels"]
    ]
    obj.ms_labels = [label[1:] for label in obj.spin_labels]
    obj.ham_adi = CMATRIX(nstates, nstates)
    obj.nac_adi = CMATRIX(nstates, nstates)
    obj.hvib_adi = CMATRIX(nstates, nstates)
    obj.time_overlap_adi = CMATRIX(nstates, nstates)
    obj.overlap_adi = CMATRIX(nstates, nstates)
    obj.basis_transform = CMATRIX(nstates, nstates)
    obj.d1ham_adi = CMATRIXList()
    obj.dc1_adi = CMATRIXList()
    for _ in range(ndof):
        obj.d1ham_adi.append(CMATRIX(nstates, nstates))
        obj.dc1_adi.append(CMATRIX(nstates, nstates))

    offset = 0
    for result in results:
        local_nstates = len(result["labels"])
        for i in range(local_nstates):
            gi = offset + i
            energy = complex(result["energies"][i] - energy_zero)
            obj.ham_adi.set(gi, gi, energy)
            obj.hvib_adi.set(gi, gi, energy)
            obj.basis_transform.set(gi, gi, 1.0 + 0.0j)
            obj.overlap_adi.set(gi, gi, 1.0 + 0.0j)
            for j in range(local_nstates):
                obj.time_overlap_adi.set(
                    gi, offset + j, complex(result["time_overlap"][i, j])
                )

        if result["gradients"] is not None:
            for i, root in enumerate(result["root_of_state"]):
                gradient = result["gradients"][root]
                if gradient is None:
                    continue
                for atom in range(natoms):
                    for xyz in range(3):
                        obj.d1ham_adi[3 * atom + xyz].set(
                            offset + i, offset + i, complex(gradient[atom, xyz])
                        )

        if result["nac_vectors"] is not None:
            nac_vectors = np.asarray(result["nac_vectors"])
            for i, (root_i, ms_i) in enumerate(result["labels"]):
                for j, (root_j, ms_j) in enumerate(result["labels"]):
                    if ms_i != ms_j:
                        continue
                    for atom in range(natoms):
                        for xyz in range(3):
                            obj.dc1_adi[3 * atom + xyz].set(
                                offset + i,
                                offset + j,
                                complex(nac_vectors[root_i, root_j, atom, xyz]),
                            )
        offset += local_nstates

    for i in range(nstates):
        for j in range(i + 1, nstates):
            dij = (
                obj.time_overlap_adi.get(i, j)
                - obj.time_overlap_adi.get(j, i)
            ) / (2.0 * dt)
            obj.nac_adi.set(i, j, dij)
            obj.nac_adi.set(j, i, -dij.conjugate())
            obj.hvib_adi.set(i, j, -1.0j * dij)
            obj.hvib_adi.set(j, i, 1.0j * dij.conjugate())
    return obj


# =============================================================================
# Main Libra callback
# =============================================================================

def pyscf_compute_adi(q,params,full_id):

    if params.get("spin_manifolds") is not None:
        return _pyscf_compute_adi_spin_manifolds(q, params, full_id)

    # 1. Libra input -> generic ES input
    itraj = _get_trajectory_index(full_id)
    geometry = _q_to_geometry(q, itraj, params["atom_labels"] )
    request = _params_to_request(params, itraj)

    # 2. Get this trajectory's persistent ES strategy.
    #
    # ES_Strategy tracks its own previous-geometry snapshot internally
    # (get_previous_state/snapshot_state), so the *same* instance must be reused
    # across calls for a given trajectory or time-overlaps and NACs are never
    # available. Instances are therefore memoized per trajectory index -- the
    # same {itraj: ...} convention the mopac, dftbplus, and cp2k interfaces use
    # for their own per-trajectory caches.

    strategies = params.setdefault("es_strategies", {})

    if itraj not in strategies:
        strategy_spec = params.get(  "strategy_factory",params.get("es_strategy") )

        if strategy_spec is None:
            raise KeyError("Missing strategy specification: expected 'strategy_factory' or 'es_strategy'.")

        if callable(strategy_spec):
            strategies[itraj] = strategy_spec()
        elif hasattr(strategy_spec, "copy"):
            strategies[itraj] = strategy_spec.copy()
        elif hasattr(strategy_spec, "clone"):
            strategies[itraj] = strategy_spec.clone()
        else:
            strategies[itraj] = strategy_spec

    current = strategies[itraj]

    # 3. Run ES calculation.
    result = strategies[itraj].compute_result(geometry, request)

    # 4. Generic ES result -> Libra result
    obj = _es_result_to_libra(
        result,
        request,
        natoms=len(params["atom_labels"]),
        dt=float(params.get("dt", 41.0)),
    )

    return obj
