from __future__ import annotations

from typing import Any, Callable

import numpy as np

from .aux import as_result_mapping, check_matrix_shape, write_active_slice
from .basic import init_hamiltonian_storage
from ..utils.linalg import diag_embed, generalized_eigh


ADIABATIC_RESULT_FIELDS = {
    "ham_adi": "ham_adi",
    "H_adi": "ham_adi",
    "nac_adi": "nac_adi",
    "NAC_adi": "nac_adi",
    "hvib_adi": "hvib_adi",
    "Hvib_adi": "hvib_adi",
    "basis_transform": "basis_transform",
    "time_overlap_adi": "time_overlap_adi",
    "dc1_adi": "dc1_adi",
    "DC1_adi": "dc1_adi",
    "d1ham_adi": "d1ham_adi",
    "dH_adi": "d1ham_adi",
    "d2ham_adi": "d2ham_adi",
    "d2H_adi": "d2ham_adi",
}


def compute_adiabatic(
    storage: Any,
    traj: Any,
    model_fn: Callable | None = None,
    params: Any = None,
    der_lvl: int = 1,
    from_diabatic: bool = False,
) -> Any:
    """
    Compute or ingest adiabatic Hamiltonian data.

    If `from_diabatic` is True, this translates
    nHamiltonian::compute_adiabatic(der_lvl): the active diabatic Hamiltonian
    and overlap are diagonalized and derivative couplings are transformed when
    `der_lvl >= 1`. Otherwise `model_fn` is called and any returned adiabatic
    fields are validated and written into TensorStorage.
    """

    init_hamiltonian_storage(storage, der_lvl=der_lvl)
    if from_diabatic:
        compute_adiabatic_from_diabatic(storage, traj, der_lvl=der_lvl)
        return storage
    if model_fn is None:
        raise ValueError("model_fn is required unless from_diabatic=True")

    idx = traj.tbf_ids
    result = as_result_mapping(
        model_fn(
            R=storage.q[traj.id, idx],
            P=storage.p[traj.id, idx],
            storage=storage,
            traj=traj,
            params=params,
        )
    )

    for source, target in ADIABATIC_RESULT_FIELDS.items():
        value = result.get(source)
        if value is None:
            continue
        _validate_adiabatic_field(storage, value, source, target)
        write_active_slice(storage, traj, target, value)

    if "gs_kinetic_energy" in result:
        storage.gs_kinetic_energy = float(result["gs_kinetic_energy"])

    if "hvib_adi" not in result and "Hvib_adi" not in result:
        build_adiabatic_hvib(storage, traj)

    return storage


def compute_adiabatic_from_diabatic(
    storage: Any,
    traj: Any,
    der_lvl: int = 1,
    energy_gap_floor: float = 1.0e-25,
) -> Any:
    """
    Diagonalize stored diabatic data and populate adiabatic quantities.

    The derivative-coupling translation follows the C++ implementation:

    tmp = U.H dH_dia[n] U
    dtilda = U.H dc1_dia[n] U H_adi
    tmp -= dtilda + dtilda.H
    d1ham_adi[n] = diag(tmp)
    dc1_adi[n, i, j] = tmp[i, j] / (E_j - E_i)
    """

    idx = traj.tbf_ids
    ham_dia = storage.ham_dia[traj.id, idx]
    ovlp_dia = storage.ovlp_dia[traj.id, idx]

    energies, transform = generalized_eigh(ham_dia, ovlp_dia)
    storage.ham_adi[traj.id, idx] = diag_embed(energies)
    storage.basis_transform[traj.id, idx] = transform

    if der_lvl >= 1:
        _compute_adiabatic_derivatives(
            storage,
            traj,
            transform,
            energies,
            energy_gap_floor,
        )

    build_adiabatic_hvib(storage, traj)
    return storage


def build_adiabatic_hvib(storage: Any, traj: Any) -> Any:
    """Build `hvib_adi = ham_adi - i * nac_adi` for active TBF slices."""

    idx = traj.tbf_ids
    if storage.ham_adi is None or storage.hvib_adi is None:
        return storage
    value = storage.ham_adi[traj.id, idx]
    if storage.nac_adi is not None:
        value = value - 1j * storage.nac_adi[traj.id, idx]
    storage.hvib_adi[traj.id, idx] = value
    return storage


def _compute_adiabatic_derivatives(
    storage: Any,
    traj: Any,
    transform: Any,
    energies: Any,
    energy_gap_floor: float,
) -> None:
    idx = traj.tbf_ids
    u = np.asarray(transform)
    u_h = np.swapaxes(np.conjugate(u), -1, -2)
    ham_adi = diag_embed(energies)
    d1_dia = np.asarray(storage.d1ham_dia[traj.id, idx])
    dc1_dia = np.asarray(storage.dc1_dia[traj.id, idx])

    tmp = np.einsum("...ij,...djk,...kl->...dil", u_h, d1_dia, u)
    dtilda = np.einsum("...ij,...djk,...kl,...lm->...dim", u_h, dc1_dia, u, ham_adi)
    tmp = tmp - (dtilda + np.swapaxes(np.conjugate(dtilda), -1, -2))

    d1_adi = np.zeros_like(tmp)
    diag_idx = np.arange(storage.nadi)
    d1_adi[..., diag_idx, diag_idx] = tmp[..., diag_idx, diag_idx]

    dc1_adi = np.zeros_like(tmp)
    for i in range(storage.nadi):
        for j in range(i + 1, storage.nadi):
            gap = energies[..., j] - energies[..., i]
            gap = np.where(np.abs(gap) < energy_gap_floor, energy_gap_floor, gap)
            val = tmp[..., i, j] / gap[..., None]
            dc1_adi[..., i, j] = val
            dc1_adi[..., j, i] = -val

    storage.d1ham_adi[traj.id, idx] = d1_adi
    storage.dc1_adi[traj.id, idx] = dc1_adi


def _validate_adiabatic_field(storage: Any, value: Any, source: str, target: str) -> None:
    ns = storage.nadi
    nd = storage.nnucl
    if target in ("dc1_adi", "d1ham_adi"):
        check_matrix_shape(value, (nd, ns, ns), source)
    elif target == "d2ham_adi":
        check_matrix_shape(value, (nd, nd, ns, ns), source)
    else:
        check_matrix_shape(value, (ns, ns), source)
