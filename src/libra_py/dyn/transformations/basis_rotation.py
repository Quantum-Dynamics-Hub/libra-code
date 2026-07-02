from __future__ import annotations

from typing import Any


def amplitudes_dia_to_adi(
    ampl_dia: Any,
    basis_transform: Any,
    ovlp_dia: Any,
    backend: Any,
) -> Any:
    """
    Transform diabatic amplitudes to the adiabatic representation.

    If |psi_adi> = |psi_dia> U and U.H S U = I, then C_adi = U.H S C_dia.
    Leading batch dimensions are preserved.
    """

    u_h = backend.conjugate_transpose(basis_transform)
    return backend.einsum("...ij,...jk,...k->...i", u_h, ovlp_dia, ampl_dia)


def amplitudes_adi_to_dia(
    ampl_adi: Any,
    basis_transform: Any,
    backend: Any,
) -> Any:
    """
    Transform adiabatic amplitudes to the diabatic representation.

    If |psi_adi> = |psi_dia> U, then C_dia = U C_adi.
    """

    return backend.einsum("...ij,...j->...i", basis_transform, ampl_adi)


def matrix_dia_to_adi(matrix_dia: Any, basis_transform: Any, backend: Any) -> Any:
    """Transform an operator matrix from diabatic to adiabatic coordinates."""

    u_h = backend.conjugate_transpose(basis_transform)
    return backend.einsum("...ij,...jk,...kl->...il", u_h, matrix_dia, basis_transform)


def matrix_adi_to_dia(matrix_adi: Any, basis_transform: Any, backend: Any) -> Any:
    """Transform an operator matrix from adiabatic to diabatic coordinates."""

    u_h = backend.conjugate_transpose(basis_transform)
    return backend.einsum("...ij,...jk,...kl->...il", basis_transform, matrix_adi, u_h)


def rotate_operator(matrix: Any, transform: Any, backend: Any) -> Any:
    """Return T.H A T for a batched matrix/operator transform."""

    t_h = backend.conjugate_transpose(transform)
    return backend.einsum("...ij,...jk,...kl->...il", t_h, matrix, transform)


def storage_amplitudes_dia_to_adi(storage: Any, traj: Any) -> Any:
    """
    Transform active TensorStorage amplitudes from diabatic to adiabatic.

    This storage wrapper intentionally lives in transformations, not
    hamiltonians, because Hamiltonian construction should not generally depend
    on electronic amplitudes.
    """

    idx = traj.tbf_ids
    storage.ampl_adi[traj.id, idx] = amplitudes_dia_to_adi(
        storage.ampl_dia[traj.id, idx],
        storage.basis_transform[traj.id, idx],
        storage.ovlp_dia[traj.id, idx],
        storage.backend,
    )
    return storage.ampl_adi[traj.id, idx]


def storage_amplitudes_adi_to_dia(storage: Any, traj: Any) -> Any:
    """Transform active TensorStorage amplitudes from adiabatic to diabatic."""

    idx = traj.tbf_ids
    storage.ampl_dia[traj.id, idx] = amplitudes_adi_to_dia(
        storage.ampl_adi[traj.id, idx],
        storage.basis_transform[traj.id, idx],
        storage.backend,
    )
    return storage.ampl_dia[traj.id, idx]
