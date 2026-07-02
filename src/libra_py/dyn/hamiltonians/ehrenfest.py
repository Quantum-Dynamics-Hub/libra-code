from __future__ import annotations

from typing import Any

import numpy as np


def ehrenfest_energy_dia(storage: Any, traj: Any, amplitudes: Any = None) -> Any:
    """
    Compute diabatic Ehrenfest energy for active TBFs.

    E = C.H H_dia C / (C.H S_dia C)
    """

    idx = traj.tbf_ids
    c = _as_array(amplitudes, storage.ampl_dia[traj.id, idx])
    h = np.asarray(storage.ham_dia[traj.id, idx])
    s = np.asarray(storage.ovlp_dia[traj.id, idx])
    norm = _expectation(c, s)
    return _expectation(c, h) / norm


def ehrenfest_energy_adi(
    storage: Any,
    traj: Any,
    amplitudes: Any = None,
    transform: Any = None,
) -> Any:
    """
    Compute adiabatic Ehrenfest energy for active TBFs.

    E = C.H T.H H_adi T C / (C.H C). If no transform is supplied, T is the
    identity, matching the default C++ overload.
    """

    idx = traj.tbf_ids
    c = _as_array(amplitudes, storage.ampl_adi[traj.id, idx])
    h = np.asarray(storage.ham_adi[traj.id, idx])
    if transform is not None:
        t = np.asarray(transform)
        h = np.einsum("...ij,...jk,...kl->...il", _h(t), h, t)
    norm = np.einsum("...i,...i->...", np.conjugate(c), c)
    return _expectation(c, h) / norm


def ehrenfest_force_tensors_adi(
    storage: Any,
    traj: Any,
    amplitudes: Any = None,
    option: int = 0,
    transform: Any = None,
) -> Any:
    """
    Return adiabatic mean-field force tensors for active TBFs.

    Each tensor F[n] satisfies force[n] = C.H F[n] C. With `option=0`, NAC
    correction terms are included. With `option=1`, only -dH/dR is used.
    """

    idx = traj.tbf_ids
    c = _as_array(amplitudes, storage.ampl_adi[traj.id, idx])
    h = np.asarray(storage.ham_adi[traj.id, idx])
    d1 = np.asarray(storage.d1ham_adi[traj.id, idx])
    dc = np.asarray(storage.dc1_adi[traj.id, idx])
    norm = np.einsum("...i,...i->...", np.conjugate(c), c)

    if transform is None:
        t = _identity_like_hamiltonian(h)
    else:
        t = np.asarray(transform)

    t_h = _h(t)
    h_t = np.einsum("...ij,...jk,...kl->...il", t_h, h, t)
    d1_t = np.einsum("...ij,...djk,...kl->...dil", t_h, d1, t)
    if option == 0:
        dc_t = np.einsum("...ij,...djk,...kl->...dil", t_h, dc, t)
        tmp = np.einsum("...dij,...jk->...dik", _h(dc_t), h_t)
        tmp = tmp + _h(tmp)
        tensors = -(d1_t - tmp) / norm[..., None, None, None]
    elif option == 1:
        tensors = -d1_t / norm[..., None, None, None]
    else:
        raise ValueError("option must be 0 or 1")

    return tensors


def ehrenfest_force_tensors_dia(
    storage: Any,
    traj: Any,
    amplitudes: Any = None,
    option: int = 0,
) -> Any:
    """
    Return diabatic mean-field force tensors for active TBFs.

    With `option=0`, derivative-coupling corrections are included. With
    `option=1`, only -dH/dR is used.
    """

    idx = traj.tbf_ids
    c = _as_array(amplitudes, storage.ampl_dia[traj.id, idx])
    h = np.asarray(storage.ham_dia[traj.id, idx])
    s = np.asarray(storage.ovlp_dia[traj.id, idx])
    d1 = np.asarray(storage.d1ham_dia[traj.id, idx])
    dc = np.asarray(storage.dc1_dia[traj.id, idx])
    inv_s = np.linalg.inv(s)
    norm = _expectation(c, s)

    if option == 0:
        tmp = np.einsum("...dij,...jk,...kl->...dil", _h(dc), inv_s, h)
        tmp = tmp + _h(tmp)
        tensors = -(d1 - tmp) / norm[..., None, None, None]
    elif option == 1:
        tensors = -d1 / norm[..., None, None, None]
    else:
        raise ValueError("option must be 0 or 1")

    return tensors


def ehrenfest_forces_adi(
    storage: Any,
    traj: Any,
    amplitudes: Any = None,
    option: int = 0,
    transform: Any = None,
    gamma: float | Any = 0.0,
) -> Any:
    """
    Compute adiabatic Ehrenfest forces for active TBFs.

    The return shape is `(ntbf, ndof)` for a trajectory slice. The optional
    Meyer-Miller/SQC `gamma` term follows the C++ implementation.
    """

    idx = traj.tbf_ids
    c = _as_array(amplitudes, storage.ampl_adi[traj.id, idx])
    tensors = ehrenfest_force_tensors_adi(
        storage,
        traj,
        c,
        option=option,
        transform=transform,
    )
    forces = _expectation_by_dof(c, tensors)

    gamma_arr = np.asarray(gamma)
    if np.any(gamma_arr != 0.0):
        d1 = np.asarray(storage.d1ham_adi[traj.id, idx])
        trace_grad = np.trace(d1, axis1=-2, axis2=-1)
        norm = np.einsum("...i,...i->...", np.conjugate(c), c)
        forces += np.asarray(gamma)[..., None] * trace_grad / norm[..., None]

    return forces


def ehrenfest_forces_dia(
    storage: Any,
    traj: Any,
    amplitudes: Any = None,
    option: int = 0,
) -> Any:
    """Compute diabatic Ehrenfest forces for active TBFs."""

    idx = traj.tbf_ids
    c = _as_array(amplitudes, storage.ampl_dia[traj.id, idx])
    tensors = ehrenfest_force_tensors_dia(storage, traj, c, option=option)
    return _expectation_by_dof(c, tensors)


def _expectation(c: Any, matrix: Any) -> Any:
    return np.einsum("...i,...ij,...j->...", np.conjugate(c), matrix, c)


def _expectation_by_dof(c: Any, tensors: Any) -> Any:
    return np.einsum("...i,...dij,...j->...d", np.conjugate(c), tensors, c)


def _h(matrix: Any) -> Any:
    return np.swapaxes(np.conjugate(matrix), -1, -2)


def _as_array(value: Any, default: Any) -> Any:
    return np.asarray(default if value is None else value)


def _identity_like_hamiltonian(hamiltonian: Any) -> Any:
    ns = hamiltonian.shape[-1]
    eye = np.eye(ns, dtype=complex)
    return np.broadcast_to(eye, hamiltonian.shape)
