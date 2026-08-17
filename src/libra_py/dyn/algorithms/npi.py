"""Meek-Levine norm-preserving interpolation for NumPy and PyTorch arrays."""

from __future__ import annotations

import math
from typing import Any

import numpy as np


ORTHOGONALITY_TOL = 1.0e-6
SINGULARITY_TOL = 1.0e-12


def _array_module(value: Any):
    if value.__class__.__module__.startswith("torch"):
        import torch

        return torch
    return np


def _scalar(value: Any) -> float:
    if hasattr(value, "detach"):
        value = value.detach().cpu()
    return float(value.item() if hasattr(value, "item") else value)


def _is_complex(xp: Any, value: Any) -> bool:
    if xp is np:
        return np.iscomplexobj(value)
    return bool(xp.is_complex(value))


def _asarray(time_overlap: Any, backend: Any | None):
    if backend is None:
        xp = _array_module(time_overlap)
    elif backend == "numpy" or backend is np:
        xp = np
    elif backend == "torch" or getattr(backend, "__name__", None) == "torch":
        import torch

        xp = torch
    else:
        # TensorStorage exposes either NumPy/PyTorch itself or a facade with an
        # asarray operation. Infer the native array module after conversion.
        converted = backend.asarray(time_overlap)
        return _array_module(converted), converted

    if xp is np:
        return xp, np.asarray(time_overlap)
    return xp, xp.as_tensor(time_overlap)


def _validate(overlap: Any, dt: float, xp: Any, tol: float) -> Any:
    if overlap.ndim < 2 or overlap.shape[-1] != overlap.shape[-2] or overlap.shape[-1] == 0:
        raise ValueError(
            "nac_npi: time_overlap must have shape (..., nstates, nstates) with "
            f"nstates > 0; received {tuple(overlap.shape)}. Supply overlaps between "
            "the same complete state set at both time steps."
        )
    if not math.isfinite(dt) or dt <= 0.0:
        raise ValueError(
            f"nac_npi: dt must be finite and positive; received {dt}. Pass the "
            "positive time interval separating the two electronic-structure steps."
        )
    if not bool(xp.all(xp.isfinite(overlap))):
        raise ValueError(
            "nac_npi: time_overlap contains a non-finite element. Check the "
            "electronic-structure calculation and state-overlap construction."
        )

    if _is_complex(xp, overlap):
        max_imaginary = _scalar(xp.max(xp.abs(overlap.imag)))
        if max_imaginary > tol:
            raise ValueError(
                "nac_npi: this Meek-Levine implementation requires real overlaps, "
                f"but max|Im(S)| = {max_imaginary} exceeds {tol}. Fix the electronic-state "
                "gauge to make the overlaps real or use a complex-unitary TDC method."
            )
        real_overlap = overlap.real
    else:
        real_overlap = overlap

    diagonal = xp.diagonal(real_overlap, dim1=-2, dim2=-1) if xp is not np else np.diagonal(
        real_overlap, axis1=-2, axis2=-1
    )
    min_diagonal = _scalar(xp.min(diagonal))
    if min_diagonal < -tol:
        raise ValueError(
            "nac_npi: a diagonal overlap is negative (minimum = "
            f"{min_diagonal}), indicating discontinuous electronic-state phases. "
            "Phase-match and reorder the states so corresponding-state overlaps are non-negative."
        )

    transpose = real_overlap.swapaxes(-1, -2) if xp is np else real_overlap.transpose(-1, -2)
    gram = transpose @ real_overlap
    identity = xp.eye(real_overlap.shape[-1], dtype=real_overlap.dtype)
    if xp is not np:
        identity = identity.to(device=real_overlap.device)
    orthogonality_error = _scalar(xp.max(xp.abs(gram - identity)))
    if orthogonality_error > tol:
        raise ValueError(
            "nac_npi: time_overlap is not orthogonal; max|S^T S - I| = "
            f"{orthogonality_error} exceeds {tol}. Use the same complete state space "
            "at both steps and apply polar/Lowdin orthogonalization after state matching."
        )

    determinants = xp.linalg.det(real_overlap)
    min_determinant = _scalar(xp.min(determinants))
    if min_determinant <= 0.0:
        raise ValueError(
            "nac_npi: every phase-matched overlap must be a proper rotation, but "
            f"the minimum det(S) is {min_determinant}. Correct state phases/permutations "
            "so every determinant is positive before applying NPI."
        )
    return real_overlap


def nac_npi(
    time_overlap: Any,
    dt: float,
    *,
    backend: Any | None = None,
    orthogonality_tol: float = ORTHOGONALITY_TOL,
) -> Any:
    """Return interval-averaged real TDC matrices using Meek-Levine NPI.

    Parameters
    ----------
    time_overlap
        A NumPy array or PyTorch tensor with shape ``(..., nstates, nstates)``.
        This includes the ``(ntraj, ntbf, nstates, nstates)`` layout used by
        :class:`libra_py.dyn.core.storage.TensorStorage`. Complex storage is
        accepted when its imaginary component is negligible.
    dt
        Positive time interval between the two sets of electronic states.
    backend
        Optional ``numpy``, ``torch``, backend name, or dyn backend facade.
        By default the backend is inferred without moving or copying tensors.
    orthogonality_tol
        Maximum accepted element of ``|S.T @ S - I|`` and imaginary overlap.

    Returns
    -------
    numpy.ndarray or torch.Tensor
        Antisymmetric TDC matrices with the same shape, backend, dtype family,
        and (for PyTorch) device as ``time_overlap``.

    Notes
    -----
    This is the real, pairwise Meek-Levine prescription. Inputs must already
    be state-matched, phase-corrected, and orthogonalized proper rotations.
    """
    xp, source = _asarray(time_overlap, backend)
    dt = _scalar(dt)
    tol = float(orthogonality_tol)
    if not math.isfinite(tol) or tol <= 0.0:
        raise ValueError("nac_npi: orthogonality_tol must be finite and positive")
    overlap = _validate(source, dt, xp, tol)

    result = xp.zeros_like(overlap)
    zero = xp.zeros_like(overlap[..., 0, 0])
    one = xp.ones_like(zero)
    nstates = overlap.shape[-1]

    for i in range(nstates):
        for j in range(i + 1, nstates):
            w00 = xp.clip(overlap[..., i, i], -1.0, 1.0)
            w01 = xp.clip(overlap[..., i, j], -1.0, 1.0)
            w10 = xp.clip(overlap[..., j, i], -1.0, 1.0)
            w11 = xp.clip(overlap[..., j, j], -1.0, 1.0)

            a0 = xp.arccos(w00) - xp.arcsin(w01)
            b0 = xp.arccos(w00) + xp.arcsin(w01)
            c0 = xp.arccos(w11) - xp.arcsin(w10)
            d0 = xp.arccos(w11) + xp.arcsin(w10)
            a = -xp.sinc(a0 / math.pi)
            b = xp.sinc(b0 / math.pi)
            c = xp.sinc(c0 / math.pi)
            d = xp.sinc(d0 / math.pi)

            wlj = xp.sqrt(xp.clip(1.0 - w00 * w00 - w10 * w10, 0.0, None))
            nonzero_wlj = wlj != 0.0
            safe_wlj = xp.where(nonzero_wlj, wlj, one)
            wlk = xp.clip(-(w01 * w00 + w11 * w10) / safe_wlj, -1.0, 1.0)
            swlj = xp.arcsin(wlj)
            swlk = xp.arcsin(wlk)
            radicand = xp.clip((1.0 - wlj * wlj) * (1.0 - wlk * wlk), 0.0, None)
            root = xp.sqrt(radicand)
            denominator = swlj * swlj - swlk * swlk
            singular = xp.abs(denominator) <= SINGULARITY_TOL
            safe_denominator = xp.where(singular, one, denominator)
            regular_e = (
                2.0 * swlj
                * (wlj * wlk * swlj + (root - 1.0) * swlk)
                / safe_denominator
            )
            limit_e = xp.where(wlk < 0.0, -wlj * wlj, wlj * wlj)
            e = xp.where(nonzero_wlj, xp.where(singular, limit_e, regular_e), zero)

            tdc = 0.5 / dt * (
                xp.arccos(w00) * (a + b) + xp.arcsin(w10) * (c + d) + e
            )
            result[..., j, i] = tdc
            result[..., i, j] = -tdc

    if _is_complex(xp, source):
        if xp is np:
            return result.astype(source.dtype, copy=False)
        return result.to(dtype=source.dtype)
    return result
