# Copyright (C) 2022-2026 Matthew Dutra, Alexey V. Akimov
# SPDX-License-Identifier: GPL-3.0-or-later

"""Array utilities shared by Gaussian dynamics; preserve dtype/device/autograd."""
import sys
import numpy as np

def xp_of(x):
    torch = sys.modules.get("torch")
    return torch if torch is not None and isinstance(x, torch.Tensor) else np


def asarray(x, like, complex=False):
    xp = xp_of(like)
    if xp is np:
        dtype = np.result_type(like.dtype, np.complex64) if complex else like.dtype
        return np.asarray(x, dtype=dtype)
    dtype = (xp.complex128 if like.dtype in (xp.float64, xp.complex128) else xp.complex64) if complex else like.dtype
    return xp.as_tensor(x, dtype=dtype, device=like.device)


def zeros(shape, like, complex=False):
    xp = xp_of(like)
    dtype = asarray(0., like, complex).dtype
    if xp is np:
        return np.zeros(shape, dtype=dtype)
    return xp.zeros(shape, dtype=dtype, device=like.device)


def stack(values, axis=0):
    xp = xp_of(values[0])
    return xp.stack(values, axis=axis) if xp is np else xp.stack(values, dim=axis)


def copy(x):
    return x.copy() if xp_of(x) is np else x.clone()


def adjoint(x):
    return x.conj().swapaxes(-1, -2)


def indices(values, like):
    xp = xp_of(like)
    if xp is np:
        return np.asarray(values, dtype=int)
    return xp.as_tensor(values, dtype=xp.long, device=like.device)


def to_numpy(x):
    """Detach only at explicit file/plot boundaries."""
    return np.asarray(x) if xp_of(x) is np else x.detach().cpu().resolve_conj().resolve_neg().numpy()


def scalar(x):
    """Detach a scalar used ONLY for nondifferentiable control-flow decisions."""
    return float(x) if xp_of(x) is np else x.detach().item()
