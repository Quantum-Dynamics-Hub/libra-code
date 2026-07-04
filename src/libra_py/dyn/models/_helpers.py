from __future__ import annotations


def merged_params(params, defaults):
    result = dict(defaults)
    result.update(params or {})
    return result


def complex_dtype(x):
    if x.__class__.__module__.startswith("torch"):
        import torch

        return torch.complex128 if x.dtype == torch.float64 else torch.complex64
    return complex


def zeros_from(x, shape):
    if x.__class__.__module__.startswith("torch"):
        import torch

        return torch.zeros(shape, dtype=complex_dtype(x), device=x.device)
    import numpy as np

    return np.zeros(shape, dtype=complex)


def real_zeros_like(x):
    if x.__class__.__module__.startswith("torch"):
        import torch

        return torch.zeros_like(x)
    import numpy as np

    return np.zeros_like(x)


def fill_symmetric(matrix, h00, h11, h01):
    matrix[..., 0, 0] = h00
    matrix[..., 1, 1] = h11
    matrix[..., 0, 1] = h01
    matrix[..., 1, 0] = h01


def two_state_with_derivatives(model, x, h00, h11, h01, dh00, dh11, dh01, dof=0):
    shape = tuple(x.shape)
    H = zeros_from(x, (*shape, 2, 2))
    dH = zeros_from(x, (*shape, model.ndof, 2, 2))
    fill_symmetric(H, h00, h11, h01)
    fill_symmetric(dH[..., dof, :, :], dh00, dh11, dh01)
    return H, dH


def one_state_with_derivatives(model, value, derivatives):
    shape = tuple(value.shape)
    H = zeros_from(value, (*shape, 1, 1))
    dH = zeros_from(value, (*shape, model.ndof, 1, 1))
    H[..., 0, 0] = value
    for dof, derivative in enumerate(derivatives):
        dH[..., dof, 0, 0] = derivative
    return H, dH
