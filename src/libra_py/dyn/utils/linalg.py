from __future__ import annotations

from typing import Any

import numpy as np

from ..backends.arrays import xp_of, adjoint


def generalized_eigh(hamiltonian: Any, overlap: Any) -> tuple[Any, Any]:
    """
    Solve H U = S U E for one Hermitian matrix or a batch of matrices.

    Parameters
    ----------
    hamiltonian
        Matrix or batched matrices with trailing shape `(nstates, nstates)`.
    overlap
        Positive-definite overlap matrix or batched matrices with the same
        trailing shape.

    Returns
    -------
    eigenvalues, eigenvectors
        Eigenvalues are sorted ascending. Eigenvectors satisfy `U.H S U = I`
        up to numerical precision.

    Notes
    -----
    The current dyn backend facade does not expose a generalized Hermitian
    eigensolver, so this implementation uses NumPy and a Cholesky reduction:
    `S = L L.H`, `H_orth = L^{-1} H L^{-H}`.
    """

    h = np.asarray(hamiltonian)
    s = np.asarray(overlap)
    if h.shape != s.shape:
        raise ValueError(
            f"hamiltonian and overlap shapes must match; got {h.shape} and {s.shape}"
        )
    if h.shape[-1] != h.shape[-2]:
        raise ValueError("generalized_eigh expects square trailing matrices")

    batch_shape = h.shape[:-2]
    ns = h.shape[-1]
    h_flat = h.reshape((-1, ns, ns))
    s_flat = s.reshape((-1, ns, ns))
    e_flat = np.empty((h_flat.shape[0], ns), dtype=float)
    u_flat = np.empty_like(h_flat, dtype=complex)

    for ibatch, (h_mat, s_mat) in enumerate(zip(h_flat, s_flat)):
        if ns == 1:
            e_flat[ibatch, 0] = (h_mat[0, 0] / s_mat[0, 0]).real
            u_flat[ibatch, 0, 0] = 1.0 / np.sqrt(s_mat[0, 0])
            continue

        chol = np.linalg.cholesky(s_mat)
        inv_chol = np.linalg.inv(chol)
        ortho_ham = inv_chol @ h_mat @ inv_chol.conj().T
        values, vectors = np.linalg.eigh(ortho_ham)
        u_flat[ibatch] = inv_chol.conj().T @ vectors
        e_flat[ibatch] = values

    return e_flat.reshape((*batch_shape, ns)), u_flat.reshape((*batch_shape, ns, ns))


def diag_embed(values: Any, dtype=complex):
    """Return batched diagonal matrices with `values` on the last-axis diagonal."""

    values = np.asarray(values)
    ns = values.shape[-1]
    out = np.zeros((*values.shape, ns), dtype=dtype)
    idx = np.arange(ns)
    out[..., idx, idx] = values
    return out


# Spectral metric operations for coherent nonorthogonal expansions.

def _orthogonalizer(S, rcond):
    xp = xp_of(S)
    if not 0 <= rcond < 1:
        raise ValueError("overlap_rcond must lie in [0,1)")
    vals, vecs = xp.linalg.eigh((S+adjoint(S))/2)
    keep = vals > rcond*vals.max()
    if not bool(keep.any()):
        raise ValueError("overlap has no positive retained eigenvalues")
    return vecs[:, keep]/xp.sqrt(vals[keep])[None, :]


def metric_solve(S, rhs, rcond=1e-12):
    """Solve S*x=rhs in its positive spectral subspace (relative cutoff)."""
    W = _orthogonalizer(S, rcond)
    return W @ (adjoint(W) @ rhs)
