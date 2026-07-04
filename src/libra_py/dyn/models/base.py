from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np


@dataclass
class AnalyticalHamiltonianModel:
    """
    Base class for backend-aware analytical Hamiltonian models.

    Conventions
    -----------
    Subclasses implement ``diabatic_with_derivatives(Q, params)`` with the
    exact-dynamics coordinate convention ``Q.shape == (ndof, *grid_shape)``.
    The returned tensors use trailing electronic dimensions:

    ``H_dia.shape == (*grid_shape, nstates, nstates)``
        Diabatic potential/Hamiltonian matrix.

    ``dH_dia.shape == (*grid_shape, ndof, nstates, nstates)``
        First spatial derivatives ``dH_dia[..., a, i, j] = dH_ij/dQ_a``.

    The same object is callable by ``HamiltonianEngine`` for trajectory slices.
    In that mode ``R.shape == (ntbf, ndof)`` is converted internally to the
    exact-style ``Q.shape == (ndof, ntbf)``.

    Analytical models define only their native diabatic representation.
    Adiabatic energies, derivative couplings, velocity-projected NACs,
    vibronic Hamiltonians, gauge fixes, and basis rotations belong to the
    Hamiltonian and transformation layers.

    Model result dictionary
    -----------------------
    ``evaluate`` returns ``H_dia``, ``dH_dia``, zero ``DC1_dia``, and identity
    ``S_dia``. This mirrors the legacy ``libra_py.models`` objects
    ``ham_dia``, ``d1ham_dia``, ``dc1_dia``, and ``ovlp_dia`` but uses
    backend-friendly arrays rather than Libra C++ matrix wrappers.

    Units
    -----
    Parameters and coordinates are in atomic units unless a model docstring
    explicitly states otherwise. Constants translated from older model files
    are converted during parameter construction.
    """

    params: dict[str, Any] = field(default_factory=dict)
    backend: str = "numpy"

    ndof: int = 1
    nstates: int = 2

    def __post_init__(self):
        if self.backend not in ("numpy", "torch"):
            raise ValueError("backend must be 'numpy' or 'torch'")

    @property
    def xp(self):
        if self.backend == "torch":
            try:
                import torch
            except ImportError as exc:  # pragma: no cover - optional dependency
                raise ImportError("backend='torch' requires PyTorch") from exc
            required = ("as_tensor", "zeros_like", "exp", "where", "linalg", "diag_embed")
            missing = [name for name in required if not hasattr(torch, name)]
            if missing:
                joined = ", ".join(missing)
                raise ImportError(
                    "backend='torch' requires a usable PyTorch installation; "
                    f"the imported torch module is missing: {joined}"
                )
            return torch
        return np

    def __call__(self, R, P=None, storage=None, traj=None, params=None):
        """
        Evaluate a trajectory batch for ``HamiltonianEngine``.

        ``R`` is expected to have shape ``(ntbf, ndof)``. Only diabatic
        quantities are returned. The Hamiltonian machinery derives adiabatic
        quantities from these fields when needed.
        """

        effective_params = self._merged_params(params)
        Q = self._trajectory_to_exact_q(R)
        return self.evaluate(Q, params=effective_params)

    def evaluate(self, Q, params=None) -> dict[str, Any]:
        """
        Return native diabatic model data for exact-style ``Q``.

        Parameters
        ----------
        Q
            Coordinate tensor with first axis enumerating nuclear degrees of
            freedom. Remaining axes are arbitrary grid or batch dimensions.
        params
            Optional parameter overrides merged with ``self.params``.
        """

        effective_params = self._merged_params(params)
        H_dia, dH_dia = self.diabatic_with_derivatives(Q, effective_params)
        return {
            "H_dia": H_dia,
            "dH_dia": dH_dia,
            "DC1_dia": self.zeros_like(dH_dia),
            "S_dia": self.overlap_like(H_dia),
        }

    def diabatic(self, Q, params=None):
        """Return only ``H_dia`` for exact-dynamics coordinate grids."""

        return self.diabatic_with_derivatives(Q, self._merged_params(params))[0]

    def diabatic_with_derivatives(self, Q, params):
        raise NotImplementedError

    def overlap_like(self, H):
        xp = self.xp
        shape = H.shape[:-2]
        if self.backend == "torch":
            eye = xp.eye(self.nstates, dtype=H.dtype, device=H.device)
            return eye.expand(*shape, self.nstates, self.nstates).clone()
        return np.broadcast_to(np.eye(self.nstates, dtype=H.dtype), (*shape, self.nstates, self.nstates)).copy()

    def zeros_like(self, value):
        return self.xp.zeros_like(value)

    def asarray(self, value):
        if self.backend == "torch":
            return self.xp.as_tensor(value)
        return np.asarray(value)

    def _merged_params(self, params=None):
        merged = dict(self.params)
        if params:
            merged.update(params)
        return merged

    def _trajectory_to_exact_q(self, R):
        xp = self.xp
        arr = self.asarray(R)
        if len(arr.shape) == 1:
            arr = arr.reshape(1, self.ndof)
        return xp.swapaxes(arr, 0, 1)
