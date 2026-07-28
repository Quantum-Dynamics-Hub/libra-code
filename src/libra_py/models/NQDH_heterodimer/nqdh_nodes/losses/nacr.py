"""Loss helpers for nonadiabatic coupling vector training."""

from __future__ import annotations

import itertools
from collections.abc import Sequence

import numpy as np
import torch
from hippynn.graphs import IdxType
from hippynn.graphs.nodes.base import SingleNode


def build_pair_sign_matrix(
    pairs: Sequence[tuple[int, int]],
    states: Sequence[int] | None = None,
) -> np.ndarray:
    """Build all state-sign-consistent pair signs for real electronic states.

    If state ``i`` has sign ``s_i`` and state ``j`` has sign ``s_j``, the NACR
    pair ``d_ij`` transforms as ``s_i * s_j * d_ij``. The global sign is
    redundant, so the first selected state is fixed to ``+1``.
    """
    pair_list = [(int(i), int(j)) for i, j in pairs]
    if not pair_list:
        raise ValueError("At least one NACR pair is required")

    if states is None:
        state_list = sorted({state for pair in pair_list for state in pair})
    else:
        state_list = [int(state) for state in states]

    if not state_list:
        raise ValueError("At least one electronic state is required")
    missing = sorted({state for pair in pair_list for state in pair} - set(state_list))
    if missing:
        raise ValueError(f"Pair states {missing} are missing from the sign-state list")

    anchor = state_list[0]
    free_states = state_list[1:]
    rows: list[list[float]] = []
    for free_signs in itertools.product((-1.0, 1.0), repeat=len(free_states)):
        state_signs = {anchor: 1.0, **dict(zip(free_states, free_signs))}
        rows.append([state_signs[i] * state_signs[j] for i, j in pair_list])

    return np.asarray(rows, dtype=np.float32)


class StateConsistentPhaselessLoss(torch.nn.Module):
    """Phaseless vector loss with one consistent state-sign assignment.

    The input tensors are expected to have shape ``(batch, n_pairs, n_values)``.
    ``pair_indices`` selects which predicted pairs are active in this metric.
    ``true_pair_indices`` can select corresponding pairs from a larger stored
    target tensor.
    ``pair_signs`` has shape ``(n_sign_assignments, n_selected_pairs)``.
    """

    def __init__(
        self,
        pair_indices: Sequence[int],
        pair_signs: np.ndarray,
        mode: str,
        pair_weights: Sequence[float] | None = None,
        true_pair_indices: Sequence[int] | None = None,
        eps: float = 1e-12,
    ) -> None:
        super().__init__()
        if mode not in {"mse", "rmse", "mae"}:
            raise ValueError(f"Unsupported NACR loss mode {mode!r}")

        pair_indices_array = np.asarray(pair_indices, dtype=np.int64)
        if true_pair_indices is None:
            true_pair_indices_array = pair_indices_array
        else:
            true_pair_indices_array = np.asarray(true_pair_indices, dtype=np.int64)
        pair_signs_array = np.asarray(pair_signs, dtype=np.float32)
        if pair_indices_array.ndim != 1 or pair_indices_array.size == 0:
            raise ValueError("pair_indices must be a non-empty 1D sequence")
        if true_pair_indices_array.shape != pair_indices_array.shape:
            raise ValueError(
                f"true_pair_indices must have shape {pair_indices_array.shape}, got {true_pair_indices_array.shape}"
            )
        if pair_signs_array.ndim != 2:
            raise ValueError("pair_signs must have shape (n_assignments, n_pairs)")
        if pair_signs_array.shape[1] != pair_indices_array.size:
            raise ValueError(
                "pair_signs second dimension must match selected pair count: "
                f"{pair_signs_array.shape[1]} != {pair_indices_array.size}"
            )

        if pair_weights is None:
            pair_weights_array = np.ones(pair_indices_array.size, dtype=np.float32)
        else:
            pair_weights_array = np.asarray(pair_weights, dtype=np.float32)
            if pair_weights_array.shape != (pair_indices_array.size,):
                raise ValueError(
                    f"pair_weights must have shape {(pair_indices_array.size,)}, got {pair_weights_array.shape}"
                )
        if np.any(pair_weights_array < 0):
            raise ValueError("pair_weights must be non-negative")
        if float(pair_weights_array.sum()) <= 0:
            raise ValueError("At least one pair weight must be positive")

        self.mode = mode
        self.eps = eps
        self.register_buffer("pair_indices", torch.as_tensor(pair_indices_array, dtype=torch.long))
        self.register_buffer("true_pair_indices", torch.as_tensor(true_pair_indices_array, dtype=torch.long))
        self.register_buffer("pair_signs", torch.as_tensor(pair_signs_array, dtype=torch.float32))
        self.register_buffer("pair_weights", torch.as_tensor(pair_weights_array, dtype=torch.float32))

    def forward(self, predicted: torch.Tensor, true: torch.Tensor) -> torch.Tensor:
        if predicted.ndim != 3 or true.ndim != 3:
            raise ValueError(f"Expected NACR tensors with shape (batch, n_pairs, n_values), got {predicted.shape}")
        if predicted.shape[0] != true.shape[0] or predicted.shape[2] != true.shape[2]:
            raise ValueError(f"Predicted/true NACR target shapes differ: {predicted.shape} != {true.shape}")

        selected_predicted = predicted.index_select(1, self.pair_indices)
        selected_true = true.index_select(1, self.true_pair_indices)
        signs = self.pair_signs.to(dtype=selected_true.dtype).view(1, -1, selected_true.shape[1], 1)
        weights = self.pair_weights.to(dtype=selected_true.dtype).view(1, 1, selected_true.shape[1])

        signed_true = signs * selected_true.unsqueeze(1)
        diff = selected_predicted.unsqueeze(1) - signed_true

        if self.mode in {"mse", "rmse"}:
            per_pair = diff.square().mean(dim=-1)
        else:
            per_pair = diff.abs().mean(dim=-1)

        per_assignment = (per_pair * weights).sum(dim=-1) / weights.sum()
        best_per_sample = per_assignment.min(dim=1).values
        if self.mode == "rmse":
            return torch.sqrt(best_per_sample.mean() + self.eps)
        return best_per_sample.mean()


class StateConsistentPhaselessLossNode(SingleNode):
    """hippy-nn graph node wrapping :class:`StateConsistentPhaselessLoss`."""

    _index_state = IdxType.Scalar

    def __init__(
        self,
        name: str,
        nacr_node,
        pair_indices: Sequence[int],
        pair_signs: np.ndarray,
        mode: str,
        pair_weights: Sequence[float] | None = None,
        true_pair_indices: Sequence[int] | None = None,
    ) -> None:
        module = StateConsistentPhaselessLoss(
            pair_indices=pair_indices,
            pair_signs=pair_signs,
            mode=mode,
            pair_weights=pair_weights,
            true_pair_indices=true_pair_indices,
        )
        parents = (nacr_node.main_output.pred, nacr_node.main_output.true)
        super().__init__(name, parents, module=module)


class ScaledNACRCurvatureLoss(torch.nn.Module):
    """Gauge-invariant curvature loss on the per-pair scaled-NACR norm.

    For pair ``ij`` the curvature (Kubo sum-over-states numerator) is

    ``K_ij = sum_mu (d~_ij^mu)^2 = (E_j - E_i)^2 * sum_mu (d_ij^mu)^2``

    Because ``K_ij`` is quadratic in the scaled coupling, it is invariant under
    electronic-state sign flips (``d_ij -> s_i s_j d_ij`` leaves ``K_ij``
    unchanged), so no phaseless minimization is required. In the energy-ordered
    adiabatic gauge the pair labels are canonical, so this is a direct per-pair
    comparison of predicted and reference curvature.

    Inputs have shape ``(batch, n_pairs, n_values)``. ``pair_indices`` selects
    predicted pairs and ``true_pair_indices`` selects the matching pairs from a
    possibly larger stored target tensor.
    """

    def __init__(
        self,
        pair_indices: Sequence[int],
        mode: str,
        pair_weights: Sequence[float] | None = None,
        true_pair_indices: Sequence[int] | None = None,
        eps: float = 1e-12,
    ) -> None:
        super().__init__()
        if mode not in {"mse", "rmse", "mae"}:
            raise ValueError(f"Unsupported NACR curvature loss mode {mode!r}")

        pair_indices_array = np.asarray(pair_indices, dtype=np.int64)
        if true_pair_indices is None:
            true_pair_indices_array = pair_indices_array
        else:
            true_pair_indices_array = np.asarray(true_pair_indices, dtype=np.int64)
        if pair_indices_array.ndim != 1 or pair_indices_array.size == 0:
            raise ValueError("pair_indices must be a non-empty 1D sequence")
        if true_pair_indices_array.shape != pair_indices_array.shape:
            raise ValueError(
                f"true_pair_indices must have shape {pair_indices_array.shape}, got {true_pair_indices_array.shape}"
            )

        if pair_weights is None:
            pair_weights_array = np.ones(pair_indices_array.size, dtype=np.float32)
        else:
            pair_weights_array = np.asarray(pair_weights, dtype=np.float32)
            if pair_weights_array.shape != (pair_indices_array.size,):
                raise ValueError(
                    f"pair_weights must have shape {(pair_indices_array.size,)}, got {pair_weights_array.shape}"
                )
        if np.any(pair_weights_array < 0):
            raise ValueError("pair_weights must be non-negative")
        if float(pair_weights_array.sum()) <= 0:
            raise ValueError("At least one pair weight must be positive")

        self.mode = mode
        self.eps = eps
        self.register_buffer("pair_indices", torch.as_tensor(pair_indices_array, dtype=torch.long))
        self.register_buffer("true_pair_indices", torch.as_tensor(true_pair_indices_array, dtype=torch.long))
        self.register_buffer("pair_weights", torch.as_tensor(pair_weights_array, dtype=torch.float32))

    def forward(self, predicted: torch.Tensor, true: torch.Tensor) -> torch.Tensor:
        if predicted.ndim != 3 or true.ndim != 3:
            raise ValueError(f"Expected NACR tensors with shape (batch, n_pairs, n_values), got {predicted.shape}")
        if predicted.shape[0] != true.shape[0] or predicted.shape[2] != true.shape[2]:
            raise ValueError(f"Predicted/true NACR target shapes differ: {predicted.shape} != {true.shape}")

        selected_predicted = predicted.index_select(1, self.pair_indices)
        selected_true = true.index_select(1, self.true_pair_indices)
        weights = self.pair_weights.to(dtype=selected_true.dtype).view(1, -1)

        curvature_pred = selected_predicted.square().sum(dim=-1)
        curvature_true = selected_true.square().sum(dim=-1)
        diff = curvature_pred - curvature_true

        if self.mode in {"mse", "rmse"}:
            per_pair = diff.square()
        else:
            per_pair = diff.abs()

        per_sample = (per_pair * weights).sum(dim=-1) / weights.sum()
        if self.mode == "rmse":
            return torch.sqrt(per_sample.mean() + self.eps)
        return per_sample.mean()


class ScaledNACRCurvatureLossNode(SingleNode):
    """hippy-nn graph node wrapping :class:`ScaledNACRCurvatureLoss`."""

    _index_state = IdxType.Scalar

    def __init__(
        self,
        name: str,
        nacr_node,
        pair_indices: Sequence[int],
        mode: str,
        pair_weights: Sequence[float] | None = None,
        true_pair_indices: Sequence[int] | None = None,
    ) -> None:
        module = ScaledNACRCurvatureLoss(
            pair_indices=pair_indices,
            mode=mode,
            pair_weights=pair_weights,
            true_pair_indices=true_pair_indices,
        )
        parents = (nacr_node.main_output.pred, nacr_node.main_output.true)
        super().__init__(name, parents, module=module)
