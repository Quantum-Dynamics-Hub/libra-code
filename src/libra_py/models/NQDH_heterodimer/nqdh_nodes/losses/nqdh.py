"""Losses for the Neural Quasi-Diabatic Hamiltonian (NQDH) model.

Two terms specific to NQDH (gradient and coupling losses reuse the existing
adiabatic / scaled-NACR losses on the diagonal/off-diagonal of ``M``):

  PowerSumEnergyLoss   compares the predicted power sums ``p_k = Tr(W^k)`` to the
                       power sums of the reference energies, ``sum_i (E_i^ref)^k``.
                       Label-free (depends only on the unordered spectrum) and
                       needs no eigendecomposition, so it is smooth through
                       interior conical intersections.

  SmoothnessLoss       penalizes ``||dW/dR||^2`` (Frobenius, over nuclear DOF).
                       W's smoothness is the physical content of the model, so it
                       is imposed rather than assumed.

See ``docs/nqdh_model.md``.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import torch
from hippynn.graphs import IdxType
from hippynn.graphs.nodes.base import SingleNode

from .nacr import StateConsistentPhaselessLoss


def reference_power_sums(energies: torch.Tensor, n_states: int) -> torch.Tensor:
    """Power sums ``p_k = sum_i E_i^k`` for ``k = 1..K`` from reference energies.

    ``energies`` has shape ``(batch, K)`` (the K band-state reference energies, in
    any order; the result is order-independent).
    """
    powers = [energies.pow(k).sum(dim=-1) for k in range(1, n_states + 1)]
    return torch.stack(powers, dim=-1)


class PowerSumEnergyLoss(torch.nn.Module):
    """Energy loss in power-sum (characteristic-polynomial) space.

    Compares predicted power sums ``p_k(W)`` to the power sums of the reference
    energies. Each ``k`` is normalized by ``energy_scale**k`` so that, for energies
    of magnitude ~``energy_scale``, all ``K`` residuals are O(1) and no single high
    order dominates. The reductions match the other losses in this repo
    (``mse`` / ``rmse`` / ``mae``).

    Inputs to ``forward``:
        predicted_power_sums: ``(batch, K)``  =  ``Tr(W^k)``, k=1..K
        reference_energies:   ``(batch, K)``  the K band reference energies
    """

    def __init__(self, n_states: int, energy_scale: float = 1.0, mode: str = "mse",
                 eps: float = 1e-12, max_k: int | None = None) -> None:
        super().__init__()
        if n_states < 1:
            raise ValueError(f"n_states must be >= 1, got {n_states}")
        if mode not in {"mse", "rmse", "mae"}:
            raise ValueError(f"Unsupported energy loss mode {mode!r}")
        if energy_scale <= 0:
            raise ValueError(f"energy_scale must be positive, got {energy_scale}")
        self.n_states = n_states
        self.mode = mode
        self.eps = eps
        # max_k caps the power-sum order used in the loss. The high orders Tr(W^k) of a
        # K-state band have an enormous dynamic range (dominated by the largest |E|) and
        # are ill-conditioned at large K, while the low orders are the well-conditioned
        # spectral moments (k=1 sum, k=2 spread, k=3 skew) -- smooth through degeneracies
        # and a stable spectral backbone alongside the per-state eigenvalue loss. max_k=K
        # (default) keeps the full characteristic-polynomial readout.
        self.max_k = n_states if max_k is None else int(max_k)
        if not (1 <= self.max_k <= n_states):
            raise ValueError(f"max_k must be in [1, {n_states}], got {self.max_k}")
        # per-k normalization 1/scale^k keeps Tr(W^k) residuals comparable across k
        inv = np.array([energy_scale ** (-k) for k in range(1, self.max_k + 1)], dtype=np.float32)
        self.register_buffer("inv_scale_k", torch.as_tensor(inv))

    def forward(self, predicted_power_sums: torch.Tensor, *reference_energy_columns: torch.Tensor) -> torch.Tensor:
        if predicted_power_sums.shape[-1] != self.n_states:
            raise ValueError(f"predicted_power_sums last dim {predicted_power_sums.shape[-1]} != K {self.n_states}")
        # reference energies arrive as either one (batch, K) tensor or K separate
        # (batch, 1)/(batch,) database columns (the per-state E{s} .true values).
        if len(reference_energy_columns) == 1 and reference_energy_columns[0].shape[-1] == self.n_states:
            reference_energies = reference_energy_columns[0]
        else:
            if len(reference_energy_columns) != self.n_states:
                raise ValueError(
                    f"expected {self.n_states} reference energy columns, got {len(reference_energy_columns)}"
                )
            reference_energies = torch.cat([c.reshape(c.shape[0], 1) for c in reference_energy_columns], dim=-1)

        ref_pk = reference_power_sums(reference_energies, self.n_states)
        # use only the first max_k orders (the predicted high orders, if computed, get
        # no gradient -- the ill-conditioned tail is dropped from the objective).
        pred_k = predicted_power_sums[..., : self.max_k]
        ref_k = ref_pk[..., : self.max_k]
        scale = self.inv_scale_k.to(dtype=predicted_power_sums.dtype)
        diff = (pred_k - ref_k) * scale

        if self.mode in {"mse", "rmse"}:
            per_sample = diff.square().mean(dim=-1)
        else:
            per_sample = diff.abs().mean(dim=-1)

        if self.mode == "rmse":
            return torch.sqrt(per_sample.mean() + self.eps)
        return per_sample.mean()


class WRoughness(torch.nn.Module):
    """Per-molecule spatial roughness of ``W``: ``r_b = ||dW_b/dR||_F^2``.

    This is a MODEL-side module (it needs ``positions`` and autograd), analogous to
    how the eigenframe gradient lives in the model graph. It returns a
    ``(batch,)`` tensor; the actual penalty is the mean of this, applied by a
    trivial reduction loss (``MeanScalarLoss``) so the loss graph never holds a
    raw input node.

    Computes ``dW/dR`` for every matrix entry by autograd (one grad per upper-tri
    entry, mirrored), and sums the squared derivatives over all entries and
    nuclear DOF. ``positions`` must have ``requires_grad``.
    """

    def __init__(self) -> None:
        super().__init__()

    def forward(self, W: torch.Tensor, positions: torch.Tensor) -> torch.Tensor:
        if W.ndim != 3 or W.shape[-1] != W.shape[-2]:
            raise ValueError(f"Expected W shape (batch, K, K), got {tuple(W.shape)}")
        if not positions.requires_grad:
            raise ValueError("positions must have requires_grad=True for WRoughness")
        n_batch, K, _ = W.shape
        if not W.requires_grad:
            return W.new_zeros(n_batch)
        ones = torch.ones(n_batch, dtype=W.dtype, device=W.device)
        per_sample = W.new_zeros(n_batch)
        for m in range(K):
            for n in range(m, K):
                g = torch.autograd.grad(
                    W[:, m, n], positions, grad_outputs=ones,
                    create_graph=True, retain_graph=True, allow_unused=True,
                )[0]
                if g is None:
                    continue  # this entry is constant in R -> contributes 0
                g = g.reshape(n_batch, -1)
                sq = g.square().sum(dim=-1)            # per-sample ||dW_mn/dR||^2
                per_sample = per_sample + (sq if m == n else 2.0 * sq)
        return per_sample


# Backwards-compatible alias used by the unit tests (returns the mean roughness).
class SmoothnessLoss(torch.nn.Module):
    """Mean over the batch of :class:`WRoughness` -> scalar ``||dW/dR||_F^2``."""

    def __init__(self) -> None:
        super().__init__()
        self.roughness = WRoughness()

    def forward(self, W: torch.Tensor, positions: torch.Tensor) -> torch.Tensor:
        if not W.requires_grad:
            return W.new_zeros(())
        return self.roughness(W, positions).mean()


class MeanScalar(torch.nn.Module):
    """Reduce a per-molecule tensor ``(batch,)`` to its mean (a regularizer loss)."""

    def forward(self, per_sample: torch.Tensor) -> torch.Tensor:
        return per_sample.mean()


# --------------------------------------------------------------------------- #
# hippy-nn graph-node wrappers
# --------------------------------------------------------------------------- #


class PowerSumEnergyLossNode(SingleNode):
    """Graph node wrapping :class:`PowerSumEnergyLoss`.

    Parents: ``(power_sums_node, reference_energy_node)`` where the reference
    energy node carries the ``(batch, K)`` band reference energies (e.g. a
    MoleculeCatNode over the per-state energy db targets, or a single stacked
    energy db target).
    """

    _index_state = IdxType.Scalar

    def __init__(
        self,
        name: str,
        power_sums_node,
        reference_energy_nodes,
        n_states: int,
        energy_scale: float = 1.0,
        mode: str = "mse",
        max_k: int | None = None,
    ) -> None:
        module = PowerSumEnergyLoss(n_states, energy_scale=energy_scale, mode=mode, max_k=max_k)
        # predicted power sums (.pred) + the K reference energies straight from the
        # database (.true), so the reference is never recomputed by the model.
        ref_parents = tuple(node.main_output.true for node in reference_energy_nodes)
        parents = (power_sums_node.main_output.pred, *ref_parents)
        super().__init__(name, parents, module=module)


class WRoughnessNode(SingleNode):
    """MODEL-side node: (W, positions) -> per-molecule roughness ``(batch,)``.

    Lives in the model graph (it holds the ``dW/dR`` autograd and needs the raw
    positions input), so it is paired with :class:`MeanScalarLossNode` in the loss
    graph rather than depending on an input there directly.
    """

    _index_state = IdxType.Molecules

    def __init__(self, name: str, W_node, positions, **kwargs) -> None:
        positions.requires_grad = True
        super().__init__(name, (W_node.main_output, positions), module=WRoughness(), **kwargs)


class MeanScalarLossNode(SingleNode):
    """Loss-graph node: reduce a per-molecule model output to its batch mean."""

    _index_state = IdxType.Scalar

    def __init__(self, name: str, roughness_node, **kwargs) -> None:
        super().__init__(name, (roughness_node.main_output.pred,), module=MeanScalar(), **kwargs)


class WFrobeniusAnchorLoss(torch.nn.Module):
    """Eigh-free derivative-magnitude anchor.

    Compares the model's per-sample ``||dW/dR||_F^2`` (``predicted``) to the
    reference ``sum_i ||G_i||^2 + 2 sum_{i<j} ||q_ij||^2`` assembled from the
    database gradients and scaled couplings. Gauge/label invariant, needs no
    eigendecomposition, smooth through CIs -> anchors W's derivative magnitude even
    at special points where the individual eigenframe entries are unsupervised.

    ``forward(predicted, *refs)`` where ``refs`` are the K gradient ``.true``
    tensors (each ``(batch, n_atoms, 3)``) followed by the coupling ``.true``
    tensor ``(batch, n_pairs_total, 3N)``. ``coupling_pair_indices`` selects which
    stored pairs are summed (the supervised/active pairs); if None, all are used.
    Normalized by ``scale`` to keep it comparable to the other losses.
    """

    def __init__(self, n_grad: int, coupling_pair_indices: Sequence[int] | None = None,
                 scale: float = 1.0, mode: str = "mse", eps: float = 1e-12) -> None:
        super().__init__()
        if mode not in {"mse", "rmse", "mae"}:
            raise ValueError(f"Unsupported anchor loss mode {mode!r}")
        self.n_grad = n_grad
        self.mode = mode
        self.scale = float(scale)
        self.eps = eps
        if coupling_pair_indices is not None:
            self.register_buffer("pair_idx", torch.as_tensor(list(coupling_pair_indices), dtype=torch.long))
        else:
            self.pair_idx = None

    def forward(self, predicted: torch.Tensor, *refs: torch.Tensor) -> torch.Tensor:
        grads = refs[: self.n_grad]
        ref = torch.zeros_like(predicted)
        for g in grads:
            ref = ref + g.reshape(g.shape[0], -1).square().sum(dim=-1)
        if len(refs) > self.n_grad:
            coupling = refs[self.n_grad]                      # (batch, n_pairs_total, 3N)
            if self.pair_idx is not None:
                coupling = coupling.index_select(1, self.pair_idx.to(coupling.device))
            ref = ref + 2.0 * coupling.reshape(coupling.shape[0], coupling.shape[1], -1).square().sum(dim=(1, 2))
        diff = (predicted - ref) / (self.scale * self.scale)
        per_sample = diff.abs() if self.mode == "mae" else diff.square()
        out = per_sample.mean()
        return torch.sqrt(out + self.eps) if self.mode == "rmse" else out


class WFrobeniusAnchorLossNode(SingleNode):
    """Loss node: model ||dW||_F^2 vs the data-derived reference.

    Parents: model anchor node (.pred) + the K gradient nodes (.true) + the
    coupling node (.true).
    """

    _index_state = IdxType.Scalar

    def __init__(self, name: str, anchor_node, gradient_nodes, coupling_node,
                 coupling_pair_indices=None, scale: float = 1.0, mode: str = "mse") -> None:
        module = WFrobeniusAnchorLoss(
            n_grad=len(gradient_nodes), coupling_pair_indices=coupling_pair_indices,
            scale=scale, mode=mode,
        )
        parents = (
            anchor_node.main_output.pred,
            *(gn.main_output.true for gn in gradient_nodes),
            coupling_node.main_output.true,
        )
        super().__init__(name, parents, module=module)


# --------------------------------------------------------------------------- #
# Eigen-gap-masked eigenframe losses
# --------------------------------------------------------------------------- #
#
# The gradient and coupling losses depend on the eigenframe U (from eigh of W).
# Near a predicted eigenvalue degeneracy the eigenvector gradient blows up
# (1/(lambda_i - lambda_j)), which destabilizes training. These variants accept a
# per-sample detached weight (the eigen-gap mask) and down-weight such samples.


class MaskedGradientLoss(torch.nn.Module):
    """Per-sample-weighted gradient loss: weighted mean of ``||g_pred - g_true||``.

    Inputs:
        predicted: ``(batch, n_atoms, 3)``   diag(M) for one state
        true:      ``(batch, n_atoms, 3)``   stored gradient
        weight:    ``(batch,)``              detached eigen-gap mask
    """

    def __init__(self, mode: str = "mse", weight_index: int | None = None, eps: float = 1e-12) -> None:
        super().__init__()
        if mode not in {"mse", "rmse", "mae"}:
            raise ValueError(f"Unsupported gradient loss mode {mode!r}")
        self.mode = mode
        self.weight_index = weight_index   # column of a per-state mask; None => mask is (batch,)
        self.eps = eps

    def forward(self, predicted: torch.Tensor, true: torch.Tensor, weight: torch.Tensor) -> torch.Tensor:
        diff = predicted - true
        per_sample = diff.reshape(diff.shape[0], -1)
        per_sample = per_sample.square().mean(-1) if self.mode in {"mse", "rmse"} else per_sample.abs().mean(-1)
        w = weight.detach().to(dtype=per_sample.dtype)
        if w.ndim == 2:
            w = w[:, self.weight_index] if self.weight_index is not None else w.min(dim=-1).values
        out = (w * per_sample).sum() / w.sum().clamp_min(1e-12)
        return torch.sqrt(out + self.eps) if self.mode == "rmse" else out


class MaskedPhaselessCouplingLoss(StateConsistentPhaselessLoss):
    """:class:`StateConsistentPhaselessLoss` with a per-(sample, pair) eigen-gap weight.

    Reuses the full sign-consistency machinery of the base class, but weights each
    pair's per-sample error by its detached eigen-gap mask so that a near-degenerate
    pair contributes ~0 (its frame is unsupervised). ``forward(predicted, true,
    pair_weight)`` where ``pair_weight`` is ``(batch, n_selected_pairs)`` aligned
    with ``pair_indices`` (the supervised pairs). A scalar/``(batch,)`` weight is
    also accepted and broadcast across pairs.
    """

    def forward(self, predicted: torch.Tensor, true: torch.Tensor, pair_weight: torch.Tensor) -> torch.Tensor:
        selected_predicted = predicted.index_select(1, self.pair_indices)
        selected_true = true.index_select(1, self.true_pair_indices)
        n_pairs = selected_true.shape[1]
        signs = self.pair_signs.to(dtype=selected_true.dtype).view(1, -1, n_pairs, 1)
        base_w = self.pair_weights.to(dtype=selected_true.dtype).view(1, 1, n_pairs)

        gate = pair_weight.detach().to(dtype=selected_true.dtype)
        mask_idx = getattr(self, "mask_pair_indices", None)
        if gate.ndim == 2 and mask_idx is not None:
            gate = gate.index_select(1, mask_idx.to(gate.device))   # select supervised pairs
        if gate.ndim == 1:
            gate = gate.unsqueeze(-1).expand(-1, n_pairs)
        gate = gate.view(gate.shape[0], 1, n_pairs)          # (batch, 1, n_pairs)
        eff_w = base_w * gate                                 # gap-weighted per pair

        signed_true = signs * selected_true.unsqueeze(1)
        diff = selected_predicted.unsqueeze(1) - signed_true
        per_pair = diff.square().mean(dim=-1) if self.mode in {"mse", "rmse"} else diff.abs().mean(dim=-1)
        # weighted mean over pairs; denominator is the per-sample effective weight
        denom = eff_w.sum(dim=-1).clamp_min(1e-12)            # (batch, n_assign)
        per_assignment = (per_pair * eff_w).sum(dim=-1) / denom
        best_per_sample = per_assignment.min(dim=1).values

        # also down-weight whole samples by their total gate mass (so a sample with
        # all pairs degenerate contributes ~0 rather than being renormalized up)
        sample_w = gate.view(gate.shape[0], n_pairs).mean(dim=-1)
        out = (sample_w * best_per_sample).sum() / sample_w.sum().clamp_min(1e-12)
        return torch.sqrt(out + self.eps) if self.mode == "rmse" else out


class MaskedGradientLossNode(SingleNode):
    """Loss node: (pred grad, true grad, mask) -> masked gradient loss.

    ``weight_index`` selects this state's column from a per-state mask
    ``(batch, K)``; pass None if the mask is a plain ``(batch,)``.
    """

    _index_state = IdxType.Scalar

    def __init__(self, name: str, grad_node, mask_node, mode: str = "mse", weight_index: int | None = None) -> None:
        parents = (grad_node.main_output.pred, grad_node.main_output.true, mask_node.main_output.pred)
        super().__init__(name, parents, module=MaskedGradientLoss(mode=mode, weight_index=weight_index))


class CharacterLoss(torch.nn.Module):
    """MSE on diabatic characters: pred/true ``(batch, 2, 2)``.

    ``permutation_min=False`` (default): plain MSE. This is what PINS the model's
    internal diabat assignment to the fragment labels -- the whole point of the
    loss. Safe at degeneracies without any gating: both rows -> 1/2 there, so a
    state-order mismatch costs ~0 (self-quenching).

    ``permutation_min=True``: min over reference row (state) permutations -- the
    orbit-space form for systems with nontrivial eigenbundle holonomy where NO
    global assignment exists. CAUTION: for a 2-state complementary block
    ([[p,1-p],[1-p,p]]) the row swap is IDENTICAL to the fragment (column) swap,
    so this form is structurally blind to the A/B assignment and cannot pin it at
    any weight (measured: sign agreement stayed 48-52% at weight 10 and 200).
    Use it only where the assignment genuinely must stay unpinned.

    NaN-labeled samples are masked out. Deliberately NOT gap-masked: the
    information lives at the seam; the character node's eps regularizer is the
    numerical guard instead.
    """

    def __init__(self, permutation_min: bool = False) -> None:
        super().__init__()
        self.permutation_min = bool(permutation_min)

    def forward(self, predicted: torch.Tensor, true: torch.Tensor) -> torch.Tensor:
        finite = torch.isfinite(true).all(dim=-1).all(dim=-1)
        safe = torch.where(torch.isfinite(true), true, torch.zeros_like(true))
        e_id = (predicted - safe).square().mean(dim=(-1, -2))
        if self.permutation_min:
            e_sw = (predicted - safe.flip(-2)).square().mean(dim=(-1, -2))  # swap reference rows
            e_id = torch.minimum(e_id, e_sw)
        per_sample = e_id * finite.to(predicted.dtype)
        return per_sample.sum() / finite.sum().clamp_min(1)


class CharacterLossNode(SingleNode):
    """Loss node: (predicted characters, stored characters) -> permutation-min MSE."""

    _index_state = IdxType.Scalar

    def __init__(self, name: str, char_node) -> None:
        parents = (char_node.main_output.pred, char_node.main_output.true)
        super().__init__(name, parents, module=CharacterLoss())


class MaskedPhaselessCouplingLossNode(SingleNode):
    """Loss node: (pred coupling, true coupling, per-pair mask) -> masked loss.

    ``mask_node`` carries the per-pair gap mask ``(batch, n_pairs_total)``;
    ``mask_pair_indices`` selects the columns for the supervised pairs (aligned
    with ``pair_indices``).
    """

    _index_state = IdxType.Scalar

    def __init__(
        self, name: str, coupling_node, mask_node, pair_indices: Sequence[int],
        pair_signs: np.ndarray, mode: str, true_pair_indices: Sequence[int] | None = None,
        mask_pair_indices: Sequence[int] | None = None,
    ) -> None:
        module = MaskedPhaselessCouplingLoss(
            pair_indices=pair_indices, pair_signs=pair_signs, mode=mode,
            true_pair_indices=true_pair_indices,
        )
        module.mask_pair_indices = (
            None if mask_pair_indices is None
            else torch.as_tensor(list(mask_pair_indices), dtype=torch.long)
        )
        parents = (coupling_node.main_output.pred, coupling_node.main_output.true, mask_node.main_output.pred)
        super().__init__(name, parents, module=module)


class DecodedEnergyError(torch.nn.Module):
    """Energy-ordered decoded-energy error (eV) for monitoring.

    Compares the model's decoded spectrum ``eigvalsh(W)`` (already ascending) to the
    reference energies SORTED ascending, then reports the MAE/RMSE for one
    ``state_index`` (or, with ``state_index=None``, the mean over all states). Both
    sides are energy-ordered, so this is robust to the reference state labels
    swapping near degeneracies -- it answers "how well does the model spectrum match
    the electronic-structure spectrum", which is the right 1 eV / 1 meV check.

    ``forward(eigenvalues, *reference_energy_columns)`` where ``eigenvalues`` is
    ``(batch, K)`` ascending and the references are K ``(batch, 1)`` database
    columns. Result is in eV.
    """

    def __init__(self, n_states: int, state_index: int | None = None, mode: str = "mae") -> None:
        super().__init__()
        if mode not in {"mae", "rmse"}:
            raise ValueError(f"Unsupported decoded-energy mode {mode!r}")
        self.n_states = n_states
        self.state_index = state_index
        self.mode = mode

    def forward(self, eigenvalues: torch.Tensor, *reference_energy_columns: torch.Tensor) -> torch.Tensor:
        ref = torch.cat([c.reshape(c.shape[0], 1) for c in reference_energy_columns], dim=-1)
        ref_sorted, _ = torch.sort(ref, dim=-1)              # energy-order the references
        diff = eigenvalues - ref_sorted                       # both ascending -> aligned
        if self.state_index is not None:
            diff = diff[:, self.state_index]
        if self.mode == "rmse":
            return torch.sqrt(diff.square().mean())
        return diff.abs().mean()


class DecodedEnergyErrorNode(SingleNode):
    """Monitoring loss node: decoded-energy error (eV), references sorted ascending.

    Parents: the eigenvalues node (.pred) + the K reference energy nodes (.true).
    """

    _index_state = IdxType.Scalar

    def __init__(self, name: str, eigenvalues_node, reference_energy_nodes,
                 n_states: int, state_index: int | None = None, mode: str = "mae") -> None:
        module = DecodedEnergyError(n_states, state_index=state_index, mode=mode)
        ref_parents = tuple(node.main_output.true for node in reference_energy_nodes)
        parents = (eigenvalues_node.main_output.pred, *ref_parents)
        super().__init__(name, parents, module=module)


class EigenvalueEnergyLoss(torch.nn.Module):
    """Direct, state-resolved eigenvalue energy loss (trained term).

    Compares the energy-ordered model spectrum ``sort(eigvalsh(W))`` to the
    energy-ordered references, per state, in eV. Unlike the power-sum loss (which
    only constrains symmetric functions of the spectrum, so per-state error is
    shared and the highest/most-coupled state lags), this puts each ``E_i``
    DIRECTLY in the objective -- the term needed for state-resolved 1 meV.

    Stability: only the *eigenvalue* gradient is used here (``d lambda_i = v_i^T dW
    v_i``, bounded), NOT the eigenvector gradient (the ``1/(lambda_i - lambda_j)``
    blow-up). So this is safe at near-degeneracies and needs no mask for stability.
    The optional per-sample ``mask`` (detached eigen-gap "min" gate) is used only to
    soften the *sorting kink* at exact crossings -- there the power-sum term (kink
    free, smooth through CIs) carries the seam. Pass ``use_mask=False`` to disable.

    ``forward(eigenvalues, *refs[, mask])``: ``eigenvalues`` is ``(batch, K)``
    ascending; ``refs`` are K ``(batch, 1)`` database columns; if ``use_mask`` the
    LAST positional arg is the ``(batch,)`` gate.
    """

    def __init__(self, n_states: int, mode: str = "mse", use_mask: bool = True,
                 state_weights: Sequence[float] | None = None, eps: float = 1e-12) -> None:
        super().__init__()
        if mode not in {"mse", "rmse", "mae"}:
            raise ValueError(f"Unsupported eigenvalue energy loss mode {mode!r}")
        self.n_states = n_states
        self.mode = mode
        self.use_mask = use_mask
        self.eps = eps
        if state_weights is None:
            sw = torch.ones(n_states)
        else:
            sw = torch.as_tensor(list(state_weights), dtype=torch.float32)
            if sw.shape != (n_states,):
                raise ValueError(f"state_weights must have shape ({n_states},), got {tuple(sw.shape)}")
        self.register_buffer("state_weights", sw)

    def forward(self, eigenvalues: torch.Tensor, *rest: torch.Tensor) -> torch.Tensor:
        if self.use_mask:
            *refs, mask = rest
        else:
            refs, mask = rest, None
        ref = torch.cat([c.reshape(c.shape[0], 1) for c in refs], dim=-1)
        ref_sorted, _ = torch.sort(ref, dim=-1)
        diff = eigenvalues - ref_sorted                       # (batch, K), both ascending
        sw = self.state_weights.to(dtype=diff.dtype)
        per_state = diff.square() if self.mode in {"mse", "rmse"} else diff.abs()
        per_sample = (per_state * sw).sum(dim=-1) / sw.sum()  # (batch,)
        if mask is not None:
            w = mask.detach().to(dtype=per_sample.dtype)
            out = (w * per_sample).sum() / w.sum().clamp_min(1e-12)
        else:
            out = per_sample.mean()
        return torch.sqrt(out + self.eps) if self.mode == "rmse" else out


class EigenvalueEnergyLossNode(SingleNode):
    """Trained loss node: direct state-resolved eigenvalue energy loss.

    Parents: eigenvalues node (.pred) + K reference energy nodes (.true)
    [+ eigen-gap "min" mask node (.pred) if ``use_mask``].
    """

    _index_state = IdxType.Scalar

    def __init__(self, name: str, eigenvalues_node, reference_energy_nodes, n_states: int,
                 mode: str = "mse", mask_node=None, state_weights: Sequence[float] | None = None) -> None:
        use_mask = mask_node is not None
        module = EigenvalueEnergyLoss(n_states, mode=mode, use_mask=use_mask, state_weights=state_weights)
        ref_parents = tuple(node.main_output.true for node in reference_energy_nodes)
        parents = (eigenvalues_node.main_output.pred, *ref_parents)
        if use_mask:
            parents = (*parents, mask_node.main_output.pred)
        super().__init__(name, parents, module=module)


class EigenvalueGapLoss(torch.nn.Module):
    """Direct eigenvalue-GAP loss: differential spectroscopy in the objective.

    Compares energy-ordered spectrum DIFFERENCES ``lambda_j - lambda_i`` to the
    reference gaps ``E_j - E_i`` (references sorted ascending) for a set of state
    pairs. The absolute-eigenvalue loss constrains the gaps only indirectly --
    per-state errors that shift states rigidly are invisible to it as far as
    differences go, and (with excited-state down-weighting) the optimizer can
    park at a large gap error while every absolute term looks converged. This
    term puts the *reported* observables -- excitation energies (0-1, 0-2) and
    the dynamics gap (1-2) -- directly in the objective. Differential targets
    also cancel common-mode model/label error (the same mechanism that makes
    gaps the naturally-accurate output of a shared-backbone model).

    Stability: identical to ``EigenvalueEnergyLoss`` -- only bounded eigenvalue
    gradients (``d lambda_i = v_i^T dW v_i``), no eigenvector term, so safe at
    near-degeneracies. The optional detached min-gap ``mask`` softens the
    sorting kink at exact crossings, where the smooth power-sum term carries
    the seam.

    ``forward(eigenvalues, *refs[, mask])``: ``eigenvalues`` is ``(batch, K)``
    ascending; refs are K ``(batch, 1)`` database columns; if ``use_mask`` the
    last positional arg is the ``(batch,)`` gate. Result in eV.
    """

    def __init__(self, n_states: int, pairs: Sequence[tuple[int, int]] | None = None,
                 mode: str = "mae", use_mask: bool = True, eps: float = 1e-12) -> None:
        super().__init__()
        if mode not in {"mse", "rmse", "mae"}:
            raise ValueError(f"Unsupported eigenvalue gap loss mode {mode!r}")
        if pairs is None:
            pairs = [(i, j) for i in range(n_states) for j in range(i + 1, n_states)]
        for i, j in pairs:
            if not (0 <= i < j < n_states):
                raise ValueError(f"Invalid gap pair ({i}, {j}) for n_states={n_states}")
        self.n_states = n_states
        self.pairs = list(pairs)
        self.mode = mode
        self.use_mask = use_mask
        self.eps = eps

    def forward(self, eigenvalues: torch.Tensor, *rest: torch.Tensor) -> torch.Tensor:
        if self.use_mask:
            *refs, mask = rest
        else:
            refs, mask = rest, None
        ref = torch.cat([c.reshape(c.shape[0], 1) for c in refs], dim=-1)
        ref_sorted, _ = torch.sort(ref, dim=-1)
        ii = torch.tensor([p[0] for p in self.pairs], device=eigenvalues.device)
        jj = torch.tensor([p[1] for p in self.pairs], device=eigenvalues.device)
        gap_pred = eigenvalues[:, jj] - eigenvalues[:, ii]     # (batch, n_pairs)
        gap_ref = ref_sorted[:, jj] - ref_sorted[:, ii]
        diff = gap_pred - gap_ref
        per_pair = diff.square() if self.mode in {"mse", "rmse"} else diff.abs()
        per_sample = per_pair.mean(dim=-1)                     # (batch,)
        if mask is not None:
            w = mask.detach().to(dtype=per_sample.dtype)
            out = (w * per_sample).sum() / w.sum().clamp_min(1e-12)
        else:
            out = per_sample.mean()
        return torch.sqrt(out + self.eps) if self.mode == "rmse" else out


class EigenvalueGapLossNode(SingleNode):
    """Loss node: eigenvalue-gap loss (trained with mask_node, or unmasked monitor).

    Parents: eigenvalues node (.pred) + K reference energy nodes (.true)
    [+ eigen-gap "min" mask node (.pred) if given].
    """

    _index_state = IdxType.Scalar

    def __init__(self, name: str, eigenvalues_node, reference_energy_nodes, n_states: int,
                 pairs: Sequence[tuple[int, int]] | None = None, mode: str = "mae",
                 mask_node=None) -> None:
        use_mask = mask_node is not None
        module = EigenvalueGapLoss(n_states, pairs=pairs, mode=mode, use_mask=use_mask)
        ref_parents = tuple(node.main_output.true for node in reference_energy_nodes)
        parents = (eigenvalues_node.main_output.pred, *ref_parents)
        if use_mask:
            parents = (*parents, mask_node.main_output.pred)
        super().__init__(name, parents, module=module)
