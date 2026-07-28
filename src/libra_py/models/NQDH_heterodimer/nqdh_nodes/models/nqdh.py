"""Core modules and graph nodes for the Neural Quasi-Diabatic Hamiltonian (NQDH).

NQDH learns a smooth, single-valued, symmetric matrix field ``W(R)`` whose
eigen-decomposition gives the excited-state observables, instead of learning the
adiabatic surfaces and couplings directly. See ``docs/nqdh_model.md`` for the full
design; the short version:

    W(R) = D(R) + V(R)                          symmetric K x K
    energies      E_i = eigenvalues(W)
    M            = U^T (dW/dR) U                 dW rotated into the eigenframe
    gradients     g_i = M_ii                     (diagonal of M)
    scaled NACR   q_ij = M_ij                     (off-diagonal of M)

This module provides three plain ``torch.nn.Module`` building blocks (testable
without the hippy-nn graph) and thin graph-node wrappers around them:

    PackSymmetric / PackSymmetricNode      diag + offdiag  ->  symmetric W
    PowerSums      / PowerSumsNode          W  ->  (Tr(W^1), ..., Tr(W^K))
    EigenframeGrad / EigenframeGradNode     W, R  ->  M = U^T dW U

The energy readout (PowerSums) deliberately avoids any eigendecomposition: the
power sums ``p_k = Tr(W^k)`` are smooth polynomials in the entries of ``W`` and
are in bijection with the unordered spectrum (Newton's identities), so the energy
loss is label-free and smooth through interior conical intersections. Only
``EigenframeGrad`` diagonalizes, and only the gradient/coupling losses depend on
it.
"""

from __future__ import annotations

import warnings

import numpy as np
import torch
from torch import Tensor

from hippynn.graphs import IdxType
from hippynn.graphs.nodes.base import SingleNode


def _safe_eigvalsh(W: Tensor) -> Tensor:
    """``eigvalsh`` robust to clustered/degenerate spectra.

    The symmetry-adapted W carries exact and near-exact eigenvalue degeneracies by
    construction (symmetry-forced E-blocks), on which cuSOLVER's float32 ``syevd``
    intermittently fails to converge (``LinAlgError`` "ill-conditioned or too many
    repeated eigenvalues"). The matrices are tiny (K x K), so run the decomposition in
    float64 -- far more robust -- and fall back to backward-stable CPU LAPACK if the
    GPU solver still raises. The float/device casts are autograd-safe.
    """
    try:
        return torch.linalg.eigvalsh(W.double()).to(W.dtype)
    except torch._C._LinAlgError:
        return torch.linalg.eigvalsh(W.double().cpu()).to(device=W.device, dtype=W.dtype)


def _safe_eigh(W: Tensor) -> tuple[Tensor, Tensor]:
    """Eigenvalue/eigenvector counterpart of :func:`_safe_eigvalsh`."""
    try:
        vals, vecs = torch.linalg.eigh(W.double())
    except torch._C._LinAlgError:
        vals, vecs = torch.linalg.eigh(W.double().cpu())
        vals, vecs = vals.to(W.device), vecs.to(W.device)
    return vals.to(W.dtype), vecs.to(W.dtype)


class SymmetryAdaptedW(torch.nn.Module):
    """Symmetry-adapted W from a fixed equivariant basis and learned per-atom channels.

        W[m,n] = sum_c sum_A B[c, m, n, A] * psi[A, c]

    ``B`` (n_channels, K, K, n_atoms) is the equivariant basis (see
    ``src/symmetry/equivariant_w.py``); ``psi`` are per-atom learned scalars (a charge
    head with ``n_channels`` outputs). W is symmetric and equivariant by construction
    and forces the symmetry-required degeneracies at symmetric geometries.

    hippy-nn delivers ``charges`` flat as ``(n_real_atoms, n_channels)`` and positions
    padded as ``(n_mol, n_atoms, 3)``; for a single-molecule-type dataset (benzene: 12
    atoms, no padding) the flat charges reshape to ``(n_mol, n_atoms, n_channels)`` in
    atom order, which matches B's atom axis.
    """

    def __init__(self, B) -> None:
        super().__init__()
        B = np.asarray(B, dtype=float)
        self.n_channels, self.K, K2, self.n_atoms = B.shape
        if K2 != self.K:
            raise ValueError(f"B must be (n_channels, K, K, n_atoms), got {B.shape}")
        self.register_buffer("B", torch.as_tensor(B, dtype=torch.get_default_dtype()))

    def forward(self, charges: Tensor, positions: Tensor) -> Tensor:
        n_mol = positions.shape[0]
        if charges.shape[0] != n_mol * self.n_atoms:
            raise ValueError(
                f"expected {n_mol * self.n_atoms} atoms ({n_mol} mol x {self.n_atoms}); "
                f"got {charges.shape[0]} -- padded/variable-size systems not supported yet"
            )
        psi = charges.reshape(n_mol, self.n_atoms, self.n_channels)
        return torch.einsum("cmna,bac->bmn", self.B.to(psi.dtype), psi)


class SymmetryAdaptedWNode(SingleNode):
    """Graph node for :class:`SymmetryAdaptedW`. Parents: (psi_charges, positions)."""

    _index_state = IdxType.Molecules

    def __init__(self, name, parents, B, **kwargs):
        charges, positions = parents
        positions.requires_grad = True
        super().__init__(name, (charges.main_output, positions), module=SymmetryAdaptedW(B), **kwargs)


def zero_charge_head(charge_node) -> int:
    """Zero a hippy-nn HCharge head so it outputs exactly 0 at init.

    The off-diagonal V head is built from a charge head whose untrained output is
    NOT zero: with random weights plus a bias, summed over all atoms, it produces a
    large molecular value (e.g. ~5 eV for benzene's 12 atoms). That corrupts W at
    init (W = D + V with V ~ 5 eV instead of ~0), so the model does not start as the
    pretrained adiabatic energy model and must first unlearn a spurious coupling.

    Zeroing every Linear weight and bias in the head makes the per-atom charge 0,
    hence V = 0 and W = diag(D) exactly at init -- the intended ``W = D + V`` start,
    with the coupling grown from zero. Returns the number of tensors zeroed.
    """
    module = charge_node.torch_module
    n = 0
    with torch.no_grad():
        for layer in getattr(module, "layers", []):
            if getattr(layer, "weight", None) is not None:
                layer.weight.zero_(); n += 1
            if getattr(layer, "bias", None) is not None:
                layer.bias.zero_(); n += 1
    return n


def init_symmetry_psi(psi_node, basis, ref_energies, weight_jitter: float = 1e-2,
                      seed: int = 0) -> int:
    """Init the symmetry-adapted W psi head: small random last-layer weights plus the
    diagonal (site-energy) channel biases, so W starts near ``diag(ref_energies)`` --
    the symmetry-adapted analog of hierarchical-D + (near-)zero-V init.

    The diagonal bias is least-squares fit so ``W_mm = sum_c bias_c * S[c,m]`` with
    ``S[c,m] = sum_A B[c,m,m,A]`` reproduces the reference energies as closely as the
    (rank-deficient) symmetric diagonal allows.

    The ``weight_jitter`` is load-bearing, not cosmetic. With the last layer fully
    zeroed, psi would be geometry-INDEPENDENT: W is then constant in R (``dW/dR = 0``)
    and -- because the symmetric diagonal is rank-deficient at uniform psi -- exactly
    degenerate (K=9 benzene: all 8 excited states collapse onto one eigenvalue). That
    point is singular for every backward path (degenerate eigh eigenvectors, 0/0 in the
    gap-masked losses), so the first optimizer step yields NaN gradients that clipping
    can't tame. A small random last-layer weight gives psi a tiny per-atom,
    geometry-dependent variation that breaks the exact degeneracy and makes ``dW/dR``
    nonzero, so training starts in a regular region and moves away. (Genuine near-D6h
    frames encountered later are only NEAR-degenerate -- handled by detach_u and the
    float64/CPU eig fallback.) Returns the number of diagonal channels initialized.
    """
    module = psi_node.torch_module
    n_diag = basis.kind.count("diag_a1g")
    K = basis.K
    B = np.asarray(basis.B)                                   # (n_channels, K, K, n_atoms)
    S = B[:n_diag][:, np.arange(K), np.arange(K), :].sum(-1)  # (n_diag, K) diagonal contributions
    # rcond truncates the (numerically tiny) singular values of S.T: at a uniform psi
    # only the totally-symmetric A1g channels move the diagonal, so S.T is rank-deficient
    # (e.g. K=9 benzene: 4 diag channels but effective rank 2). rcond=None would keep the
    # ~1e-12 singulars and return a ~1e12 bias, which amplifies the basis's ~1e-15 off-diag
    # residual into a corrupt W. A modest cutoff gives the best reachable diagonal warm-start
    # (a few-eV residual the energy loss then closes) with an O(eV)-scale bias.
    bias_diag, *_ = np.linalg.lstsq(S.T, np.asarray(ref_energies, float), rcond=1e-6)
    with torch.no_grad():
        for layer in getattr(module, "layers", []):
            if getattr(layer, "weight", None) is not None:
                layer.weight.zero_()
            if getattr(layer, "bias", None) is not None:
                layer.bias.zero_()
        last = module.layers[-1]
        if weight_jitter > 0 and getattr(last, "weight", None) is not None:
            g = torch.Generator().manual_seed(seed)
            last.weight.copy_(torch.empty_like(last.weight).normal_(0.0, weight_jitter, generator=g))
        last.bias[:n_diag] = torch.as_tensor(bias_diag, dtype=last.bias.dtype)
    return n_diag


def _offdiag_indices(K: int) -> tuple[Tensor, Tensor]:
    """Row/col indices of the upper triangle (excluding diagonal), row-major.

    For K=3 this is rows=[0,0,1], cols=[1,2,2], i.e. (0,1),(0,2),(1,2): the same
    pair ordering used by the scaled-NACR targets.
    """
    iu = torch.triu_indices(K, K, offset=1)
    return iu[0], iu[1]


def all_pair_indices(K: int) -> list[tuple[int, int]]:
    """All upper-triangle (i, j), i<j, row-major: the full coupling pair list."""
    return [(int(i), int(j)) for i in range(K) for j in range(i + 1, K)]


class PackSymmetric(torch.nn.Module):
    """Assemble a symmetric ``(batch, K, K)`` matrix from diagonal + off-diagonal.

    Inputs:
        diag: ``(batch, K)``                     the diagonal entries (D site energies)
        off:  ``(batch, n_active_pairs)``        the ACTIVE off-diagonal entries

    ``active_pairs`` lists which ``(i, j)`` (i<j) off-diagonals ``off`` provides, in
    order. Any pair not listed is fixed to structural zero (e.g. a decoupled
    ground state). With ``active_pairs=None`` all ``K(K-1)/2`` pairs are active
    (the fully-coupled default). The result is exactly symmetric, and reduces to
    ``diag(diag)`` when no pairs are active.
    """

    def __init__(self, n_states: int, active_pairs: list[tuple[int, int]] | None = None) -> None:
        super().__init__()
        if n_states < 1:
            raise ValueError(f"n_states must be >= 1, got {n_states}")
        self.n_states = n_states
        if active_pairs is None:
            active_pairs = all_pair_indices(n_states)
        for i, j in active_pairs:
            if not (0 <= i < j < n_states):
                raise ValueError(f"active pair {(i, j)} invalid for K={n_states} (need 0<=i<j<K)")
        self.active_pairs = [(int(i), int(j)) for i, j in active_pairs]
        self.n_active_pairs = len(self.active_pairs)
        rows = torch.tensor([i for i, _ in self.active_pairs], dtype=torch.long)
        cols = torch.tensor([j for _, j in self.active_pairs], dtype=torch.long)
        self.register_buffer("_rows", rows)
        self.register_buffer("_cols", cols)

    def forward(self, diag: Tensor, off: Tensor) -> Tensor:
        K = self.n_states
        if diag.shape[-1] != K:
            raise ValueError(f"diag last dim {diag.shape[-1]} != n_states {K}")
        if off.shape[-1] != self.n_active_pairs:
            raise ValueError(f"off last dim {off.shape[-1]} != n_active_pairs {self.n_active_pairs}")
        if diag.shape[:-1] != off.shape[:-1]:
            raise ValueError(f"diag/off batch shapes differ: {diag.shape[:-1]} vs {off.shape[:-1]}")

        batch = diag.shape[:-1]
        W = diag.new_zeros((*batch, K, K))
        W[..., torch.arange(K), torch.arange(K)] = diag
        if self.n_active_pairs > 0:
            W[..., self._rows, self._cols] = off
            W[..., self._cols, self._rows] = off
        return W


class Eigenvalues(torch.nn.Module):
    """Decode the full spectrum of ``W`` -> ``(batch, K)`` ascending eigenvalues.

    These are the model's predicted state energies, energy-ordered. For monitoring
    only -- the energy *training* loss uses power sums (no eigensolver). The
    decoded-energy metric sorts the references ascending too and compares
    energy-ordered spectrum to energy-ordered spectrum, so it is robust to the
    reference state labels swapping near degeneracies (S1/S2 cross at small gaps).
    """

    def __init__(self, n_states: int) -> None:
        super().__init__()
        self.n_states = n_states

    def forward(self, W: Tensor) -> Tensor:
        return _safe_eigvalsh(W)                     # (batch, K) ascending


class EigenvaluesNode(SingleNode):
    """Monitoring node: W -> ``(batch, K)`` ascending eigenvalues (model energies)."""

    _index_state = IdxType.Molecules

    def __init__(self, name: str, W_node, n_states: int, **kwargs) -> None:
        super().__init__(name, (W_node.main_output,), module=Eigenvalues(n_states), **kwargs)


class MFrobeniusSq(torch.nn.Module):
    """Per-sample ``||M||_F^2`` from the eigenframe gradient ``M`` -> ``(batch,)``.

    By orthogonal invariance ``||M||_F = ||dW/dR||_F``, so this equals the eigh-free
    derivative magnitude WITHOUT recomputing the ``dW/dR`` Jacobian: ``M`` (from
    :class:`EigenframeGrad`) already contains it. Reused by BOTH the anchor loss and
    the smoothness penalty so the expensive Jacobian is built only once per batch.
    """

    def forward(self, M: Tensor) -> Tensor:
        # M: (batch, K, K, n_R) -> sum of squares over (state, state, nuclear DOF)
        return M.reshape(M.shape[0], -1).square().sum(dim=-1)


class MFrobeniusSqNode(SingleNode):
    """Model-side node: M -> per-sample ``||M||_F^2 = ||dW/dR||_F^2`` ``(batch,)``."""

    _index_state = IdxType.Molecules

    def __init__(self, name: str, M_node, **kwargs) -> None:
        super().__init__(name, (M_node.main_output,), module=MFrobeniusSq(), **kwargs)


class PowerSums(torch.nn.Module):
    """Power sums of the spectrum, computed directly from ``W`` (no eigensolver).

    Returns ``(batch, K)`` with entries ``p_k = Tr(W^k)`` for ``k = 1 .. K``. By
    Newton's identities these determine the unordered eigenvalues, so they are the
    energy-loss target (a stable, label-free characteristic-polynomial readout).
    """

    def __init__(self, n_states: int) -> None:
        super().__init__()
        if n_states < 1:
            raise ValueError(f"n_states must be >= 1, got {n_states}")
        self.n_states = n_states

    def forward(self, W: Tensor) -> Tensor:
        if W.shape[-2:] != (self.n_states, self.n_states):
            raise ValueError(f"W last dims {tuple(W.shape[-2:])} != (K, K)=({self.n_states},{self.n_states})")
        powers = []
        Wk = torch.eye(self.n_states, dtype=W.dtype, device=W.device).expand_as(W).clone()
        for _k in range(1, self.n_states + 1):
            Wk = Wk @ W
            powers.append(torch.diagonal(Wk, dim1=-2, dim2=-1).sum(-1))
        return torch.stack(powers, dim=-1)


class EigenframeGrad(torch.nn.Module):
    """Rotate ``dW/dR`` into the eigenframe: ``M = U^T (dW/dR) U``.

    Inputs:
        W:         ``(batch, K, K)`` symmetric, must depend on ``positions``
        positions: the nuclear coordinate tensor ``W`` was built from
                   (``requires_grad`` must be enabled)

    Output:
        M: ``(batch, K, K, n_R)`` where ``n_R`` is the flattened nuclear DOF.
           ``M[..., i, i, :]`` is ``d lambda_i / dR`` (gradients) and
           ``M[..., i, j, :]`` for ``i != j`` is the scaled coupling
           ``q_ij = (E_j - E_i) d_ij``. ``M`` is symmetric in (i, j).

    Eigenvectors come from ``torch.linalg.eigh`` (symmetric solver). Their
    backward pass carries ``1/(lambda_i - lambda_j)`` factors; near degeneracies
    use float64 and/or down-weight these outputs with a gap mask (handled by the
    loss, not here). ``dW/dR`` is taken with ``create_graph=True`` so the result is
    itself differentiable for training.

    ``active_pairs`` (same convention as :class:`PackSymmetric`) lists the
    off-diagonal entries that are actually nonzero; structurally-zero entries
    (e.g. a decoupled ground state's couplings) are SKIPPED in the Jacobian loop.
    This matters: a structurally-zero ``W_mn`` is still autograd-connected (it is a
    slice of a tensor with grad history), so differentiating it costs a full
    network backward pass that returns all zeros. For the K=3 ground-state-decoupled
    band that waste is 2 of 6 backward passes (~1/3 of the dominant Jacobian cost).
    """

    def __init__(self, n_states: int, active_pairs: list[tuple[int, int]] | None = None,
                 detach_u: bool = False, jacobian_mode: str = "auto") -> None:
        super().__init__()
        if n_states < 1:
            raise ValueError(f"n_states must be >= 1, got {n_states}")
        if jacobian_mode not in {"auto", "vectorized", "loop"}:
            raise ValueError(f"jacobian_mode must be auto|vectorized|loop, got {jacobian_mode!r}")
        self.n_states = n_states
        # detach_u: build M = U^T dW U with U DETACHED, so the loss backprops only
        # through dW/dR (stable) and not through eigh's eigenvectors (which carry the
        # 1/(lambda_i - lambda_j) blow-up). Essential at K=9 where near-degeneracies
        # are pervasive; the model still fits gradients/couplings via dW, and the
        # eigenframe is determined by W (pinned by the energy losses).
        self.detach_u = detach_u
        # jacobian_mode: how dW/dR is built. "vectorized" = one batched VJP (vmap over
        # the entries); "loop" = one autograd.grad per entry; "auto" probes the
        # vectorized path on the first forward and falls back to the loop if it raises.
        # The batched VJP is faster where it works (CPU / pure-pytorch kernels) but is
        # INCOMPATIBLE with hippy-nn's custom CUDA (Triton) kernels ("Cannot access data
        # pointer of Tensor that doesn't have storage"), so the device/kernel config
        # decides -- hence "auto" rather than a hard-coded choice.
        self.jacobian_mode = jacobian_mode
        self._use_vectorized: bool | None = {"vectorized": True, "loop": False}.get(jacobian_mode)
        if active_pairs is None:
            active_pairs = all_pair_indices(n_states)
        for i, j in active_pairs:
            if not (0 <= i < j < n_states):
                raise ValueError(f"active pair {(i, j)} invalid for K={n_states}")
        # entries to differentiate: the K diagonals + the active off-diagonals
        self.entries = [(m, m) for m in range(n_states)] + [(int(i), int(j)) for i, j in active_pairs]
        self.register_buffer("_entry_rows", torch.tensor([m for m, _ in self.entries], dtype=torch.long))
        self.register_buffer("_entry_cols", torch.tensor([n for _, n in self.entries], dtype=torch.long))

    def forward(self, W: Tensor, positions: Tensor) -> Tensor:
        K = self.n_states
        if W.shape[-2:] != (K, K):
            raise ValueError(f"W last dims {tuple(W.shape[-2:])} != (K, K)")
        if not positions.requires_grad:
            raise ValueError("positions must have requires_grad=True for EigenframeGrad")

        # eigenvectors U: (batch, K, K), columns are eigenvectors, ascending eig
        _eigvals, U = _safe_eigh(W)

        n_batch = W.shape[0]
        n_R = positions.reshape(n_batch, -1).shape[1]

        dW = self._jacobian(W, positions, n_batch, n_R)

        # M = U^T dW U, applied per nuclear DOF:  M[...,a] = U^T dW[...,a] U
        # einsum over state indices, keeping batch and nuclear-DOF axes. Detaching U
        # (when requested) removes the unstable eigenvector-sensitivity gradient.
        Ueff = U.detach() if self.detach_u else U
        M = torch.einsum("bmi,bmnr,bnj->bijr", Ueff, dW, Ueff)
        return M

    def _jacobian(self, W: Tensor, positions: Tensor, n_batch: int, n_R: int) -> Tensor:
        """dW/dR for the active entries, dispatching on ``jacobian_mode`` (with the
        first-forward probe for ``auto``)."""
        # checkpoints pickled before the vectorized path lack _use_vectorized and the
        # _entry_rows/_cols buffers; for those, only the loop is available.
        if getattr(self, "_use_vectorized", None) is None and not hasattr(self, "_entry_rows"):
            self._use_vectorized = False
        if self._use_vectorized is None:                      # auto: probe once
            try:
                dW = self._jacobian_vectorized(W, positions, n_batch, n_R)
                self._use_vectorized = True
                return dW
            except Exception as exc:                          # noqa: BLE001 - any backend failure
                warnings.warn(
                    "EigenframeGrad: batched-VJP (vectorized) Jacobian unavailable on this "
                    f"device/kernel config ({type(exc).__name__}: {exc}); using the per-entry "
                    "loop. This is expected with hippy-nn's custom CUDA kernels.",
                    RuntimeWarning, stacklevel=2,
                )
                self._use_vectorized = False
        if self._use_vectorized:
            return self._jacobian_vectorized(W, positions, n_batch, n_R)
        return self._jacobian_loop(W, positions, n_batch, n_R)

    def _jacobian_loop(self, W: Tensor, positions: Tensor, n_batch: int, n_R: int) -> Tensor:
        """One ``autograd.grad`` per entry. Universally compatible; the graph is shared
        across calls via ``retain_graph``. Structural zeros are skipped (``allow_unused``)."""
        K = self.n_states
        dW = W.new_zeros((n_batch, K, K, n_R))
        ones = torch.ones(n_batch, dtype=W.dtype, device=W.device)
        for m, n in self.entries:
            gmn = torch.autograd.grad(
                W[:, m, n], positions, grad_outputs=ones,
                create_graph=True, retain_graph=True, allow_unused=True,
            )[0]
            if gmn is None:
                continue  # this entry is constant in R -> zero derivative
            gmn = gmn.reshape(n_batch, n_R)
            dW[:, m, n, :] = gmn
            if n != m:
                dW[:, n, m, :] = gmn
        return dW

    def _jacobian_vectorized(self, W: Tensor, positions: Tensor, n_batch: int, n_R: int) -> Tensor:
        """One batched vector-Jacobian product (``is_grads_batched`` vmaps the backward
        over the E entries). The per-entry cotangent sums over the batch, which is exact
        because ``position[b]`` only influences ``W[b]`` (the batch decouples)."""
        K = self.n_states
        rows, cols = self._entry_rows, self._entry_cols
        n_entries = rows.numel()
        w_entries = W[:, rows, cols]                                  # (n_batch, E)
        cot = torch.eye(n_entries, dtype=W.dtype, device=W.device)
        cot = cot[:, None, :].expand(n_entries, n_batch, n_entries)   # (E, n_batch, E)
        jac = torch.autograd.grad(
            w_entries, positions, grad_outputs=cot, is_grads_batched=True,
            create_graph=True, retain_graph=True,
        )[0]                                                          # (E, n_batch, n_atoms, 3)
        jac = jac.reshape(n_entries, n_batch, n_R).permute(1, 0, 2)   # (n_batch, E, n_R)
        dW = W.new_zeros((n_batch, K, K, n_R))
        dW[:, rows, cols, :] = jac
        off = rows != cols
        dW[:, cols[off], rows[off], :] = jac[:, off, :]              # symmetrize off-diagonals
        return dW


def _smooth_gate(gap: Tensor, eps: float) -> Tensor:
    """Smooth shoulder ``g^2/(g^2+eps^2)``: ~1 for gap>>eps, ~0 for gap<<eps."""
    e2 = eps * eps
    return (gap * gap) / (gap * gap + e2)


class EigenGapMask(torch.nn.Module):
    """Smooth, DETACHED gate(s) on the *predicted* eigenvalue gaps of ``W``.

    The eigenvector gradient of ``torch.linalg.eigh`` carries ``1/(lambda_i -
    lambda_j)`` factors that blow up at near-degeneracies, and physically those are
    the special points where the adiabatic gauge is left undetermined (see
    docs/nqdh_model.md sections 5b/7b). This module produces detached per-sample
    gate weights ``g^2/(g^2+eps^2)`` to down-weight the eigenframe losses there.

    ``mode`` selects the granularity:
      - ``"min"``       -> ``(batch,)``: gate on the smallest adjacent gap (global).
      - ``"per_state"`` -> ``(batch, K)``: gate state ``i`` on its nearest-neighbour
                           eigenvalue gap (for the per-state gradient losses).
      - ``"per_pair"``  -> ``(batch, K(K-1)/2)``: gate pair ``(i,j)`` on
                           ``|lambda_i - lambda_j|`` in upper-triangle order (for the
                           per-pair coupling loss).

    Gating on the *predicted* gap is intentional: the spectrum is pinned everywhere
    by the eigh-free losses (power sums + the ||dW||_F^2 anchor), so the predicted
    gap converges to the true gap while remaining the quantity that actually governs
    where the frame is singular (and the quantity NAMD sees at runtime).
    """

    def __init__(self, eps: float = 0.05, mode: str = "min") -> None:
        super().__init__()
        if eps <= 0:
            raise ValueError(f"eps must be positive, got {eps}")
        if mode not in {"min", "per_state", "per_pair"}:
            raise ValueError(f"Unsupported mask mode {mode!r}")
        self.eps = float(eps)
        self.mode = mode

    def forward(self, W: Tensor) -> Tensor:
        K = W.shape[-1]
        with torch.no_grad():
            eig = _safe_eigvalsh(W)                              # (batch, K) ascending
            if K < 2:
                if self.mode == "per_state":
                    return torch.ones(W.shape[0], 1, dtype=W.dtype, device=W.device)
                return torch.ones(W.shape[0], dtype=W.dtype, device=W.device)

            if self.mode == "min":
                gaps = eig[:, 1:] - eig[:, :-1]
                return _smooth_gate(gaps.min(dim=-1).values, self.eps)

            if self.mode == "per_state":
                # nearest-neighbour gap for each state (ascending eig -> neighbours
                # are i-1 and i+1); endpoints use their single neighbour.
                left = torch.full_like(eig, float("inf"))
                right = torch.full_like(eig, float("inf"))
                left[:, 1:] = eig[:, 1:] - eig[:, :-1]
                right[:, :-1] = eig[:, 1:] - eig[:, :-1]
                nn = torch.minimum(left, right)                  # (batch, K)
                return _smooth_gate(nn, self.eps)

            # per_pair: |lambda_i - lambda_j| for upper-triangle (i<j), row-major
            rows, cols = _offdiag_indices(K)
            pair_gap = (eig[:, cols] - eig[:, rows]).abs()       # (batch, n_pairs)
            return _smooth_gate(pair_gap, self.eps)


class EigenGapMaskNode(SingleNode):
    """Model-side node: W -> eigen-gap gate(s). See :class:`EigenGapMask`.

    ``mode`` is one of ``"min"`` (batch,), ``"per_state"`` (batch, K), or
    ``"per_pair"`` (batch, n_pairs).
    """

    _index_state = IdxType.Molecules

    def __init__(self, name: str, W_node, eps: float = 0.05, mode: str = "min", **kwargs) -> None:
        super().__init__(name, (W_node.main_output,), module=EigenGapMask(eps, mode=mode), **kwargs)


class SpecialPointFlag(torch.nn.Module):
    """Mark special points: where the minimum predicted gap is below a threshold.

    Returns ``(batch,)`` boolean-as-float (1.0 = special). This is the crisp
    ``is_special`` marker exposed to NAMD (separate from the smooth training gate);
    the threshold ``delta_special`` is a model/physics property, not a numerical
    knob. See docs/nqdh_model.md section 7b for the runtime contract.
    """

    def __init__(self, delta_special: float = 0.05) -> None:
        super().__init__()
        if delta_special <= 0:
            raise ValueError(f"delta_special must be positive, got {delta_special}")
        self.delta_special = float(delta_special)

    def forward(self, W: Tensor) -> Tensor:
        with torch.no_grad():
            eig = _safe_eigvalsh(W)
            if eig.shape[-1] < 2:
                return torch.zeros(W.shape[0], dtype=W.dtype, device=W.device)
            min_gap = (eig[:, 1:] - eig[:, :-1]).min(dim=-1).values
            return (min_gap < self.delta_special).to(dtype=W.dtype)


class WFrobeniusDerivative(torch.nn.Module):
    """Eigh-free anchor: per-sample ``||dW/dR||_F^2`` (the model side).

    Equals ``sum_i ||dE_i/dR||^2 + 2 sum_{i<j} ||q_ij||^2`` by orthogonal
    invariance, so it is supervised everywhere from existing data WITHOUT an
    eigendecomposition and is smooth through conical intersections. ``positions``
    must have ``requires_grad`` (it holds the dW/dR autograd, like WRoughness).
    Returns ``(batch,)``.
    """

    def __init__(self) -> None:
        super().__init__()

    def forward(self, W: Tensor, positions: Tensor) -> Tensor:
        if not positions.requires_grad:
            raise ValueError("positions must have requires_grad=True for WFrobeniusDerivative")
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
                    continue
                sq = g.reshape(n_batch, -1).square().sum(dim=-1)
                per_sample = per_sample + (sq if m == n else 2.0 * sq)
        return per_sample


class SpecialPointFlagNode(SingleNode):
    """Inference/diagnostic node: W -> is_special ``(batch,)`` (1.0 = special)."""

    _index_state = IdxType.Molecules

    def __init__(self, name: str, W_node, delta_special: float = 0.05, **kwargs) -> None:
        super().__init__(name, (W_node.main_output,), module=SpecialPointFlag(delta_special), **kwargs)


class CoupledBlockCharacter(torch.nn.Module):
    """Analytic diabatic characters of a 2-state coupled block of ``W``.

    For block positions (i, j): delta = (W_jj - W_ii)/2, v = W_ij,
    r = sqrt(delta^2 + v^2 + eps^2). Returns ``(batch, 2, 2)``:
    rows = (lower, upper) block eigenstate (ascending), cols = (d_i, d_j):

        P[lower, d_i] = 1/2 (1 + delta/r)   (v=0, delta>0 -> lower IS d_i -> 1)

    Closed form in W's entries -- differentiable with NO eigensolver, so it teaches
    the eigenframe mixing angle even with detach-U on, and stays finite at exact
    degeneracy (eps caps the 1/r gradient; the frame genuinely rotates fast there).
    Squared characters carry no sign/phase gauge (Z2 quotiented structurally).
    """

    def __init__(self, i: int, j: int, eps: float = 1e-3) -> None:
        super().__init__()
        self.i, self.j = int(i), int(j)
        self.eps = float(eps)

    def forward(self, W: Tensor) -> Tensor:
        a = W[..., self.i, self.i]
        b = W[..., self.j, self.j]
        v = W[..., self.i, self.j]
        delta = 0.5 * (b - a)
        c2 = delta / torch.sqrt(delta ** 2 + v ** 2 + self.eps ** 2)
        lo = torch.stack([0.5 * (1 + c2), 0.5 * (1 - c2)], dim=-1)
        hi = torch.stack([0.5 * (1 - c2), 0.5 * (1 + c2)], dim=-1)
        return torch.stack([lo, hi], dim=-2)                        # (batch, 2 states, 2 diabats)


class CoupledBlockCharacterNode(SingleNode):
    """Graph node wrapping :class:`CoupledBlockCharacter`. Carries ``db_name`` so the
    loss binds ``.true`` to the stored reference characters (e.g. ``DiabChar``)."""

    _index_state = IdxType.Molecules

    def __init__(self, name: str, w_node, i: int, j: int, eps: float = 1e-3,
                 db_name: str | None = None) -> None:
        super().__init__(name, (w_node.main_output,),
                         module=CoupledBlockCharacter(i, j, eps=eps), db_name=db_name)


class WEntry(torch.nn.Module):
    """Read one entry of ``W`` -> ``(batch, 1)``; optionally the smoothed absolute
    value ``sqrt(x^2 + eps^2)`` (for sign-free targets like |J|)."""

    def __init__(self, i: int, j: int, absolute: bool = False, eps: float = 1e-4) -> None:
        super().__init__()
        self.i, self.j = int(i), int(j)
        self.absolute = bool(absolute)
        self.eps2 = float(eps) ** 2

    def forward(self, W: Tensor) -> Tensor:
        x = W[..., self.i, self.j]
        if self.absolute:
            x = torch.sqrt(x * x + self.eps2)
        return x.unsqueeze(-1)


class WEntryNode(SingleNode):
    """Graph node: one (optionally |.|-smoothed) W entry with a ``db_name`` target --
    direct supervision of the diabatic Hamiltonian's entries (e.g. TDM-diabatized
    site energies DiabE_A/DiabE_B and coupling magnitude AbsJ). Unsorted and
    eigensolver-free: the seam zero-crossing of the site-energy splitting is a
    smooth signed regression instead of a folded adiabatic gap."""

    _index_state = IdxType.Molecules

    def __init__(self, name: str, w_node, i: int, j: int, absolute: bool = False,
                 db_name: str | None = None) -> None:
        super().__init__(name, (w_node.main_output,),
                         module=WEntry(i, j, absolute=absolute), db_name=db_name)


class WBlockDelta(torch.nn.Module):
    """Signed half-splitting of a 2-state block: (W_jj - W_ii)/2 -> ``(batch, 1)``.

    Shift-invariant (the per-frame common energy error cancels), so it can be
    supervised to the ~tens-of-meV class where absolute site energies cannot --
    and its SIGN carries which diabat is lower, the seam's zero-crossing field.
    """

    def __init__(self, i: int, j: int) -> None:
        super().__init__()
        self.i, self.j = int(i), int(j)

    def forward(self, W: Tensor) -> Tensor:
        return (0.5 * (W[..., self.j, self.j] - W[..., self.i, self.i])).unsqueeze(-1)


class WBlockDeltaNode(SingleNode):
    """Graph node wrapping :class:`WBlockDelta` with a ``db_name`` target (DiabDelta)."""

    _index_state = IdxType.Molecules

    def __init__(self, name: str, w_node, i: int, j: int, db_name: str | None = None) -> None:
        super().__init__(name, (w_node.main_output,), module=WBlockDelta(i, j), db_name=db_name)


class WFrobeniusDerivativeNode(SingleNode):
    """Model-side node: (W, positions) -> per-sample ||dW/dR||_F^2 ``(batch,)``."""

    _index_state = IdxType.Molecules

    def __init__(self, name: str, W_node, positions, **kwargs) -> None:
        positions.requires_grad = True
        super().__init__(name, (W_node.main_output, positions), module=WFrobeniusDerivative(), **kwargs)


class MDiagonal(torch.nn.Module):
    """Extract the diagonal of ``M`` -> per-state gradients ``(batch, K, n_R)``.

    ``M[..., i, i, :] = d lambda_i / dR`` are the adiabatic energy gradients in the
    energy-ordered (ascending eigenvalue) convention.
    """

    def __init__(self, n_states: int) -> None:
        super().__init__()
        self.n_states = n_states

    def forward(self, M: Tensor) -> Tensor:
        # M: (batch, K, K, n_R) -> diagonal over the two state axes -> (batch, n_R, K) -> (batch, K, n_R)
        return torch.diagonal(M, dim1=1, dim2=2).transpose(-1, -2)


class MDiagonalStateNode(SingleNode):
    """Graph node: M -> one state's gradient, shaped for the gradient loss.

    Returns ``(batch, n_atoms, 3)`` for ``state_index`` (ascending-eigenvalue
    order), matching the stored ``G{s}`` target layout. ``M_ii = d lambda_i / dR``
    is the energy gradient; ``sign`` matches the dataset convention (``+1`` for
    NEXMD ``dE/dR`` targets, as used by the combined trainer's GradientNode).
    Carry a ``db_name`` so the loss binds ``.true`` to the stored gradient.
    """

    _index_state = IdxType.Molecules

    def __init__(
        self, name: str, M_node, state_index: int, n_states: int,
        n_atoms: int, sign: float = 1.0, db_name: str | None = None,
    ) -> None:
        module = _MColumn(state_index, state_index, n_states, n_atoms=n_atoms, sign=sign)
        super().__init__(name, (M_node.main_output,), module=module, db_name=db_name)


class MOffDiagonalPairs(torch.nn.Module):
    """Extract the off-diagonal pairs of ``M`` -> scaled couplings.

    Returns ``(batch, n_pairs, n_R)`` for the upper-triangle pairs (row-major,
    matching the scaled-NACR target ordering): entry ``p`` is ``M[..., i, j, :]``
    for the p-th ``(i, j)`` with ``i < j``. This is the predicted scaled NACR
    ``q_ij = (E_j - E_i) d_ij``, ready for the existing phaseless coupling loss.
    """

    def __init__(self, n_states: int) -> None:
        super().__init__()
        rows, cols = _offdiag_indices(n_states)
        self.register_buffer("_rows", rows)
        self.register_buffer("_cols", cols)

    def forward(self, M: Tensor) -> Tensor:
        # M: (batch, K, K, n_R) -> (batch, n_pairs, n_R)
        return M[:, self._rows, self._cols, :]


class _MColumn(torch.nn.Module):
    """Extract a single ``M[:, i, j, :]`` entry.

    If ``n_atoms`` is given, reshapes the flattened nuclear axis ``n_R = 3*n_atoms``
    to ``(n_atoms, 3)`` (gradient layout) and applies ``sign``; otherwise returns
    the flat ``(batch, n_R)``.
    """

    def __init__(self, i: int, j: int, n_states: int, n_atoms: int | None = None, sign: float = 1.0) -> None:
        super().__init__()
        self.i, self.j = i, j
        self.n_atoms = n_atoms
        self.sign = float(sign)

    def forward(self, M: Tensor) -> Tensor:
        col = M[:, self.i, self.j, :]
        if self.n_atoms is not None:
            col = self.sign * col.reshape(col.shape[0], self.n_atoms, 3)
        return col


class MOffDiagonalPairsNode(SingleNode):
    """Graph node wrapping :class:`MOffDiagonalPairs` (M -> scaled-NACR pairs).

    Carries ``db_name`` so the phaseless coupling loss binds ``.true`` to the
    stored ``ScaledNACR`` target. Output index state is Molecules.
    """

    _index_state = IdxType.Molecules

    def __init__(self, name: str, M_node, n_states: int, db_name: str | None = None) -> None:
        super().__init__(name, (M_node.main_output,), module=MOffDiagonalPairs(n_states), db_name=db_name)


# --------------------------------------------------------------------------- #
# hippy-nn graph-node wrappers
# --------------------------------------------------------------------------- #


class PackSymmetricNode(SingleNode):
    """Graph node wrapping :class:`PackSymmetric`.

    Parents: ``(diag_node, offdiag_node)`` where ``diag_node`` carries the K
    diagonal site energies stacked to ``(batch, K)`` and ``offdiag_node`` carries
    the ``(batch, K(K-1)/2)`` off-diagonal entries. Output index state is
    Molecules. Pass ``n_states``.
    """

    _index_state = IdxType.Molecules

    def __init__(self, name: str, parents, n_states: int, active_pairs=None, **kwargs) -> None:
        diag_node, offdiag_node = parents
        parents = (diag_node.main_output, offdiag_node.main_output)
        super().__init__(name, parents, module=PackSymmetric(n_states, active_pairs=active_pairs), **kwargs)


class PowerSumsNode(SingleNode):
    """Graph node wrapping :class:`PowerSums` (W -> Tr(W^k))."""

    _index_state = IdxType.Molecules

    def __init__(self, name: str, parents, n_states: int, **kwargs) -> None:
        (w_node,) = parents
        super().__init__(name, (w_node.main_output,), module=PowerSums(n_states), **kwargs)


class EigenframeGradNode(SingleNode):
    """Graph node wrapping :class:`EigenframeGrad` (W, positions -> M).

    Parents: ``(w_node, positions_node)``. Sets ``requires_grad`` on positions, as
    the NACR nodes do, so the internal ``dW/dR`` autograd works. Pass the same
    ``active_pairs`` as the PackSymmetric node so structurally-zero entries are
    skipped in the Jacobian loop.
    """

    _index_state = IdxType.Molecules

    def __init__(self, name: str, parents, n_states: int, active_pairs=None,
                 detach_u: bool = False, jacobian_mode: str = "auto", **kwargs) -> None:
        w_node, positions = parents
        positions.requires_grad = True
        super().__init__(name, (w_node.main_output, positions),
                         module=EigenframeGrad(n_states, active_pairs=active_pairs,
                                               detach_u=detach_u, jacobian_mode=jacobian_mode),
                         **kwargs)
