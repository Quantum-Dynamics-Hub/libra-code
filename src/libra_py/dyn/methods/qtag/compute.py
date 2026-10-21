# Copyright (C) 2022-2026 Matthew Dutra, Alexey V. Akimov
# SPDX-License-Identifier: GPL-3.0-or-later

"""QTAG moving-basis propagation using only NumPy or PyTorch arrays."""
from dataclasses import dataclass
import numpy as np
from ...gwp import integrals as gwp
from . import core
from ._backend import xp_of, asarray, zeros, stack, copy, indices, scalar


@dataclass
class QTAGResult:
    """Final state plus history; iteration yields the legacy six final values.

    history includes t=0 and every completed step (nsteps+1 snapshots).
    Tensor histories retain graphs; use torch.no_grad() for production runs
    that do not require parameter derivatives. File/plot helpers detach only
    their copies, without changing this object.
    """
    q: object
    p: object
    alp: object
    s: object
    states: list
    coeff: object
    history: dict

    def __iter__(self):
        return iter((self.q, self.p, self.alp, self.s, self.states, self.coeff))


def time_overlap(nQ, nP, nA, nS, nstate, oQ, oP, oA, oS, ostate):
    """Return <G_new|G_old>, including state masks and QTAG half-widths."""
    return gwp.gwp_overlap_matrix(nQ, nP, nS, nA/2, nstate, oQ, oP, oS, oA/2, ostate)


def qtag_momentum(q, p, a, s, coeff_on_surf, dyn_params=None):
    """Return momentum, real log-gradient, and their slopes (ndof,ntraj).

    mom_calc_type=0 returns raw ratios (zero slopes, matching legacy updates);
    1 uses a density-weighted linear least-squares fit per coordinate;
    2 returns raw ratios with their analytic spatial derivatives.
    d_weight=0 disables density weighting. fit_rcond controls a spectral
    pseudoinverse of the 2x2 normal matrix. Unlike the iterative legacy solver,
    linfit_beta is not needed; a direct solve gives deterministic accuracy.
    """
    pars = dyn_params or {}
    q, p, a, s = core._basis(q, p, a, s)
    ratio, slope, density = core.momentum_fields(q, p, a, s, coeff_on_surf,
                                                pars.get("node_tol", 1e-14))
    kind = pars.get("mom_calc_type", 1)
    xp = xp_of(q)
    if kind == 0:
        return ratio.imag, ratio.real, xp.zeros_like(q), xp.zeros_like(q)
    if kind == 2:
        return ratio.imag, ratio.real, slope.imag, slope.real
    if kind != 1:
        raise ValueError("mom_calc_type must be 0, 1 or 2")
    weight = density if pars.get("d_weight", 1) else xp.ones_like(density)
    if scalar(weight.max()) == 0:
        return (xp.zeros_like(q),)*4
    weight = weight/weight.max()
    fitted, slopes = [], []
    for d in range(q.shape[0]):
        design = stack([xp.ones_like(q[d]), q[d]], axis=1)
        normal = design.T @ (weight[:, None]*design)
        rhs = (design.T+0j) @ (weight*ratio[d])
        # normal is real but rhs is complex; promote before PyTorch matmul.
        fit = core.metric_solve(normal+0j, rhs, pars.get("fit_rcond", 1e-12))
        fitted.append((design+0j) @ fit)
        slopes.append(xp.ones_like(q[d])*fit[1])
    fitted, slopes = stack(fitted), stack(slopes)
    return fitted.imag, fitted.real, slopes.imag, slopes.real


def qtag_pops(surf_ids, coeff, S, target_states):
    """Return state populations C_n^H S_nn C_n in target_states order."""
    c = asarray(coeff, S, True).reshape(-1)
    labels = indices(surf_ids, c)
    values = []
    for state in target_states:
        masked = c*(labels == state)
        values.append((masked.conj() @ S @ masked).real)
    return stack(values)


def qtag_energy(coeff, H):
    """Return the unnormalized expectation C^H H C (a real scalar)."""
    c = asarray(coeff, H, True).reshape(-1)
    return (c.conj() @ H @ c).real


def coefficient_matrix(coeff, labels, target_states):
    """Map one coefficient per labeled packet to (len(target_states),ntraj)."""
    c = coeff.reshape(-1)
    st = indices(labels, c)
    return stack([c*(st == n) for n in target_states])


def propagate_basis(q, p, alp, s, states, coeff, dyn_params, surf_pops):
    """Euler-update the Gaussian basis, then project coefficients onto it.

    q_update_method=1 uses q += dt*invM*momentum; p_update_method=1 replaces
    p with quantum momentum; a_update_method=1 uses alp -= 2*dt*invM*alp*gmom.
    Zero freezes a parameter. Phase evolution is not implemented in legacy
    QTAG: nonzero s_update_method raises instead of silently doing nothing.
    Low-population surfaces synchronize selected parameters to the largest
    population surface (or freeze them); synchronization needs equal counts.
    Stable sorting handles ties. The most populated nonzero state is always
    evolved, even if below decpl_den. A negative new width raises with advice
    to reduce dt; widths and populations are not silently clipped.
    """
    pars = dyn_params
    dt = pars["dt"]
    invM = asarray(pars["iM"], q).reshape(-1, 1)
    target = list(pars.get("states", sorted(set(states))))
    if set(target) != set(states) or len(target) != len(set(target)):
        raise ValueError("states must list each occupied basis label exactly once")
    if len(surf_pops) != len(target):
        raise ValueError("surf_pops must correspond to states")
    if pars.get("s_update_method", 0) != 0:
        raise NotImplementedError("QTAG phase evolution is not implemented; use s_update_method=0")
    updates = [pars.get(k+"_update_method", default) for k, default in (("q", 1), ("p", 1), ("a", 0))]
    sync = [pars.get(k+"_sync_method", default) for k, default in (("q", 1), ("p", 1), ("a", 0), ("s", 0))]
    if any(v not in (0, 1) for v in updates+sync):
        raise ValueError("update and sync methods must be 0 or 1")
    new = [copy(x) for x in (q, p, alp, s)]
    order = sorted(range(len(target)), key=lambda k: -scalar(surf_pops[k]))
    reference = None
    for rank, k in enumerate(order):
        state = target[k]
        idx = indices([i for i, n in enumerate(states) if n == state], q)
        old = [x[:, idx] for x in (q, p, alp, s)]
        pop = scalar(surf_pops[k])
        if pop > pars.get("decpl_den", .25) or (rank == 0 and pop > 0):
            mom, r, gmom, gr = qtag_momentum(*old, coeff.reshape(-1)[idx], pars)
            block = [old[0]+dt*invM*mom if updates[0] else old[0],
                     mom if updates[1] else old[1],
                     old[2]-2*dt*invM*old[2]*gmom if updates[2] else old[2], old[3]]
            if not bool((block[2] > 0).all()):
                raise ValueError("basis update produced nonpositive widths; reduce dt or freeze widths")
            if rank == 0:
                reference = block
        elif reference is not None:
            if any(sync) and old[0].shape != reference[0].shape:
                raise ValueError("synchronization requires equal packet counts on the two states")
            block = [reference[j] if sync[j] else old[j] for j in range(4)]
        else:
            block = old
        for dest, value in zip(new, block):
            dest[:, idx] = value
    nq, np_, na, ns = new
    Snew = time_overlap(nq, np_, na, ns, states, nq, np_, na, ns, states)
    cross = time_overlap(nq, np_, na, ns, states, q, p, alp, s, states)
    cnew = core.metric_solve(Snew, cross @ coeff.reshape(-1, 1), pars.get("overlap_rcond", 1e-12))
    return nq, np_, na, ns, list(states), cnew


def run_qtag(q, p, alp, s, states, coeff, invM, dyn_params, compute_model, model_params=None):
    """Run symmetric electronic-half-step / basis-update / electronic-half-step.

    Input q,p,alp,s: (ndof,ntraj), coeff: (ntraj,) or (ntraj,1), invM: one
    positive inverse mass per DOF. Parameters: dt>0, nsteps>=0; see
    propagate_basis and core.qtag_hamiltonian_and_overlap for update choices.
    overlap_rcond defaults to 1e-12. target_states controls saved populations,
    not the propagation state set. Outputs contain the INITIAL snapshot and
    nsteps final snapshots, with observables evaluated at the saved time.

    Returns QTAGResult (iterable as q,p,alp,s,states,coeff). Saving is explicit:
    save.save_qtag_data(result,path). No files, printouts, or input mutations
    occur here. Singular overlap directions are truncated; moving-basis
    projection can lose norm. The history reports it rather than renormalizing.
    """
    q, p, alp, s = [copy(x) for x in core._basis(q, p, alp, s)]
    states = core._states(states, q.shape[1])
    c = asarray(coeff, q, True).reshape(-1, 1)
    if c.shape[0] != q.shape[1]:
        raise ValueError("one coefficient per packet is required")
    mass = asarray(invM, q).reshape(-1)
    if mass.shape[0] != q.shape[0] or not bool((mass > 0).all()):
        raise ValueError("one positive inverse mass per DOF is required")
    pars = dict(dyn_params)
    dt, nsteps = pars.get("dt", 1.), pars.get("nsteps", 1)
    if not np.isfinite(dt) or dt <= 0 or int(nsteps) != nsteps or nsteps < 0:
        raise ValueError("dt must be positive and nsteps a nonnegative integer")
    if not bool(xp_of(q).isfinite(c).all()) or not bool(xp_of(q).isfinite(mass).all()):
        raise ValueError("coefficients and inverse masses must be finite")
    all_states = list(pars.get("states", sorted(set(states))))
    if set(all_states) != set(states) or len(all_states) != len(set(all_states)):
        raise ValueError("states must list each occupied basis label exactly once")
    targets = list(pars.get("target_states", all_states))
    if not targets or any(int(n) != n or n < 0 for n in targets):
        raise ValueError("target_states must contain nonnegative integer labels")
    pars.update(dt=dt, iM=mass, states=all_states)
    rc = pars.get("overlap_rcond", 1e-12)
    history = {k: [] for k in ("time", "Etot", "dEtot", "norm", "pops", "q", "p", "a", "s", "coeffs")}
    S, H = core.qtag_hamiltonian_and_overlap(q, p, alp, s, states, mass, compute_model, model_params, pars)
    initial_energy = qtag_energy(c, H)
    def record(step):
        energy = qtag_energy(c, H)
        vals = dict(time=asarray(step*dt, q), Etot=energy, dEtot=energy-initial_energy,
                    norm=(c.conj().T @ S @ c)[0, 0].real,
                    pops=qtag_pops(states, c, S, targets), q=q, p=p, a=alp, s=s, coeffs=c[:, 0])
        for key, value in vals.items():
            history[key].append(copy(value))
    record(0)
    for step in range(int(nsteps)):
        c = core.propagate_electronic_qtag(dt/2, c, H, S, rc)
        pops = qtag_pops(states, c, S, all_states)
        q, p, alp, s, states, c = propagate_basis(q, p, alp, s, states, c, pars, pops)
        step_params = dict(model_params or {}, timestep=step+1)
        S, H = core.qtag_hamiltonian_and_overlap(q, p, alp, s, states, mass, compute_model, step_params, pars)
        c = core.propagate_electronic_qtag(dt/2, c, H, S, rc)
        record(step+1)
    history = {key: stack(value) for key, value in history.items()}
    history.update(states=indices(states, q), target_states=indices(targets, q))
    return QTAGResult(q, p, alp, s, states, c, history)


def wf_calc_nD(result, params):
    """Reconstruct complex wavefunctions/densities on a Cartesian grid.

    result is QTAGResult or a loaded history dict. params contains xmin/xmax,
    npoints (one per DOF), snaps (saved snapshot indices), and optionally
    which_states. Returns axes, points, psi, density, snaps, states, times.
    psi/density have shape (nsnaps,nstates,*npoints). This reconstructs all
    phases and does not assume equal counts or sorted state labels.
    """
    h = result.history if isinstance(result, QTAGResult) else result
    like = h["q"]
    d = like.shape[1]
    lo, hi, counts = params.get("xmin", [-4.]*d), params.get("xmax", [4.]*d), params.get("npoints", [100]*d)
    if not len(lo) == len(hi) == len(counts) == d or any(int(n) != n or n < 2 for n in counts):
        raise ValueError("grid bounds/counts must contain one entry per DOF (at least two points)")
    if any(b <= a for a, b in zip(lo, hi)):
        raise ValueError("xmax must exceed xmin")
    xp = xp_of(like)
    axes = [asarray(np.linspace(a, b, int(n)), like) for a, b, n in zip(lo, hi, counts)]
    points = stack(xp.meshgrid(*axes, indexing="ij")).reshape(d, -1)
    snaps = list(params.get("snaps", [0, len(h["time"])-1]))
    states = list(params.get("which_states", [int(x) for x in h["target_states"]]))
    psi = []
    for k in snaps:
        if not 0 <= k < len(h["time"]):
            raise ValueError("snapshot index is outside recorded history")
        c = coefficient_matrix(h["coeffs"][k], h["states"], states)
        psi.append(core.qtag_psi(points, h["q"][k], h["p"][k], h["a"][k], h["s"][k], c)
                   .reshape(len(states), *counts))
    values = stack(psi)
    return dict(axes=axes, points=points, psi=values, density=abs(values)**2,
                snaps=snaps, states=states, times=h["time"][indices(snaps, like)])
