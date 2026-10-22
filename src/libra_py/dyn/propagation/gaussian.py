# Copyright (C) 2022-2026 Matthew Dutra, Alexey V. Akimov
# SPDX-License-Identifier: GPL-3.0-or-later

"""Adaptive Gaussian basis evolution and coefficient projection.

Kernels accept (ndof, npackets) arrays without owning trajectory state.
Quantum-momentum estimation is provided by algorithms.quantum_potential.
"""
from ..gwp import integrals as gwp
from .nuclear import drift_coordinates
from ..algorithms.quantum_potential import quantum_momentum
from ..utils.linalg import metric_solve
from ..backends.arrays import asarray, copy, indices, scalar

def time_overlap(nQ, nP, nA, nS, nstate, oQ, oP, oA, oS, ostate):
    """Return <G_new|G_old>, including state masks and QTAG half-widths."""
    return gwp.gwp_overlap_matrix(nQ, nP, nS, nA/2, nstate, oQ, oP, oS, oA/2, ostate)


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
            mom, r, gmom, gr = quantum_momentum(*old, coeff.reshape(-1)[idx], pars)
            block = [drift_coordinates(old[0], mom, invM, dt) if updates[0] else old[0],
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
    cnew = metric_solve(Snew, cross @ coeff.reshape(-1, 1), pars.get("overlap_rcond", 1e-12))
    return nq, np_, na, ns, list(states), cnew
