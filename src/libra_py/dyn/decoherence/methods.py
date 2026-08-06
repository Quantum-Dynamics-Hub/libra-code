"""Wavefunction decoherence operations from ``dyn_decoherence_methods.cpp``."""

from __future__ import annotations
import numpy as np


def project_out(amplitudes, state, trajectory=None):
    """Project one state out and renormalize the remaining superposition.

    ``amplitudes`` is either ``(nstates,)`` or trajectory-major
    ``(ntraj,nstates)``. For a batch, ``trajectory`` is required. State
    ``state`` is set to zero and the remaining coefficients are divided by
    ``sqrt(1-|c_state|**2)``. The array is modified in place and returned.

    If the projected state carries all population, the whole vector becomes
    zero. This matches ``project_out`` in ``dyn_decoherence_methods.cpp``.
    """
    c = amplitudes if isinstance(amplitudes, np.ndarray) else np.asarray(amplitudes)
    row = c if c.ndim == 1 else c[int(trajectory)]
    norm = max(0.0, 1.0 - float(abs(row[state]) ** 2))
    row *= 1.0 / np.sqrt(norm) if norm > 0.0 else 0.0
    row[state] = 0.0j
    return amplitudes


def collapse(amplitudes, state, trajectory=None, collapse_option=0):
    """Collapse a coefficient vector onto one electronic state in place.

    ``collapse_option=0`` preserves ``c_state/|c_state|`` when nonzero;
    ``collapse_option=1`` assigns ``1+0j`` and therefore discards phase.
    Batched amplitudes use shape ``(ntraj,nstates)`` and require the trajectory
    index. The modified input is returned for convenient composition.
    """
    c = amplitudes if isinstance(amplitudes, np.ndarray) else np.asarray(amplitudes)
    row = c if c.ndim == 1 else c[int(trajectory)]
    old = row[state]; row[:] = 0.0j
    row[state] = old / abs(old) if collapse_option == 0 and abs(old) > 0.0 else 1.0 + 0.0j
    return amplitudes


def collapse_dm(density_matrix, state):
    """Collapse a square density matrix onto ``|state><state|`` in place.

    All elements are cleared and the selected diagonal element is set to one.
    The modified matrix is returned.
    """
    density_matrix[...] = 0.0j; density_matrix[state, state] = 1.0 + 0.0j
    return density_matrix


def sdm(amplitudes, dt, active_states, rates, tolerance=0.0):
    """Apply Simplified Decay of Mixing (SDM/mSDM).

    Inactive amplitude ``i`` decays by ``exp(-dt*rate[i,active])`` and the
    active amplitude is rescaled to conserve unit norm and preserve phase.

    Parameters use trajectory-major shapes ``(ntraj,nstates)`` and
    ``(ntraj,nstates,nstates)``; a single rate matrix is broadcast across
    trajectories. ``dt`` is in atomic units. ``tolerance`` has the same role
    as C++ ``sdm_norm_tolerance`` when active population slightly exceeds one.

    A new complex array is returned; the input coefficients are not modified.
    """
    c = np.array(amplitudes, dtype=complex, copy=True)
    single = c.ndim == 1; c = c[None] if single else c
    states = np.broadcast_to(np.asarray(active_states, dtype=int), (len(c),))
    r = np.asarray(rates, dtype=float); r = np.broadcast_to(r if r.ndim==3 else r[None], (len(c),)+r.shape[-2:])
    for t, active in enumerate(states):
        paa = abs(c[t, active]) ** 2
        if paa > 1.0 + tolerance:
            c[t] /= np.sqrt(paa); paa = 1.0
        for i in range(c.shape[1]):
            if i != active: c[t, i] *= np.exp(-dt * r[t, i, active])
        inactive = np.sum(np.abs(c[t]) ** 2) - abs(c[t, active]) ** 2
        if inactive > 1.0 + 1e-12:
            raise ValueError("SDM inactive-state population exceeds one")
        if paa > 0.0: c[t, active] *= np.sqrt(max(0.0, 1.0-inactive) / paa)
    return c[0] if single else c


def instantaneous_decoherence(amplitudes, accepted_states, proposed_states,
                              initial_states, variant=1, collapse_option=0):
    """Apply an instantaneous-decoherence correction in place.

    The state arrays contain one entry per trajectory and correspond to the
    states before proposal, after proposal, and after acceptance. Variant
    numbering is identical to C++:

    ``0`` ID-S, collapse only after successful hops;
    ``1`` ID-A, collapse after every nontrivial attempt;
    ``2`` ID-C, collapse onto every accepted state;
    ``3`` IDN, project out a rejected proposed state;
    ``4`` ID-F, project only after frustrated hops.

    ``collapse_option`` is passed to :func:`collapse`. The modified amplitude
    batch is returned.
    """
    c = amplitudes
    accepted=np.asarray(accepted_states); proposed=np.asarray(proposed_states); initial=np.asarray(initial_states)
    for t in range(len(accepted)):
        attempted = proposed[t] != initial[t]
        success = accepted[t] == proposed[t]
        if variant == 0 and accepted[t] != initial[t]: collapse(c, accepted[t], t, collapse_option)
        elif variant == 1 and attempted: collapse(c, accepted[t] if success else initial[t], t, collapse_option)
        elif variant == 2: collapse(c, accepted[t], t, collapse_option)
        elif variant == 3 and attempted: (collapse(c, accepted[t], t, collapse_option) if success else project_out(c, proposed[t], t))
        elif variant == 4 and attempted and not success: project_out(c, proposed[t], t)
    return c


def bcsh(amplitudes, active_states, reversal_events):
    """Apply branching-corrected surface-hopping coefficient resets.

    ``reversal_events[t,i] > 0`` identifies reflection of the wavepacket on
    state ``i`` of trajectory ``t``. Reflected inactive branches are projected
    out; an active branch is retained. A corrected copy of amplitudes is
    returned. Event detection itself belongs to the nuclear BCSH workflow.
    """
    c=np.array(amplitudes, dtype=complex, copy=True); rev=np.asarray(reversal_events)
    for t,a in enumerate(np.asarray(active_states, dtype=int)):
        for i in range(c.shape[1]):
            if i != a and rev[t,i] > 0: project_out(c, i, t)
    return c


def afssh_dzdt(dz, hvib, force_matrix, amplitudes, mass, active_state):
    """Evaluate the AFSSH first-moment differential equations.

    ``dz`` is a ``(2*nstates,2*nstates)`` block matrix containing ``dR`` in
    its upper-left and ``dP`` in its lower-right block. ``hvib`` and the
    diagonal state-force matrix are ``(nstates,nstates)``; ``amplitudes`` is a
    state vector. ``mass`` and ``active_state`` refer to one DOF and one
    trajectory. The returned derivative is a new block matrix.

    This is equations 15--16 of J. Chem. Theory Comput. 2016, 12, 5256 and a
    direct analogue of the C++ ``afssh_dzdt`` overload.
    """
    dz=np.asarray(dz); n=dz.shape[0]//2; dR=dz[:n,:n]; dP=dz[n:,n:]
    h=np.asarray(hvib); f=np.asarray(force_matrix); c=np.asarray(amplitudes).reshape(n,1)
    sigma=c@c.conj().T; dF=f-np.eye(n)*f[active_state,active_state]
    out=np.zeros_like(dz); out[:n,:n]=-1j*(h@dR-dR@h)-dP/mass
    out[n:,n:]=-1j*(h@dP-dP@h)-0.5*(dF@sigma+sigma@dF)
    return out


def integrate_afssh_moments(dR, dP, hvib, force_matrix, amplitudes, mass,
                            active_state, dt, nsteps=1):
    """Integrate AFSSH position and momentum moments using classical RK4.

    This single-DOF, single-trajectory operation advances square ``dR`` and
    ``dP`` matrices through ``nsteps`` substeps of size ``dt``. Both matrices
    are modified in place and returned as ``(dR, dP)``. Other arguments are
    passed unchanged to :func:`afssh_dzdt`.
    """
    n=len(dR); z=np.zeros((2*n,2*n),complex); z[:n,:n]=dR; z[n:,n:]=dP
    rhs=lambda x: afssh_dzdt(x,hvib,force_matrix,amplitudes,mass,active_state)
    for _ in range(nsteps):
        k1=rhs(z); k2=rhs(z+.5*dt*k1); k3=rhs(z+.5*dt*k2); k4=rhs(z+dt*k3)
        z += dt*(k1+2*k2+2*k3+k4)/6
    dR[...] = z[:n,:n]; dP[...] = z[n:,n:]
    return dR, dP
