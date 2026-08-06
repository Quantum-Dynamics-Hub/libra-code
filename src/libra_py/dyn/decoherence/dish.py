"""DISH operations translated from ``src/dyn/dyn_methods_dish.cpp``."""
from __future__ import annotations
import numpy as np
from .methods import collapse, project_out
from .times import coherence_intervals

def decoherence_event(coherence_time, coherence_interval, option=0, rng=None):
    """Select at most one decohering state per trajectory.

    ``coherence_time`` and ``coherence_interval`` are trajectory-major arrays
    shaped ``(ntraj,nstates)``. Option 0 compares clocks directly with the
    interval; option 1 samples an exponential waiting time whose mean is the
    interval, as intended by the original DISH paper. If several states are
    eligible, one is sampled uniformly.

    Returns an integer array of length ``ntraj``. ``-1`` denotes no event.
    """
    rng=rng or np.random.default_rng(); clocks=np.asarray(coherence_time); tau=np.asarray(coherence_interval)
    out=np.full(clocks.shape[0],-1,dtype=int)
    for t in range(len(clocks)):
        threshold=tau[t] if option==0 else rng.exponential(tau[t])
        possible=np.flatnonzero(clocks[t] >= threshold)
        if len(possible): out[t]=rng.choice(possible)
    return out

def dish_hop_proposal(active_states, amplitudes, coherence_time, rates, rng=None, event_option=0):
    """Propose old-style DISH hops and handle unselected events.

    Coherence intervals are computed from ``amplitudes`` and pair ``rates``.
    At an event, the affected state is proposed with probability ``|c_i|^2``;
    otherwise it is projected out. Event clocks and amplitudes are modified in
    place. The returned proposed states default to the current active states.
    """
    rng=rng or np.random.default_rng(); interval=coherence_intervals(amplitudes,rates)
    events=decoherence_event(coherence_time,interval,event_option,rng); proposed=np.array(active_states,copy=True)
    for t,state in enumerate(events):
        if state >= 0:
            coherence_time[t,state]=0.0
            if rng.random() <= abs(amplitudes[t,state])**2: proposed[t]=state
            else: project_out(amplitudes,state,t)
    return proposed

def dish_project_out_collapse(old_states, proposed_states, new_states, amplitudes, collapse_option=0):
    """Finalize old-style DISH after energetic hop acceptance.

    An accepted nontrivial proposal collapses onto the new state; a rejected
    proposal is projected out. A trivial proposal collapses onto the original
    state, following C++ ``dish_project_out_collapse``. Amplitudes are modified
    in place and returned.
    """
    for t,(old,prop,new) in enumerate(zip(old_states,proposed_states,new_states)):
        if prop != old:
            collapse(amplitudes,new,t,collapse_option) if new==prop else project_out(amplitudes,prop,t)
        else: collapse(amplitudes,old,t,collapse_option)
    return amplitudes

def dish_rev2023(amplitudes, active_states, coherence_time, rates, dt,
                 times_type=-1, event_option=1, collapse_option=0, rng=None):
    """Apply revised-2023 DISH collapse/project events without hopping.

    All state clocks are advanced by ``dt``. For Schwartz-1 time types 2 and
    4, diagonal inverse times are used directly; other time models use
    :func:`coherence_intervals`. At an event, state ``i`` is collapsed with
    probability ``|c_i|^2`` and projected out otherwise.

    ``amplitudes`` and ``coherence_time`` are modified in place. The returned
    integer event array uses ``-1`` for trajectories without an event.
    ``active_states`` is retained in the signature for correspondence with the
    C++ workflow, although revised DISH does not use it in the event decision.
    """
    rng=rng or np.random.default_rng(); coherence_time += dt
    if times_type in (2,4):
        diag=np.diagonal(rates,axis1=1,axis2=2); intervals=np.full_like(diag,1e10)
        np.divide(1.0,diag,out=intervals,where=np.abs(diag)>0)
    else: intervals=coherence_intervals(amplitudes,rates)
    events=decoherence_event(coherence_time,intervals,event_option,rng)
    for t,state in enumerate(events):
        if state >= 0:
            coherence_time[t,state]=0.0
            if rng.random() <= abs(amplitudes[t,state])**2: collapse(amplitudes,state,t,collapse_option)
            else: project_out(amplitudes,state,t)
    return events
