"""
Small generic integration helpers.
"""

from __future__ import annotations


def normalize_dt(dt):
    """Return a floating timestep and reject non-positive values."""

    dt = float(dt)
    if dt <= 0.0:
        raise ValueError("dt must be positive")
    return dt


def run_steps(step_fn, nsteps: int, dt):
    """Call `step_fn(dt)` `nsteps` times and collect returned summaries."""

    if nsteps < 0:
        raise ValueError("nsteps must be non-negative")
    dt = normalize_dt(dt)
    return [step_fn(dt) for _ in range(nsteps)]
