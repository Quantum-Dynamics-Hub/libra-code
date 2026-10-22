"""QTAG specialization of the common dynamics engine.

Numerical kernels live in dyn.gwp, propagation, initialization and observables.
The former standalone implementation is retained only as a test oracle.
"""
from .method import QTAGMethod


def QTAGEngine(manager, model, params=None):
    """Compatibility constructor returning DynamicsEngine(method='qtag')."""
    from ...engine import DynamicsEngine
    return DynamicsEngine(manager=manager, model_fn=model, method='qtag', params=params)


def run_qtag(manager, model, params, savers=None, keep_history=True):
    """Convenience wrapper around the common engine's run loop."""
    return QTAGEngine(manager, model, params).run(savers=savers, keep_history=keep_history)


__all__ = ['QTAGMethod', 'QTAGEngine', 'run_qtag']
