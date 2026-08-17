"""Tests for the gradient_state / active-root handling in
libra_py.packages.pyscf.methods (strategy_compute_adi and friends).

These exercise model_params["gradient_state"] in {None, int, "all", "active"}
and the "active" resolution against model_params["act_state"][itraj], which
libra_py.dynamics.tsh.compute writes every step.
"""

from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[4]
SRC = ROOT / "src"

if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from liblibra_core import CMATRIX, MATRIX

from libra_py.packages.pyscf.interfaces import ES_Strategy, MolecularGeometry
from libra_py.packages.pyscf.methods import (
    _params_to_request,
    _resolve_gradient_state,
    strategy_compute_adi,
)


# =============================================================================
# _resolve_gradient_state / _params_to_request: pure-python, no ES backend
# =============================================================================

def test_gradient_state_none_means_no_gradient():
    assert _resolve_gradient_state({"gradient_state": None}, itraj=0) is None


def test_gradient_state_all_means_all_roots():
    assert _resolve_gradient_state({"gradient_state": "all"}, itraj=0) == "all"


def test_gradient_state_explicit_int_passes_through():
    assert _resolve_gradient_state({"gradient_state": 1}, itraj=0) == 1


def test_gradient_state_active_resolves_per_trajectory():
    params = {"gradient_state": "active", "act_state": {0: 2, 1: 0}}
    assert _resolve_gradient_state(params, itraj=0) == 2
    assert _resolve_gradient_state(params, itraj=1) == 0


def test_gradient_state_active_without_act_state_raises():
    with pytest.raises(KeyError):
        _resolve_gradient_state({"gradient_state": "active"}, itraj=0)


def test_params_to_request_defaults_to_all_roots():
    request = _params_to_request({"nstates": 3}, itraj=0)
    assert request.gradient_state == "all"
    assert request.n_singlets == 3
    assert request.n_triplets == 0


def test_params_to_request_active_root_only():
    params = {"nstates": 3, "gradient_state": "active", "act_state": {0: 1}}
    request = _params_to_request(params, itraj=0)
    assert request.gradient_state == 1


# =============================================================================
# strategy_compute_adi: end-to-end through a fake ES_Strategy backend
# =============================================================================

class _FakeStrategy(ES_Strategy):
    """Two-state backend: energies are fixed, gradient(root) = (root + 1) * ones."""

    nstates = 2

    def __init__(self):
        self._geom = None

    def snapshot_state(self):
        pass

    def get_state(self):
        return None

    def get_previous_state(self):
        return None

    def set_geom(self, geom):
        self._geom = geom

    def get_geom(self):
        return self._geom

    def compute_H_el(self):
        # Real backends (e.g. CISD.compute_H_el) return a 1D energy vector.
        return np.array([-1.0, -0.5])

    def compute_gradient(self, root=0):
        natoms = len(self._geom.atom_labels)
        return np.full((natoms, 3), float(root + 1))

    def compute_all_gradients(self):
        return [self.compute_gradient(root) for root in range(self.nstates)]


def _make_q(coords_bohr):
    # Nuclear coordinates are real-valued (MATRIX), matching what the
    # dynamics driver passes to compute_model in practice.
    natoms = len(coords_bohr)
    q = MATRIX(3 * natoms, 1)
    for atom in range(natoms):
        for xyz in range(3):
            q.set(3 * atom + xyz, 0, float(coords_bohr[atom][xyz]))
    return q


def _base_params():
    return {
        "atom_labels": ["H", "H"],
        "es_strategy": _FakeStrategy(),
        "nstates": 2,
        "time_overlap": False,
    }


def test_strategy_compute_adi_active_root_only_gradient():
    q = _make_q([[0.0, 0.0, 0.0], [0.0, 0.0, 1.4]])
    params = _base_params()
    params["gradient_state"] = "active"
    params["act_state"] = {0: 1}  # trajectory 0's active root is state 1

    obj = strategy_compute_adi(q, params, full_id=[0, 0])

    # Only state 1's diagonal block should be nonzero; state 0's should be zero.
    for dof in range(6):
        block = obj.d1ham_adi[dof]
        assert block.get(1, 1).real == pytest.approx(2.0)
        assert block.get(0, 0).real == pytest.approx(0.0)


def test_strategy_compute_adi_all_roots_gradient():
    q = _make_q([[0.0, 0.0, 0.0], [0.0, 0.0, 1.4]])
    params = _base_params()
    params["gradient_state"] = "all"

    obj = strategy_compute_adi(q, params, full_id=[0, 0])

    for dof in range(6):
        block = obj.d1ham_adi[dof]
        assert block.get(0, 0).real == pytest.approx(1.0)
        assert block.get(1, 1).real == pytest.approx(2.0)


def test_strategy_compute_adi_no_gradient():
    q = _make_q([[0.0, 0.0, 0.0], [0.0, 0.0, 1.4]])
    params = _base_params()
    params["gradient_state"] = None

    obj = strategy_compute_adi(q, params, full_id=[0, 0])

    assert not hasattr(obj, "d1ham_adi")


def test_strategy_compute_adi_active_root_differs_per_trajectory():
    q = _make_q([[0.0, 0.0, 0.0], [0.0, 0.0, 1.4]])

    strategies = {0: _FakeStrategy(), 1: _FakeStrategy()}
    act_state = {0: 0, 1: 1}

    for itraj in (0, 1):
        params = {
            "atom_labels": ["H", "H"],
            "es_strategy": strategies[itraj],
            "nstates": 2,
            "time_overlap": False,
            "gradient_state": "active",
            "act_state": act_state,
        }
        obj = strategy_compute_adi(q, params, full_id=[0, itraj])
        block = obj.d1ham_adi[0]
        expected_state = act_state[itraj]
        other_state = 1 - expected_state
        assert block.get(expected_state, expected_state).real == pytest.approx(
            expected_state + 1.0
        )
        assert block.get(other_state, other_state).real == pytest.approx(0.0)
