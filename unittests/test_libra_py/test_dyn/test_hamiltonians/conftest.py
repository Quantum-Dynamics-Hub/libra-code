from pathlib import Path
import sys

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[4]
SRC = ROOT / "src"

if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from libra_py.dyn.backends import backend
from libra_py.dyn.core.storage import TensorStorage
from libra_py.dyn.core.trajectory import Trajectory


@pytest.fixture
def hamiltonian_storage():
    storage = TensorStorage(
        backend=backend,
        ntraj=1,
        ndof=1,
        nstates=2,
        ntbf_initial=1,
        ntbf_capacity=2,
    )
    storage.allocate_hamiltonian_derivatives(der_lvl=2)
    storage.q[0, 0] = [0.25]
    storage.p[0, 0] = [0.5]
    storage.ampl_dia[0, 0] = [1.0 + 0.0j, 0.0 + 0.0j]
    storage.ampl_adi[0, 0] = [1.0 + 0.0j, 0.0 + 0.0j]
    return storage


@pytest.fixture
def trajectory():
    traj = Trajectory(0)
    traj.tbf_ids = [0]
    return traj


def diabatic_model(R, P, storage, traj, params=None):
    return {
        "ham_dia": np.array(
            [
                [
                    [0.0, 0.01],
                    [0.01, 0.1],
                ]
            ],
            dtype=complex,
        ),
        "ovlp_dia": np.array([np.eye(2)], dtype=complex),
        "nac_dia": np.array(
            [
                [
                    [0.0, 0.02],
                    [-0.02, 0.0],
                ]
            ],
            dtype=complex,
        ),
        "d1ham_dia": np.array(
            [
                [
                    [
                        [0.1, 0.0],
                        [0.0, 0.2],
                    ]
                ]
            ],
            dtype=complex,
        ),
        "dc1_dia": np.zeros((1, 1, 2, 2), dtype=complex),
        "d2ham_dia": np.zeros((1, 1, 1, 2, 2), dtype=complex),
    }


def adiabatic_model(R, P, storage, traj, params=None):
    return {
        "H_adi": np.array([np.diag([0.0, 0.2])], dtype=complex),
        "NAC_adi": np.array(
            [
                [
                    [0.0, 0.03],
                    [-0.03, 0.0],
                ]
            ],
            dtype=complex,
        ),
        "dH_adi": np.array(
            [
                [
                    [
                        [0.1, 0.0],
                        [0.0, 0.2],
                    ]
                ]
            ],
            dtype=complex,
        ),
    }


@pytest.fixture
def diabatic_model_fn():
    return diabatic_model


@pytest.fixture
def adiabatic_model_fn():
    return adiabatic_model
