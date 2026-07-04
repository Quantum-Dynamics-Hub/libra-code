from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[4]
SRC = ROOT / "src"

if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from libra_py.dyn.initialization import (
    make_independent_trajectory_ensemble,
    sample_gaussian_initial_conditions,
)
from libra_py.dyn.observables import ObservableConfig, compute_observables
from libra_py.dyn.propagation import update_density


def test_sample_gaussian_initial_conditions_is_reproducible_and_shaped():
    rng = np.random.default_rng(12)
    q, p = sample_gaussian_initial_conditions(
        ntraj=4,
        q_mean=[-1.0, 2.0],
        q_sigma=[0.1, 0.2],
        p_mean=[5.0, -3.0],
        p_sigma=[0.5, 0.7],
        ndof=2,
        rng=rng,
    )

    assert q.shape == (4, 2)
    assert p.shape == (4, 2)
    assert np.allclose(q[0], [-1.00068268, 2.20922987])
    assert np.allclose(p[0], [4.94612375, -2.30086544])


def test_make_independent_trajectory_ensemble_normalizes_amplitudes():
    storage, traj = make_independent_trajectory_ensemble(
        q=[-1.0, 0.0, 1.0],
        p=[2.0, 2.0, 2.0],
        masses=4.0,
        amplitudes=[2.0 + 0.0j, 0.0 + 0.0j],
    )

    assert traj.tbf_ids == [0, 1, 2]
    assert storage.q.shape[:2] == (1, 3)
    np.testing.assert_allclose(storage.iM[0, :, 0], 0.25)
    np.testing.assert_allclose(storage.ampl_adi[0, :, 0], 1.0)
    np.testing.assert_allclose(storage.ampl_adi[0, :, 1], 0.0)


def test_initialized_ensemble_can_be_observed_with_observables_module():
    storage, traj = make_independent_trajectory_ensemble(
        q=[-1.0, 1.0],
        p=[2.0, 4.0],
        masses=2.0,
        amplitudes=[[1.0, 0.0], [0.0, 1.0]],
    )
    storage.ham_adi[0, 0] = np.diag([0.1, 0.2])
    storage.ham_adi[0, 1] = np.diag([0.3, 0.4])
    update_density(storage, traj)

    snapshot = compute_observables(
        storage,
        traj,
        ObservableConfig(
            potential="ehrenfest",
            active_counts=False,
            include_coordinates=True,
            include_momenta=True,
        ),
        time=0.5,
    )

    np.testing.assert_allclose(snapshot["mean_populations"], [0.5, 0.5])
    np.testing.assert_allclose(snapshot["kinetic_energy"], [1.0, 4.0])
    np.testing.assert_allclose(snapshot["potential_energy"], [0.1, 0.4])
    assert snapshot["q"].shape == (2, 1)
    assert snapshot["p"].shape == (2, 1)
