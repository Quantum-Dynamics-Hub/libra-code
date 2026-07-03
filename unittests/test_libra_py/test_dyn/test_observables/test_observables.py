from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[4]
SRC = ROOT / "src"

if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from libra_py.dyn.backends import backend
from libra_py.dyn.core.storage import TensorStorage
from libra_py.dyn.core.trajectory import Trajectory
from libra_py.dyn.observables import (
    ObservableConfig,
    active_potential_energy,
    active_state_counts,
    amplitude_populations,
    compute_observables,
    kinetic_energy,
    legacy_observable_keywords,
    total_energy,
)


def _system():
    storage = TensorStorage(
        backend=backend,
        ntraj=1,
        ndof=2,
        nstates=2,
        ntbf_initial=1,
    )
    traj = Trajectory(0)
    traj.tbf_ids = [0]

    storage.p[0, 0] = [2.0, -1.0]
    storage.iM[0, 0] = [0.5, 0.25]
    storage.ampl_adi[0, 0] = [0.6 + 0.0j, 0.8 + 0.0j]
    storage.ampl_dia[0, 0] = [1.0 + 0.0j, 0.0 + 0.0j]
    storage.act_states[0, 0] = 1
    storage.act_states_dia[0, 0] = 0
    storage.ham_adi[0, 0] = np.diag([0.1, 0.4])
    storage.hvib_adi[0, 0] = np.diag([0.1, 0.4])
    return storage, traj


def test_population_and_energy_observables():
    storage, traj = _system()

    np.testing.assert_allclose(amplitude_populations(storage, traj), [[0.36, 0.64]])
    np.testing.assert_array_equal(active_state_counts(storage, traj), [0, 1])
    np.testing.assert_allclose(kinetic_energy(storage, traj), [1.125])
    np.testing.assert_allclose(active_potential_energy(storage, traj), [0.4])
    np.testing.assert_allclose(total_energy(storage, traj), [1.525])


def test_compute_observables_respects_config():
    storage, traj = _system()
    config = ObservableConfig(
        populations=True,
        active_counts=False,
        energies=False,
        include_coordinates=True,
    )

    data = compute_observables(storage, traj, config, step=7, time=0.35)

    assert data["step"] == 7
    assert data["time"] == 0.35
    assert "populations" in data
    assert "active_state_counts" not in data
    assert "kinetic_energy" not in data
    assert "q" in data


def test_legacy_observable_keywords_follow_save_levels():
    names = legacy_observable_keywords(4)

    assert "Ekin_ave" in names
    assert "se_pop_adi" in names
    assert "D_adi" in names
    assert "hvib_adi" in names
    assert "dc1_adi" not in names


def test_compute_observables_supports_old_save_keyword_names():
    storage, traj = _system()
    data = compute_observables(
        storage,
        traj,
        ObservableConfig(
            keywords=(
                "timestep",
                "time",
                "Ekin_ave",
                "dEkin_ave",
                "Epot_ave",
                "Etot_ave",
                "states",
                "se_pop_adi",
                "sh_pop_adi",
                "D_adi",
                "coherence_adi",
                "q",
                "Cadi",
                "hvib_adi",
                "energy_gaps",
            )
        ),
        step=4,
        time=0.2,
    )

    assert data["timestep"] == 4
    assert data["time"] == 0.2
    assert data["Ekin_ave"] == 1.125
    assert data["dEkin_ave"] == 0.0
    assert data["Epot_ave"] == 0.4
    assert data["Etot_ave"] == 1.525
    np.testing.assert_array_equal(data["states"], [1])
    np.testing.assert_allclose(data["se_pop_adi"], [0.36, 0.64])
    np.testing.assert_allclose(data["sh_pop_adi"], [0.0, 1.0])
    np.testing.assert_allclose(data["D_adi"], [[0.36, 0.48], [0.48, 0.64]])
    np.testing.assert_allclose(data["coherence_adi"], [[0.0, 1.0], [1.0, 0.0]])
    np.testing.assert_allclose(data["q"], [[0.0, 0.0]])
    np.testing.assert_allclose(data["Cadi"], [[0.6 + 0.0j, 0.8 + 0.0j]])
    np.testing.assert_allclose(data["hvib_adi"], [[[0.1, 0.0], [0.0, 0.4]]])
    np.testing.assert_allclose(data["energy_gaps"], [[[0.0, -0.3], [0.3, 0.0]]])
