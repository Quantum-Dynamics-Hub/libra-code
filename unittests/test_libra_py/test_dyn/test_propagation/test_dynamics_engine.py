from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[4]
SRC = ROOT / "src"

if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from libra_py.dyn.backends import backend
from libra_py.dyn.control_params import DynControlParams
from libra_py.dyn.core.storage import TensorStorage
from libra_py.dyn.core.trajectory import Trajectory
from libra_py.dyn.engine import DynamicsEngine, dynamics_defaults
from libra_py.dyn.models import TullyModel1


def _storage():
    storage = TensorStorage(
        backend=backend,
        ntraj=1,
        ndof=1,
        nstates=2,
        ntbf_initial=1,
        ntbf_capacity=2,
    )
    storage.allocate_hamiltonian_derivatives()
    storage.iM[0, 0] = [1.0]
    storage.p[0, 0] = [1.0]
    storage.ampl_adi[0, 0] = [1.0 + 0.0j, 0.0 + 0.0j]
    storage.act_states[0, 0] = 1
    return storage


def _trajectory():
    traj = Trajectory(0)
    traj.tbf_ids = [0]
    return traj


def _adiabatic_model(R, P, storage, traj):
    return {
        "H_adi": np.array([np.diag([0.0, 0.2])], dtype=complex),
        "NAC_adi": np.zeros((1, 2, 2), dtype=complex),
        "dH_adi": np.array([[[[0.1, 0.0], [0.0, 0.3]]]], dtype=complex),
    }


def test_adiabatic_engine_uses_active_state_force_without_tdse():
    storage = _storage()
    traj = _trajectory()

    engine = DynamicsEngine(
        traj,
        storage,
        _adiabatic_model,
        method="adiabatic",
        rep="adiabatic",
    )
    result = engine.step(0.5)

    assert result.method == "adiabatic"
    assert result.timestep == 1
    np.testing.assert_allclose(result.forces, [[-0.3]])
    np.testing.assert_allclose(storage.p[0, 0], [0.85])
    np.testing.assert_allclose(storage.ampl_adi[0, 0], [1.0 + 0.0j, 0.0 + 0.0j])


def test_ehrenfest_engine_propagates_tdse_and_mean_field_force():
    storage = _storage()
    traj = _trajectory()
    storage.act_states[0, 0] = 0

    params = DynControlParams(force_method=2)
    engine = DynamicsEngine(
        traj,
        storage,
        _adiabatic_model,
        params=params,
        method="ehrenfest",
        rep="adiabatic",
    )
    result = engine.step(0.1)

    assert result.method == "ehrenfest"
    np.testing.assert_allclose(result.forces, [[-0.1]])
    np.testing.assert_allclose(abs(result.amplitudes[0, 0]), 1.0)
    np.testing.assert_allclose(storage.dm_adi[0, 0, 0, 0], 1.0 + 0.0j)


def test_dynamics_engine_tsh_runs_tdse_but_does_not_hop():
    storage = _storage()
    traj = _trajectory()
    storage.act_states[0, 0] = 1
    storage.ampl_adi[0, 0] = [0.0 + 0.0j, 1.0 + 0.0j]

    engine = DynamicsEngine(
        traj,
        storage,
        _adiabatic_model,
        method="tsh",
        rep="adiabatic",
    )
    result = engine.step(0.1)

    assert result.method == "tsh"
    np.testing.assert_array_equal(storage.act_states[0, 0:1], [1])
    np.testing.assert_allclose(abs(result.amplitudes[0, 1]), 1.0)


def _identity_propagator(coefficients, hamiltonian, dt, backend=None):
    del hamiltonian, dt, backend
    return coefficients


def _hopping_model(R, P, storage, traj):
    del R, P, storage, traj
    return {
        "H_adi": np.array([np.diag([0.0, 0.2])], dtype=complex),
        "Hvib_adi": np.array([[[0.0, -5.0j], [5.0j, 0.2]]]),
        "NAC_adi": np.zeros((1, 2, 2), dtype=complex),
        "dH_adi": np.array([[[[0.1, 0.0], [0.0, 0.3]]]], dtype=complex),
    }


def test_tsh_proposes_accepts_and_applies_active_state_force():
    storage = _storage()
    traj = _trajectory()
    storage.act_states[0, 0] = 0
    storage.ampl_adi[0, 0] = [2**-0.5, 2**-0.5]
    params = DynControlParams(tsh_method=0, dt=0.1, hop_acceptance_algo=0)

    result = DynamicsEngine(
        traj,
        storage,
        _hopping_model,
        params=params,
        propagator=_identity_propagator,
    ).step(0.1)

    np.testing.assert_allclose(result.hopping_probabilities, [[0.0, 1.0]])
    np.testing.assert_array_equal(result.proposed_states, [1])
    np.testing.assert_array_equal(result.accepted_states, [1])
    np.testing.assert_array_equal(storage.act_states[0, 0:1], [1])
    np.testing.assert_allclose(result.forces, [[-0.3]])


def test_tsh_energy_acceptance_can_reject_proposed_hop():
    storage = _storage()
    traj = _trajectory()
    storage.p[0, 0] = [0.1]
    storage.act_states[0, 0] = 0
    storage.ampl_adi[0, 0] = [2**-0.5, 2**-0.5]
    params = DynControlParams(tsh_method=0, dt=0.1, hop_acceptance_algo=10)

    result = DynamicsEngine(
        traj,
        storage,
        _hopping_model,
        params=params,
        propagator=_identity_propagator,
    ).step(0.1)

    np.testing.assert_array_equal(result.proposed_states, [1])
    np.testing.assert_array_equal(result.accepted_states, [0])
    np.testing.assert_array_equal(storage.act_states[0, 0:1], [0])


def test_engine_reports_edc_rates_from_pre_hop_phase():
    storage = _storage()
    traj = _trajectory()
    storage.act_states[0, 0] = 0
    storage.ampl_adi[0, 0] = [1.0, 0.0]
    params = DynControlParams(
        tsh_method=0,
        decoherence_times_type=1,
        decoherence_C_param=1.0,
        decoherence_eps_param=0.1,
    )
    result = DynamicsEngine(
        traj, storage, _adiabatic_model, params=params,
        propagator=_identity_propagator,
    ).step(0.1)

    # Ekin is close to 1/2 here; the important contract is that the rate phase
    # ran and exposed a symmetric, positive two-state rate matrix.
    assert result.decoherence_rates.shape == (1, 2, 2)
    assert result.decoherence_rates[0, 0, 1] > 0.0
    np.testing.assert_allclose(result.decoherence_rates,
                               result.decoherence_rates.swapaxes(1, 2))


def test_sdm_pre_hop_hook_can_be_used_with_ehrenfest_evolution():
    storage = _storage()
    traj = _trajectory()
    storage.act_states[0, 0] = 0
    storage.ampl_adi[0, 0] = [np.sqrt(0.6), 1j*np.sqrt(0.4)]
    params = DynControlParams(
        force_method=2,
        decoherence_algo=0,
        decoherence_times_type=0,
        decoherence_rates=np.array([[0.0, 1.0], [1.0, 0.0]]),
    )
    result = DynamicsEngine(
        traj, storage, _adiabatic_model, params=params, method="ehrenfest",
        propagator=_identity_propagator,
    ).step(0.5)

    np.testing.assert_allclose(abs(result.amplitudes[0, 1]),
                               np.sqrt(0.4)*np.exp(-0.5))
    np.testing.assert_allclose(np.vdot(result.amplitudes[0],
                                       result.amplitudes[0]).real, 1.0)


def test_instantaneous_decoherence_uses_rejected_hop_outcome():
    storage = _storage()
    traj = _trajectory()
    storage.p[0, 0] = [0.1]
    storage.act_states[0, 0] = 0
    storage.ampl_adi[0, 0] = [2**-0.5, 2**-0.5]
    params = DynControlParams(
        tsh_method=0, dt=0.1, hop_acceptance_algo=10,
        decoherence_algo=1, instantaneous_decoherence_variant=1,
    )
    result = DynamicsEngine(
        traj, storage, _hopping_model, params=params,
        propagator=_identity_propagator,
    ).step(0.1)

    np.testing.assert_array_equal(result.proposed_states, [1])
    np.testing.assert_array_equal(result.accepted_states, [0])
    np.testing.assert_allclose(result.amplitudes, [[1.0, 0.0]])


def test_post_hop_energy_defect_correction_is_exact():
    storage = _storage()
    traj = _trajectory()
    storage.act_states[0, 0] = 1
    storage.p[0, 0] = [0.9]
    storage.ham_adi[0, 0] = np.diag([0.0, 0.2])
    storage.dc1_adi[0, 0, 0, 0, 1] = 1.0
    params = DynControlParams(tsh_method=0, momenta_rescaling_algo=200)
    engine = DynamicsEngine(traj, storage, _adiabatic_model, params=params)
    engine._step_energy_before_hops = np.array([0.5])

    engine._correct_hop_energy_defect(np.array([1]), np.array([0]))

    np.testing.assert_allclose(engine._trajectory_total_energies(), [0.5],
                               atol=1.0e-14)


def test_engine_builds_time_overlap_for_local_diabatization_integrator():
    storage = _storage()
    traj = _trajectory()
    storage.q[0, 0] = [-0.2]
    storage.p[0, 0] = [0.4]
    storage.act_states[0, 0] = 0
    storage.ampl_adi[0, 0] = [1.0, 0.0]
    params = DynControlParams(
        force_method=0,
        electronic_integrator=0,
        time_overlap_method=1,
        state_tracking_algo=-1,
    )
    engine = DynamicsEngine(
        traj, storage, TullyModel1(), params=params, method="ehrenfest",
        force_mode="none",
    )
    engine.step(0.5)

    old_u = engine._previous_hamiltonian["basis_transform"]
    new_u = storage.basis_transform[0, traj.tbf_ids]
    expected = np.matmul(np.swapaxes(old_u.conj(), -1, -2), new_u)
    np.testing.assert_allclose(storage.time_overlap_adi[0, traj.tbf_ids], expected)
    assert not np.allclose(storage.proj_adi[0, traj.tbf_ids], np.eye(2))
    assert abs(storage.ampl_adi[0, 0, 1]) > 0.0


def test_engine_uses_electronic_substeps_and_legacy_defaults():
    storage = _storage()
    traj = _trajectory()
    calls = []

    def propagator(coefficients, hamiltonian, dt, backend=None):
        del hamiltonian, backend
        calls.append(dt)
        return coefficients

    params = DynControlParams(
        force_method=2,
        dt=0.3,
        num_electronic_substeps=3,
    )
    DynamicsEngine(
        traj,
        storage,
        _adiabatic_model,
        params=params,
        method="ehrenfest",
        propagator=propagator,
    ).step(0.3)

    # C++ adiabatic method 0 applies an old and a new Hamiltonian half-step
    # during each electronic substep.
    np.testing.assert_allclose(calls, [0.05] * 6)
    defaults = dynamics_defaults({"tsh_method": 0})
    assert defaults["tsh_method"] == 0
    assert defaults["rep_tdse"] == 1
    assert defaults["rep_sh"] == 1
    assert defaults["nprint"] == 1
    assert "states" in defaults["properties_to_save"]


def test_two_point_electronic_integrator_uses_old_and_new_hvib():
    storage = _storage()
    traj = _trajectory()
    storage.q[0, 0] = [0.0]
    storage.p[0, 0] = [1.0]

    def coordinate_model(R, P, storage, traj):
        del P, storage, traj
        q = float(R[0, 0])
        return {
            "H_adi": np.array([np.diag([q, 0.0])], dtype=complex),
            "NAC_adi": np.zeros((1, 2, 2), dtype=complex),
            "dH_adi": np.zeros((1, 1, 2, 2), dtype=complex),
        }

    result = DynamicsEngine(
        traj,
        storage,
        coordinate_model,
        params=DynControlParams(force_method=0, electronic_integrator=4),
        method="ehrenfest",
        force_mode="none",
    ).step(1.0)

    # Symmetric splitting: exp(-i H_new dt/2) exp(-i H_old dt/2).
    np.testing.assert_allclose(result.amplitudes[0, 0], np.exp(-0.5j))


def test_local_diabatization_reprojects_amplitude_and_active_state():
    storage = _storage()
    traj = _trajectory()
    storage.act_states[0, 0] = 0
    storage.ampl_adi[0, 0] = [1.0, 0.0]

    def model(R, P, storage, traj):
        del R, P, storage, traj
        return {
            "H_adi": np.array([np.diag([0.0, 0.2])], dtype=complex),
            "NAC_adi": np.zeros((1, 2, 2), dtype=complex),
            "dH_adi": np.array([[[[0.1, 0.0], [0.0, 0.3]]]], dtype=complex),
            "time_overlap_adi": np.array([[[0.0, 1.0], [1.0, 0.0]]]),
        }

    result = DynamicsEngine(
        traj,
        storage,
        model,
        params=DynControlParams(tsh_method=-1, electronic_integrator=0),
        method="tsh",
        propagator=_identity_propagator,
    ).step(0.1)

    np.testing.assert_allclose(abs(result.amplitudes[0]), [0.0, 1.0])
    np.testing.assert_array_equal(storage.act_states[0, 0:1], [1])


def test_engine_saver_receives_initial_and_stride_snapshots():
    storage = _storage()
    traj = _trajectory()

    class Saver:
        def __init__(self):
            self.steps = []

        def save_observables(self, storage, traj, config, step, time, metadata):
            del storage, traj, config, time, metadata
            self.steps.append(step)

    saver = Saver()
    engine = DynamicsEngine(
        traj,
        storage,
        _adiabatic_model,
        params=DynControlParams(nprint=2),
        method="adiabatic",
        saver=saver,
    )
    engine.run(3, 0.1)

    assert saver.steps == [0, 2]
