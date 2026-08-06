from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Callable, Literal

import numpy as np

from .control_params import DynControlParams
from .decoherence import (
    dephasing_informed_correction,
    dish_hop_proposal,
    dish_project_out_collapse,
    dish_rev2023,
    edc_rates,
    gu_franco,
    instantaneous_decoherence,
    schwartz_1,
    schwartz_1_interaction_width,
    schwartz_2,
    sdm,
)
from .hamiltonians import HamiltonianEngine
from .hopping import (
    accept_hops,
    handle_hops_nuclear,
    hop_proposal_probabilities,
    propose_hops,
    rescale_along_vector,
)
from .observables import ObservableConfig
from .propagation.coupled import (
    ehrenfest_forces,
    state_specific_forces,
    update_density,
)
from .propagation.electronic import exp_propagator, tdse_step
from .propagation.integrators import normalize_dt, run_steps
from .propagation.nuclear import drift, kick
from .transformations.basis_rotation import (
    storage_amplitudes_adi_to_dia,
    storage_amplitudes_dia_to_adi,
)
from .transformations.projectors import update_proj_adi

DynamicsMethod = Literal["adiabatic", "ehrenfest", "tsh"]
Representation = Literal["adiabatic", "diabatic"]
ForceMode = Literal["state-specific", "ehrenfest", "none"]


@dataclass
class StepResult:
    """Summary of one storage-backed dynamics step."""

    method: str
    rep: str
    time: float
    timestep: int
    active_states: Any = None
    forces: Any = None
    amplitudes: Any = None
    hopping_probabilities: Any = None
    proposed_states: Any = None
    accepted_states: Any = None
    decoherence_rates: Any = None


@dataclass
class DynamicsEngine:
    """
    Storage-backed nonadiabatic dynamics driver.

    Supported methods:
    - ``adiabatic``: nuclei follow the active adiabatic state; no TDSE update.
    - ``ehrenfest``: TDSE plus Ehrenfest mean-field nuclear forces.
    - ``tsh``: TDSE, hop proposal, acceptance, momentum adjustment, and
      state-specific forces using the legacy ``tsh_method`` option numbering.
    """

    traj: Any
    storage: Any
    model_fn: Callable
    ham_engine: HamiltonianEngine | None = None
    params: DynControlParams | dict[str, Any] | None = None
    method: DynamicsMethod | None = None
    rep: Representation | None = None
    force_mode: ForceMode | None = None
    hamiltonian_type: str = "vibronic"
    propagator: Callable = exp_propagator
    rng: Any = None
    time: float = 0.0
    metadata: dict[str, Any] = field(default_factory=dict)
    saver: Any = None
    observable_config: Any = None

    def __post_init__(self):
        self.params = self.params or DynControlParams()
        self.ham = self.ham_engine or HamiltonianEngine(self.storage.backend)
        self.rng = self.rng or np.random.default_rng()
        self.method = self.method or self._method_from_params()
        self.rep = self.rep or self._rep_from_params("tdse")
        self.force_mode = self.force_mode or self._force_mode_from_params()
        self.observable_config = self.observable_config or self._observable_config()
        self._previous_hamiltonian = None
        self._pending_hop_momentum_update = None
        self._initialized = False
        self._validate()

    def initialize(self, evaluate_hamiltonian: bool = True):
        """Build the initial Hamiltonian, forces, and density matrices."""

        if evaluate_hamiltonian:
            self.evaluate_hamiltonian()
        self._initialize_representations()
        if self._uses_tdse:
            update_density(self.storage, self.traj, self.rep)
        if self.method == "tsh":
            self._allocate_tsh_storage()
        self._allocate_decoherence_storage()
        self.compute_forces()
        self._setup_saver()
        self._initialized = True
        self._save_if_due(force=True)
        return self

    def step(self, dt) -> StepResult:
        """Advance one mixed quantum-classical dynamics step."""

        dt = normalize_dt(dt)
        self._step_dt = dt
        if not self._initialized:
            self.initialize()
        self._step_energy_before_hops = self._trajectory_total_energies()
        kick(self.storage, self.traj, 0.5 * dt)
        drift(self.storage, self.traj, dt)

        self._save_electronic_history()
        self._previous_hamiltonian = self._hamiltonian_snapshot()
        self._prepare_time_overlap_update()
        self.evaluate_hamiltonian()
        self._update_time_overlaps()
        ld_transform = self._state_tracking_transform()
        if self._uses_tdse:
            self.propagate_electronic(dt, T=ld_transform)
            self._apply_state_tracking_to_active_states(ld_transform)
            self._synchronize_representations()
            update_density(self.storage, self.traj, self.rep)

        decoherence_rates = None
        decoherence_algo = int(_param(self.params, "decoherence_algo", -1))
        needs_decoherence = self.method == "tsh" or (
            self._uses_tdse and decoherence_algo in (0, 3, 4, 5, 6, 7, 9)
        )
        if needs_decoherence:
            # compute_dynamics updates rates after coherent electronic
            # propagation and immediately before pre-hop corrections.
            decoherence_rates = self._compute_decoherence_rates()
            self._apply_pre_hop_decoherence(dt, decoherence_rates)

        hopping = None
        proposed = None
        accepted = None
        if self.method == "tsh":
            hopping, proposed, accepted = self.surface_hopping_step(decoherence_rates)
        forces = self.compute_forces()
        kick(self.storage, self.traj, 0.5 * dt, forces)
        self._finalize_hop_momenta()

        self.evaluate_hamiltonian()
        if self._uses_tdse:
            update_density(self.storage, self.traj, self.rep)

        self.time += dt
        self.storage.timestep += 1
        self._save_if_due()
        return self._result(forces, hopping, proposed, accepted,
                            decoherence_rates)

    def run(self, nsteps: int | None = None, dt=None):
        """Run several steps and return their summaries."""

        nsteps = int(_param(self.params, "nsteps", 1) if nsteps is None else nsteps)
        dt = _param(self.params, "dt", 41.0) if dt is None else dt
        return run_steps(self.step, nsteps, dt)

    def evaluate_hamiltonian(self):
        """Evaluate the Hamiltonian model into TensorStorage."""

        return self.ham.build_state(
            self.traj,
            self.storage,
            self.model_fn,
            rep=self.rep,
            apply_ssy=bool(_param(self.params, "do_ssy", 0)),
        )

    def propagate_electronic(self, dt, T=None):
        """Propagate active amplitudes with the TD-SE solver."""

        integrator = int(_param(self.params, "electronic_integrator", 0))
        if integrator == -1:
            return self._active_amplitudes()
        nsubsteps = int(_param(self.params, "num_electronic_substeps", 1))
        if nsubsteps < 1:
            raise ValueError("num_electronic_substeps must be positive")
        if self._previous_hamiltonian is None:
            raise RuntimeError("electronic propagation requires the previous Hamiltonian")
        result = None
        for substep in range(nsubsteps):
            substep_dt = dt / nsubsteps
            result = tdse_step(
                self.traj,
                self.storage,
                substep_dt,
                backend=self.storage.backend,
                propagator=self.propagator,
                rep=self.rep,
                hamiltonian_type=self.hamiltonian_type,
                T=T,
                previous_state=self._previous_hamiltonian,
                method=integrator,
            )
        return result

    def surface_hopping_step(self, decoherence_rates=None):
        """Propose, accept, and handle hops for active TBFs."""

        idx = self.traj.tbf_ids
        rep_sh = self._rep_from_params("sh")
        density = self._storage_slice("dm", rep_sh)
        hvib = self._storage_slice("hvib", rep_sh)
        states_field = "act_states" if rep_sh == "adiabatic" else "act_states_dia"
        initial = np.array(getattr(self.storage, states_field)[self.traj.id, idx], copy=True)
        previous_density = self._previous_density(rep_sh)
        method = int(_param(self.params, "tsh_method", -1))

        if method == 5:
            if int(_param(self.params, "decoherence_algo", -1)) != -1:
                raise ValueError("C++ convention requires tsh_method=5 (DISH) with decoherence_algo=-1")
            clocks = np.array(self.storage.coherence_time[self.traj.id, idx], copy=True)
            clocks += self._step_dt
            amplitudes = np.array(self.storage.ampl_adi[self.traj.id, idx], copy=True)
            proposed = dish_hop_proposal(
                initial, amplitudes, clocks,
                decoherence_rates, self.rng,
                int(_param(self.params, "dish_decoherence_event_option", 1)),
            )
            accepted = self._accept_hops(proposed, initial)
            dish_project_out_collapse(
                initial, proposed, accepted,
                amplitudes,
                int(_param(self.params, "collapse_option", 0)),
            )
            self.storage.coherence_time[self.traj.id, idx] = clocks
            self.storage.ampl_adi[self.traj.id, idx] = amplitudes
            self._pending_hop_momentum_update = (accepted.copy(), initial.copy())
            getattr(self.storage, states_field)[self.traj.id, idx] = accepted
            self._map_active_states(rep_sh, accepted)
            self._synchronize_representations()
            return None, proposed, accepted

        kwargs = {}
        if method in (3, 4):
            current_records, previous_records = self._hop_hamiltonian_records()
            kwargs.update(
                ham=current_records,
                ham_prev=previous_records,
                momentum=np.asarray(self.storage.p[self.traj.id, idx]),
                inverse_mass=np.asarray(self.storage.iM[self.traj.id, idx][0]),
            )
        probabilities = hop_proposal_probabilities(
            self._proposal_params(),
            density,
            hvib,
            initial,
            previous_density,
            fssh3_errors=None if self.storage.fssh3_errors is None else self.storage.fssh3_errors[self.traj.id, idx],
            **kwargs,
        )
        proposed = propose_hops(probabilities, initial, self.rng)
        accepted = self._accept_hops(proposed, initial)

        # In compute_dynamics, corrections which depend on the hop outcome
        # precede velocity rescaling and active-state assignment.
        self._apply_post_hop_decoherence(initial, proposed, accepted)
        self._pending_hop_momentum_update = (accepted.copy(), initial.copy())
        getattr(self.storage, states_field)[self.traj.id, idx] = accepted
        self._map_active_states(rep_sh, accepted)
        return probabilities, proposed, accepted

    def compute_forces(self):
        """Compute and store active nuclear forces."""

        idx = self.traj.tbf_ids
        if self.force_mode == "none":
            forces = np.zeros_like(self.storage.f[self.traj.id, idx])
        elif self.force_mode == "ehrenfest":
            forces = ehrenfest_forces(
                self.storage,
                self.traj,
                rep=self.rep,
                option=int(_param(self.params, "ehrenfest_force_option", 0)),
                gamma=_param(self.params, "sqc_gamma", 0.0),
            )
        else:
            forces = state_specific_forces(self.storage, self.traj, rep="adiabatic")

        self.storage.f[self.traj.id, idx] = forces
        return self.storage.f[self.traj.id, idx]

    @property
    def _uses_tdse(self) -> bool:
        return self.method in ("ehrenfest", "tsh")

    def _result(self, forces, hopping=None, proposed=None, accepted=None,
                decoherence_rates=None):
        idx = self.traj.tbf_ids
        amplitudes = (
            self.storage.ampl_adi[self.traj.id, idx]
            if self.rep == "adiabatic"
            else self.storage.ampl_dia[self.traj.id, idx]
        )
        active = (
            self.storage.act_states[self.traj.id, idx]
            if self.rep == "adiabatic"
            else self.storage.act_states_dia[self.traj.id, idx]
        )
        return StepResult(
            method=self.method,
            rep=self.rep,
            time=self.time,
            timestep=self.storage.timestep,
            active_states=np.array(active, copy=True),
            forces=np.array(forces, copy=True),
            amplitudes=np.array(amplitudes, copy=True),
            hopping_probabilities=None if hopping is None else np.array(hopping, copy=True),
            proposed_states=None if proposed is None else np.array(proposed, copy=True),
            accepted_states=None if accepted is None else np.array(accepted, copy=True),
            decoherence_rates=(None if decoherence_rates is None else
                               np.array(decoherence_rates, copy=True)),
        )

    def _method_from_params(self) -> DynamicsMethod:
        force_method = int(_param(self.params, "force_method", 1))
        tsh_method = int(_param(self.params, "tsh_method", -1))
        if tsh_method >= 0:
            return "tsh"
        if force_method == 2:
            return "ehrenfest"
        return "adiabatic"

    def _force_mode_from_params(self) -> ForceMode:
        if self.method == "ehrenfest":
            return "ehrenfest"
        force_method = int(_param(self.params, "force_method", 1))
        if force_method == 0:
            return "none"
        if force_method == 2:
            return "ehrenfest"
        return "state-specific"

    def _rep_from_params(self, use: str) -> Representation:
        names = {"tdse": "rep_tdse", "force": "rep_force", "sh": "rep_sh"}
        name = names[use]
        return "adiabatic" if int(_param(self.params, name, 1)) in (1, 3, 4) else "diabatic"

    def _initialize_representations(self):
        idx = self.traj.tbf_ids
        if not np.any(self.storage.ovlp_dia[self.traj.id, idx]):
            eye = np.eye(self.storage.nstates, dtype=complex)
            self.storage.ovlp_dia[self.traj.id, idx] = eye
        if not np.any(self.storage.proj_adi[self.traj.id, idx]):
            self.storage.proj_adi[self.traj.id, idx] = np.eye(self.storage.nstates)
        if not np.any(self.storage.basis_transform[self.traj.id, idx]):
            self.storage.basis_transform[self.traj.id, idx] = np.eye(self.storage.nstates)
        self._synchronize_representations()

    def _synchronize_representations(self):
        if self.rep == "adiabatic":
            storage_amplitudes_adi_to_dia(self.storage, self.traj)
        else:
            storage_amplitudes_dia_to_adi(self.storage, self.traj)
        update_density(self.storage, self.traj, "adiabatic")
        update_density(self.storage, self.traj, "diabatic")

    def _allocate_tsh_storage(self):
        method = int(_param(self.params, "tsh_method", -1))
        if method in (7, 8, 9) and self.storage.dm_adi_prev is None:
            self.storage.allocate_fssh2()
        if method == 8 and self.storage.fssh3_errors is None:
            self.storage.allocate_fssh3()
        self._save_electronic_history()

    def _allocate_decoherence_storage(self):
        """Allocate method-specific state, mirroring C++ dyn_variables setup."""
        algo = int(_param(self.params, "decoherence_algo", -1))
        method = int(_param(self.params, "tsh_method", -1))
        if (method == 5 or algo == 7) and self.storage.coherence_time is None:
            self.storage.allocate_dish()
        if algo == 2 and self.storage.dR is None:
            self.storage.allocate_afssh()
        if algo == 3 and self.storage.reversal_events is None:
            self.storage.allocate_bcsh()
        if algo in (5, 6) and self.storage.q_aux is None:
            self.storage.allocate_shxf()
        if algo == 6 and self.storage.f_xf is None:
            self.storage.allocate_mqcxf()
        if algo == 9 and self.storage.coherence_factors is None:
            self.storage.allocate_simple_decoherence()

    def _state_forces_for_decoherence(self):
        """Return adiabatic state forces as ``(trajectory,state,dof)``."""
        d1 = np.asarray(self.storage.d1ham_adi[self.traj.id, self.traj.tbf_ids])
        return -np.diagonal(d1.real, axis1=-2, axis2=-1).swapaxes(1, 2)

    def _compute_decoherence_rates(self):
        """C++ ``compute_dynamics`` phase: update rates before TSH.

        Rate construction is independent of the decoherence algorithm. This
        lets SDM, DISH, MFSD, and future methods share the same time model.
        """
        idx = self.traj.tbf_ids
        c = np.asarray(self.storage.ampl_adi[self.traj.id, idx])
        ntraj, nstates = c.shape
        option = int(_param(self.params, "decoherence_times_type", -1))
        if option == -1:
            rates = np.zeros((ntraj, nstates, nstates))
        elif option == 0:
            supplied = _param(self.params, "decoherence_rates", None)
            if supplied is None:
                raise ValueError("decoherence_times_type=0 requires decoherence_rates")
            rates = np.broadcast_to(np.asarray(supplied, float), (ntraj, nstates, nstates)).copy()
        elif option == 1:
            p = np.asarray(self.storage.p[self.traj.id, idx])
            im = np.asarray(self.storage.iM[self.traj.id, idx])
            if im.ndim == 3: im = im[:, 0]
            ekin = 0.5 * np.sum(p * p * im, axis=1)
            rates = edc_rates(self.storage.hvib_adi[self.traj.id, idx], ekin,
                              _param(self.params, "decoherence_C_param", 1.0),
                              _param(self.params, "decoherence_eps_param", 0.1))
        elif option in (2, 3, 4):
            forces = self._state_forces_for_decoherence()
            if option == 2:
                rates = schwartz_1(c, forces, _param(self.params, "schwartz_decoherence_inv_alpha"))
            elif option == 3:
                rates = schwartz_2(forces, _param(self.params, "schwartz_decoherence_inv_alpha"))
            else:
                rates = schwartz_1_interaction_width(
                    c, forces, self.storage.p[self.traj.id, idx],
                    _param(self.params, "schwartz_interaction_width"))
        elif option == 5:
            rates = gu_franco(c, _param(self.params, "reorg_energy", 0.0),
                              _param(self.params, "Temperature", 300.0))
        else:
            raise ValueError(f"unknown decoherence_times_type={option}")
        if int(_param(self.params, "dephasing_informed", 0)):
            average = _param(self.params, "ave_gaps", None)
            if average is None: raise ValueError("dephasing_informed=1 requires ave_gaps")
            rates = dephasing_informed_correction(
                rates, self.storage.hvib_adi[self.traj.id, idx], average)
        self._decoherence_rates = rates
        return rates

    def _apply_pre_hop_decoherence(self, dt, rates):
        """C++ phase: apply corrections which precede hop proposal.

        ``compute_dynamics`` places SDM, BCSH, MFSD, SHXF/MQCXF, revised
        DISH, and simple decoherence here. Instantaneous decoherence and AFSSH
        wait until the proposed and accepted states are known.
        """
        algo = int(_param(self.params, "decoherence_algo", -1))
        if algo == -1 or not self._uses_tdse: return
        if self.rep != "adiabatic":
            raise ValueError("the selected decoherence algorithm requires rep_tdse=1")
        idx=self.traj.tbf_ids
        c=np.array(self.storage.ampl_adi[self.traj.id,idx], copy=True)
        states=self.storage.act_states[self.traj.id,idx]
        if algo == 0:
            c[...] = sdm(c,dt,states,rates,_param(self.params,"sdm_norm_tolerance",0.0))
        elif algo == 7:
            clocks=np.array(self.storage.coherence_time[self.traj.id,idx], copy=True)
            dish_rev2023(c,states,clocks,rates,dt,
                         int(_param(self.params,"decoherence_times_type",-1)),
                         int(_param(self.params,"dish_decoherence_event_option",1)),
                         int(_param(self.params,"collapse_option",0)),self.rng)
            self.storage.coherence_time[self.traj.id,idx] = clocks
        elif algo not in (1, 2, 8):
            raise NotImplementedError(
                f"decoherence_algo={algo} needs auxiliary nuclear propagation not yet available in DynamicsEngine")
        self.storage.ampl_adi[self.traj.id,idx] = c
        self._synchronize_representations()

    def _apply_post_hop_decoherence(self, initial, proposed, accepted):
        """C++ phase: apply corrections requiring the hopping outcome.

        The engine calls this after acceptance but before nuclear momentum
        handling, matching ``compute_dynamics``. AFSSH and XF reset operations
        belong in this hook when their auxiliary propagation is available.
        """
        algo=int(_param(self.params,"decoherence_algo",-1))
        if algo not in (1,2,5,6,8): return
        if algo == 2:
            raise NotImplementedError("AFSSH post-hop moments/collapse are not yet implemented")
        if algo in (5,6):
            raise NotImplementedError("SHXF/MQCXF auxiliary hop reset is not yet implemented")
        if algo == 8:
            raise NotImplementedError("diabatic instantaneous decoherence requires diabatic-state hop bookkeeping")
        idx=self.traj.tbf_ids
        c=np.array(self.storage.ampl_adi[self.traj.id,idx], copy=True)
        instantaneous_decoherence(c,accepted,proposed,initial,
            int(_param(self.params,"instantaneous_decoherence_variant",1)),
            int(_param(self.params,"collapse_option",0)))
        self.storage.ampl_adi[self.traj.id,idx] = c
        self._synchronize_representations()

    def _save_electronic_history(self):
        if self.storage.dm_adi_prev is None:
            return
        idx = self.traj.tbf_ids
        self.storage.dm_adi_prev[self.traj.id, idx] = self.storage.dm_adi[self.traj.id, idx]
        self.storage.dm_dia_prev[self.traj.id, idx] = self.storage.dm_dia[self.traj.id, idx]

    def _previous_density(self, rep):
        field = "dm_adi_prev" if rep == "adiabatic" else "dm_dia_prev"
        value = getattr(self.storage, field, None)
        return None if value is None else np.asarray(value[self.traj.id, self.traj.tbf_ids])

    def _state_tracking_transform(self):
        """Build ``T_new`` with columns ordered by old tracked state labels.

        ``time_overlap_adi`` uses ``S[old, new_raw]``. The returned projector
        is passed to electronic propagation, where ``T_new @ C_old`` produces
        new raw-basis coefficients and ``T_new† @ H_new @ T_new`` produces a
        Hamiltonian in the old dynamically consistent labels.
        """

        integrator = int(_param(self.params, "electronic_integrator", 0))
        base_integrator = integrator - 100 if 100 <= integrator < 200 else integrator
        if self.rep != "adiabatic" or base_integrator not in range(0, 16):
            return None
        if int(_param(self.params, "assume_always_consistent", 0)):
            return None
        overlap = np.asarray(self.storage.time_overlap_adi[self.traj.id, self.traj.tbf_ids])
        if not np.any(overlap):
            return None
        return update_proj_adi(
            self.params,
            self.storage,
            self.traj,
            previous=self._previous_hamiltonian,
            rng=self.rng,
        )

    def _update_time_overlaps(self):
        """Build ``S[old,new] = U_old† U_new`` as C++ option 1 does.

        Analytical models normally provide consecutive diabatic-to-adiabatic
        eigenvector matrices rather than explicit time overlaps. C++
        ``update_Hamiltonian_variables`` constructs the overlap from those
        matrices before updating projectors. The same phase is required here;
        without it, local-diabatization integrators 0--2 receive an identity
        projector and cannot describe transitions for models such as Tully-1.

        ``time_overlap_method=0`` leaves an externally supplied overlap
        untouched. Option 1 uses the orthonormal-basis expression implemented
        by C++.
        """
        if self.rep != "adiabatic" or int(_param(self.params, "time_overlap_method", 1)) == 0:
            return None
        idx = self.traj.tbf_ids
        supplied = np.asarray(self.storage.time_overlap_adi[self.traj.id, idx])
        if np.any(supplied):
            return supplied
        previous = None if self._previous_hamiltonian is None else self._previous_hamiltonian.get("basis_transform")
        if previous is None:
            return None
        current = np.asarray(self.storage.basis_transform[self.traj.id, idx])
        previous = np.asarray(previous)
        overlap = np.matmul(np.swapaxes(previous.conj(), -1, -2), current)
        self.storage.time_overlap_adi[self.traj.id, idx] = overlap
        return overlap

    def _prepare_time_overlap_update(self):
        """Clear stale overlaps before allowing the model to supply new ones."""
        if self.rep == "adiabatic" and int(_param(self.params, "time_overlap_method", 1)) == 1:
            self.storage.time_overlap_adi[self.traj.id, self.traj.tbf_ids] = 0.0

    def _proposal_params(self):
        if isinstance(self.params, dict):
            values = dict(self.params)
        else:
            values = asdict(self.params)
        values["dt"] = getattr(self, "_step_dt", values.get("dt", 41.0))
        return values

    def _active_amplitudes(self):
        field = "ampl_adi" if self.rep == "adiabatic" else "ampl_dia"
        return getattr(self.storage, field)[self.traj.id, self.traj.tbf_ids]

    def _apply_state_tracking_to_active_states(self, transform):
        """Map active old labels through projector columns.

        Column ``i`` of ``transform`` represents tracked old state ``i`` in
        the raw new basis. Its largest-magnitude row is therefore the new raw
        active-state label. Unit-modulus phase corrections do not change this
        mapping.
        """

        if transform is None or int(_param(self.params, "tsh_method", -1)) in (3, 4):
            return
        idx = self.traj.tbf_ids
        states = np.asarray(self.storage.act_states[self.traj.id, idx], dtype=int)
        mapped = np.argmax(
            np.abs(transform[np.arange(len(states)), :, states]), axis=1
        )
        self.storage.act_states[self.traj.id, idx] = mapped

    def _storage_slice(self, kind, rep):
        suffix = "adi" if rep == "adiabatic" else "dia"
        return np.asarray(getattr(self.storage, f"{kind}_{suffix}")[self.traj.id, self.traj.tbf_ids])

    def _hamiltonian_snapshot(self):
        idx = self.traj.tbf_ids
        fields = (
            "ham_dia", "ham_adi", "hvib_dia", "hvib_adi", "nac_adi",
            "d1ham_dia", "d1ham_adi", "dc1_adi", "basis_transform",
        )
        return {
            field: None if getattr(self.storage, field, None) is None else np.array(getattr(self.storage, field)[self.traj.id, idx], copy=True)
            for field in fields
        }

    def _hop_hamiltonian_records(self):
        if self._previous_hamiltonian is None:
            raise RuntimeError("LZ/ZN requires a previous Hamiltonian snapshot")
        current = self._hamiltonian_snapshot()
        records, previous = [], []
        for pos in range(len(self.traj.tbf_ids)):
            record = {name: value[pos] for name, value in current.items() if value is not None}
            old = {name: value[pos] for name, value in self._previous_hamiltonian.items() if value is not None}
            if "d1ham_adi" in record:
                record["forces_adi"] = -np.diagonal(record["d1ham_adi"].real, axis1=1, axis2=2)
            records.append(record)
            previous.append(old)
        return records, previous

    def _accept_hops(self, proposed, initial):
        idx = self.traj.tbf_ids
        return accept_hops(
            self.params,
            proposed,
            initial,
            np.asarray(self.storage.ham_adi[self.traj.id, idx]),
            momenta=np.asarray(self.storage.p[self.traj.id, idx]),
            inverse_mass=np.asarray(self.storage.iM[self.traj.id, idx][0]),
            energies_dia=np.asarray(self.storage.ham_dia[self.traj.id, idx]),
            dc1_adi=None if self.storage.dc1_adi is None else np.asarray(self.storage.dc1_adi[self.traj.id, idx]),
            d1ham_adi=None if self.storage.d1ham_adi is None else np.asarray(self.storage.d1ham_adi[self.traj.id, idx]),
            rng=self.rng,
        )

    def _rescale_hop_momenta(self, accepted, initial):
        idx = self.traj.tbf_ids
        momenta = self.storage.p[self.traj.id, idx]
        handle_hops_nuclear(
            self.params,
            momenta,
            np.asarray(self.storage.iM[self.traj.id, idx][0]),
            accepted,
            initial,
            np.asarray(self.storage.ham_adi[self.traj.id, idx]),
            energies_dia=np.asarray(self.storage.ham_dia[self.traj.id, idx]),
            dc1_adi=None if self.storage.dc1_adi is None else np.asarray(self.storage.dc1_adi[self.traj.id, idx]),
            d1ham_adi=None if self.storage.d1ham_adi is None else np.asarray(self.storage.d1ham_adi[self.traj.id, idx]),
        )
        self.storage.p[self.traj.id, idx] = momenta

    def _finalize_hop_momenta(self):
        """Apply hop momentum handling after the new-surface half-kick.

        The hopping decision and active-state change occur at the nuclear
        position at ``t + dt`` so the new state is used for the second force
        half-kick. Applying the energy-conserving rescaling to the resulting
        full-step momentum prevents that half-kick from undoing conservation
        established at the hop. Trajectories without a pending TSH decision
        are unchanged.

        This differs deliberately from the ordering in C++ ``Dynamics.cpp``,
        where rescaling precedes the final half-kick. That ordering produces a
        finite energy jump at hops even when rescaling option 200 is exact.
        """
        pending = getattr(self, "_pending_hop_momentum_update", None)
        if pending is None:
            return
        accepted, initial = pending
        self._rescale_hop_momenta(accepted, initial)
        self._correct_hop_energy_defect(accepted, initial)
        self._pending_hop_momentum_update = None

    def _trajectory_total_energies(self):
        """Return active-surface kinetic plus potential energy per trajectory."""
        idx = self.traj.tbf_ids
        p = np.asarray(self.storage.p[self.traj.id, idx])
        inv_mass = np.asarray(self.storage.iM[self.traj.id, idx])
        kinetic = 0.5 * np.sum(p * p * inv_mass, axis=-1)
        states = np.asarray(self.storage.act_states[self.traj.id, idx], dtype=int)
        energies = np.diagonal(
            np.asarray(self.storage.ham_adi[self.traj.id, idx]).real,
            axis1=-2, axis2=-1,
        )
        return kinetic + energies[np.arange(len(states)), states]

    def _correct_hop_energy_defect(self, accepted, initial):
        """Remove the split-force energy defect on successful hops.

        Velocity Verlet brackets the state change with an old-surface and a
        new-surface force half-kick. Even exact hop rescaling therefore leaves
        a small integration defect. For energy-conserving rescaling options,
        this function removes that defect along the same physical direction:
        derivative coupling for 200-series options, force difference for
        210-series options, and momentum for uniform 100/110-series scaling.
        Non-hopping trajectories and option 0 are untouched.
        """
        algorithm = int(_param(self.params, "momenta_rescaling_algo", 0))
        if algorithm not in (100, 101, 110, 111, 200, 201, 210, 211):
            return
        idx = self.traj.tbf_ids
        p = np.array(self.storage.p[self.traj.id, idx], copy=True)
        inv_mass = np.asarray(self.storage.iM[self.traj.id, idx])
        current = self._trajectory_total_energies()
        target = np.asarray(self._step_energy_before_hops)
        dc = None if self.storage.dc1_adi is None else np.asarray(
            self.storage.dc1_adi[self.traj.id, idx]
        )
        deriv = None if self.storage.d1ham_adi is None else np.asarray(
            self.storage.d1ham_adi[self.traj.id, idx]
        )
        for t, (old, new) in enumerate(zip(initial, accepted)):
            if old == new or abs(current[t] - target[t]) <= 1.0e-14:
                continue
            if algorithm in (200, 201):
                direction = np.real(dc[t, :, old, new])
            elif algorithm in (210, 211):
                diagonal = np.diagonal(deriv[t].real, axis1=-2, axis2=-1)
                direction = diagonal[:, old] - diagonal[:, new]
            else:
                direction = p[t].copy()
            if np.linalg.norm(direction) <= 1.0e-14:
                direction = p[t].copy()
            rescale_along_vector(
                target[t], current[t], p[t], inv_mass[t], direction,
                which_dofs=_param(self.params, "quantum_dofs", None),
            )
        self.storage.p[self.traj.id, idx] = p

    def _map_active_states(self, rep_sh, states):
        idx = self.traj.tbf_ids
        transform = np.asarray(self.storage.basis_transform[self.traj.id, idx])
        if rep_sh == "adiabatic":
            mapped = np.argmax(np.abs(transform[np.arange(len(states)), :, states]), axis=1)
            self.storage.act_states_dia[self.traj.id, idx] = mapped
        else:
            mapped = np.argmax(np.abs(transform[np.arange(len(states)), states, :]), axis=1)
            self.storage.act_states[self.traj.id, idx] = mapped

    def _observable_config(self):
        properties = tuple(_param(self.params, "properties_to_save", []) or ())
        return ObservableConfig(
            rep=self._rep_from_params("sh"),
            keywords=properties or None,
            potential="ehrenfest" if self.method == "ehrenfest" else "active",
        )

    def _setup_saver(self):
        if self.saver is not None:
            return
        hdf5_level = int(_param(self.params, "hdf5_output_level", -1))
        other_level = max(
            int(_param(self.params, "mem_output_level", -1)),
            int(_param(self.params, "txt_output_level", -1)),
        )
        if hdf5_level < 0 and other_level < 0:
            return
        if hdf5_level >= 0:
            from .savers.disk import HDF5Saver
            self.saver = HDF5Saver(
                output_dir=_param(self.params, "prefix", "out"),
                compression="gzip" if int(_param(self.params, "use_compression", 0)) else None,
            )
        else:
            from .savers.disk import FaultTolerantSaver
            self.saver = FaultTolerantSaver(
                output_dir=_param(self.params, "prefix", "out"),
                compressed=bool(_param(self.params, "use_compression", 0)),
            )

    def _save_if_due(self, force=False):
        if self.saver is None:
            return
        stride = int(_param(self.params, "nprint", 1))
        if force or self.storage.timestep % stride == 0:
            self.saver.save_observables(
                self.storage,
                self.traj,
                self.observable_config,
                step=self.storage.timestep,
                time=self.time,
                metadata=self.metadata,
            )

    def close(self):
        """Close an attached saver when it owns external resources."""

        if self.saver is not None and hasattr(self.saver, "close"):
            self.saver.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def _validate(self):
        if self.method not in ("adiabatic", "ehrenfest", "tsh"):
            raise ValueError("method must be 'adiabatic', 'ehrenfest', or 'tsh'")
        if self.rep not in ("adiabatic", "diabatic"):
            raise ValueError("rep must be 'adiabatic' or 'diabatic'")
        if self.force_mode not in ("state-specific", "ehrenfest", "none"):
            raise ValueError("force_mode must be 'state-specific', 'ehrenfest', or 'none'")
        if self.method == "adiabatic" and self.force_mode == "ehrenfest":
            raise ValueError("adiabatic dynamics cannot use Ehrenfest forces")
        if self.force_mode == "state-specific" and self.rep != "adiabatic":
            raise NotImplementedError("state-specific forces currently require adiabatic rep")
        self._validate_decoherence_configuration()

    def _validate_decoherence_configuration(self):
        """Validate the coupled TSH/decoherence choices used by C++.

        This keeps representation and parameter constraints at the workflow
        boundary rather than scattering them through numerical routines.
        """
        algo = int(_param(self.params, "decoherence_algo", -1))
        times = int(_param(self.params, "decoherence_times_type", -1))
        tsh = int(_param(self.params, "tsh_method", -1))
        if algo not in range(-1, 10):
            raise ValueError("decoherence_algo must be one of -1 through 9")
        if times not in range(-1, 6):
            raise ValueError("decoherence_times_type must be one of -1 through 5")
        if tsh == 5 and algo != -1:
            raise ValueError("legacy DISH (tsh_method=5) requires decoherence_algo=-1")
        if algo in (1, 2, 8) and self.method != "tsh":
            raise ValueError(f"decoherence_algo={algo} requires a TSH hop outcome")
        if algo in (0, 1, 2, 3, 4, 5, 6, 7, 8, 9) and self.rep != "adiabatic":
            raise ValueError("the selected C++ decoherence algorithm requires rep_tdse=1")
        if times == 0 and _param(self.params, "decoherence_rates", None) is None:
            raise ValueError("decoherence_times_type=0 requires decoherence_rates")
        if times in (2, 3) and _param(self.params, "schwartz_decoherence_inv_alpha", None) is None:
            raise ValueError(f"decoherence_times_type={times} requires schwartz_decoherence_inv_alpha")
        if times == 4 and _param(self.params, "schwartz_interaction_width", None) is None:
            raise ValueError("decoherence_times_type=4 requires schwartz_interaction_width")
        if int(_param(self.params, "dephasing_informed", 0)) and _param(self.params, "ave_gaps", None) is None:
            raise ValueError("dephasing_informed=1 requires ave_gaps")


def run_dynamics(*args, **kwargs):
    """Convenience constructor for a storage-backed dynamics engine."""

    return DynamicsEngine(*args, **kwargs)


def dynamics_defaults(overrides: dict[str, Any] | None = None) -> dict[str, Any]:
    """Return legacy-compatible defaults for the unified dynamics workflow.

    The result is an ordinary dictionary suitable for older input-building
    code. Names follow ``libra_py.dynamics.tsh`` wherever the new storage model
    has an equivalent. Unsupported legacy keys may still be supplied because
    :class:`DynControlParams` retains permissive ``set_parameters`` behavior.
    """

    defaults = asdict(DynControlParams())
    defaults["properties_to_save"] = [
        "timestep",
        "time",
        "Ekin_ave",
        "Epot_ave",
        "Etot_ave",
        "states",
        "se_pop_adi",
        "se_pop_dia",
        "sh_pop_adi",
        "q",
        "p",
        "f",
        "Cadi",
        "Cdia",
        "hvib_adi",
        "hvib_dia",
        "St",
        "basis_transform",
        "projector",
    ]
    if overrides:
        defaults.update(overrides)
    return defaults


def _param(params, name: str, default=None):
    if isinstance(params, dict):
        return params.get(name, default)
    return getattr(params, name, default)
