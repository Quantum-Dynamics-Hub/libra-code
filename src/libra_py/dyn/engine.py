from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Callable, Literal

import numpy as np

from .control_params import DynControlParams
from .hamiltonians import HamiltonianEngine
from .hopping import (
    accept_hops,
    handle_hops_nuclear,
    hop_proposal_probabilities,
    propose_hops,
)
from .observables import ObservableConfig
from .propagation.coupled import (
    ehrenfest_forces,
    state_specific_forces,
    update_density,
)
from .propagation.electronic import exp_propagator, split_step_propagator, tdse_step
from .propagation.integrators import normalize_dt, run_steps
from .propagation.nuclear import drift, kick
from .transformations.basis_rotation import (
    storage_amplitudes_adi_to_dia,
    storage_amplitudes_dia_to_adi,
)
from .transformations.local_diabatization import orthogonalized_T

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
        kick(self.storage, self.traj, 0.5 * dt)
        drift(self.storage, self.traj, dt)

        self._save_electronic_history()
        self._previous_hamiltonian = self._hamiltonian_snapshot()
        self.evaluate_hamiltonian()
        ld_transform = self._state_tracking_transform()
        if self._uses_tdse:
            self.propagate_electronic(dt, T=ld_transform)
            self._apply_state_tracking_to_active_states(ld_transform)
            self._synchronize_representations()
            update_density(self.storage, self.traj, self.rep)

        hopping = None
        proposed = None
        accepted = None
        if self.method == "tsh":
            hopping, proposed, accepted = self.surface_hopping_step()
        forces = self.compute_forces()
        kick(self.storage, self.traj, 0.5 * dt, forces)

        self.evaluate_hamiltonian()
        if self._uses_tdse:
            update_density(self.storage, self.traj, self.rep)

        self.time += dt
        self.storage.timestep += 1
        self._save_if_due()
        return self._result(forces, hopping, proposed, accepted)

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
        result = None
        # Options 3 and 4 follow the C++ one- and two-point adiabatic-Hvib
        # schemes.  Option 4 uses a symmetric old/new split, avoiding the
        # systematic right-endpoint bias of propagating every option with the
        # newly evaluated Hamiltonian only.  The LD families (0--2, 10--12)
        # retain the projection-matrix path below.
        propagator = self.propagator
        previous = None
        if self.rep == "adiabatic" and integrator in (3, 4):
            previous = self._previous_hamiltonian
            if integrator == 4:
                propagator = split_step_propagator
        for substep in range(nsubsteps):
            substep_dt = dt / nsubsteps
            state_or_dt = previous if integrator == 3 else substep_dt
            result = tdse_step(
                self.traj,
                self.storage,
                state_or_dt,
                substep_dt if integrator == 3 else None,
                backend=self.storage.backend,
                propagator=propagator,
                rep=self.rep,
                hamiltonian_type=self.hamiltonian_type,
                T=T if substep == 0 else None,
                previous_state=previous,
            )
        return result

    def surface_hopping_step(self):
        """Propose, accept, and handle hops for active TBFs."""

        idx = self.traj.tbf_ids
        rep_sh = self._rep_from_params("sh")
        density = self._storage_slice("dm", rep_sh)
        hvib = self._storage_slice("hvib", rep_sh)
        states_field = "act_states" if rep_sh == "adiabatic" else "act_states_dia"
        initial = np.array(getattr(self.storage, states_field)[self.traj.id, idx], copy=True)
        previous_density = self._previous_density(rep_sh)
        method = int(_param(self.params, "tsh_method", -1))

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
        self._rescale_hop_momenta(accepted, initial)
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

    def _result(self, forces, hopping=None, proposed=None, accepted=None):
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
        if self.rep != "adiabatic" or int(_param(self.params, "electronic_integrator", 0)) not in (0, 1, 2, 10, 11, 12):
            return None
        if int(_param(self.params, "assume_always_consistent", 0)):
            return None
        overlap = np.asarray(self.storage.time_overlap_adi[self.traj.id, self.traj.tbf_ids])
        if not np.any(overlap):
            return None
        transforms = []
        for matrix in overlap:
            transforms.append(
                orthogonalized_T(
                    np.linalg.pinv(matrix),
                    tol=float(_param(self.params, "phase_correction_tol", 1e-3)),
                )
            )
        transform = np.asarray(transforms)
        self.storage.proj_adi[self.traj.id, self.traj.tbf_ids] = transform
        return transform

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
            "d1ham_dia", "d1ham_adi", "dc1_adi",
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
