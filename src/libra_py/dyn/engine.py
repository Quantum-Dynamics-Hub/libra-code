from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Literal

import numpy as np

from .control_params import DynControlParams
from .hamiltonians import HamiltonianEngine
from .propagation.coupled import (
    ehrenfest_forces,
    state_specific_forces,
    update_density,
)
from .propagation.electronic import exp_propagator, tdse_step
from .propagation.integrators import normalize_dt, run_steps
from .propagation.nuclear import drift, kick

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


@dataclass
class DynamicsEngine:
    """
    Storage-backed nonadiabatic dynamics driver.

    Supported first-pass methods:
    - ``adiabatic``: nuclei follow the active adiabatic state; no TDSE update.
    - ``ehrenfest``: TDSE plus Ehrenfest mean-field nuclear forces.
    - ``tsh``: TDSE on the active representation plus state-specific forces,
      but no hopping/decoherence/momentum rescaling yet.
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

    def __post_init__(self):
        self.params = self.params or DynControlParams()
        self.ham = self.ham_engine or HamiltonianEngine(self.storage.backend)
        self.rng = self.rng or np.random.default_rng()
        self.method = self.method or self._method_from_params()
        self.rep = self.rep or self._rep_from_params("tdse")
        self.force_mode = self.force_mode or self._force_mode_from_params()
        self._validate()

    def initialize(self, evaluate_hamiltonian: bool = True):
        """Build the initial Hamiltonian, forces, and density matrices."""

        if evaluate_hamiltonian:
            self.evaluate_hamiltonian()
        if self._uses_tdse:
            update_density(self.storage, self.traj, self.rep)
        self.compute_forces()
        return self

    def step(self, dt) -> StepResult:
        """Advance one mixed quantum-classical dynamics step."""

        dt = normalize_dt(dt)
        self.evaluate_hamiltonian()
        self.compute_forces()

        kick(self.storage, self.traj, 0.5 * dt)
        drift(self.storage, self.traj, dt)

        self.evaluate_hamiltonian()
        if self._uses_tdse:
            self.propagate_electronic(dt)
        forces = self.compute_forces()
        kick(self.storage, self.traj, 0.5 * dt, forces)

        self.evaluate_hamiltonian()
        if self._uses_tdse:
            update_density(self.storage, self.traj, self.rep)

        self.time += dt
        self.storage.timestep += 1
        return self._result(forces)

    def run(self, nsteps: int, dt):
        """Run several steps and return their summaries."""

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

    def propagate_electronic(self, dt):
        """Propagate active amplitudes with the TD-SE solver."""

        return tdse_step(
            self.traj,
            self.storage,
            dt,
            backend=self.storage.backend,
            propagator=self.propagator,
            rep=self.rep,
            hamiltonian_type=self.hamiltonian_type,
        )

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

    def _result(self, forces):
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
        name = "rep_tdse" if use == "tdse" else "rep_force"
        return "adiabatic" if int(_param(self.params, name, 1)) in (1, 3, 4) else "diabatic"

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


def _param(params, name: str, default=None):
    if isinstance(params, dict):
        return params.get(name, default)
    return getattr(params, name, default)
