"""
Minimal on-the-fly Libra FSSH example for HeH+ with a PySCF CASSCF backend.

This example follows the same CASSCF instantiation pattern as
``implementations/test_casscf.py``:

    CASSCF(norbcas=2, nelecas=2, nroots=NSTATES, basis="sto-3g", charge=1)

The important design point for a trajectory swarm is that each trajectory gets
its own stateful ``CASSCF`` object and its own adapter. The backend caches the
previous SCF/CASSCF state to build time overlaps, so sharing one instance
across trajectories would mix their electronic histories.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import numpy as np

# Allow running this file directly with proper package root in sys.path
#if __name__ == "__main__" and __package__ is None:
#    file_path = Path(__file__).resolve()
#    for parent in file_path.parents:
#        if parent.name == "src":
#            sys.path.insert(0, str(parent))
#            break
#    else:
#        raise RuntimeError("Could not locate src/ directory on path for libra_py import")

from liblibra_core import CMATRIX, Cpp2Py, Random, dyn_variables, nHamiltonian, update_Hamiltonian_variables

import libra_py.units as units
from libra_py.dynamics.tsh.compute import run_dynamics
from libra_py.dynamics.tsh.recipes.fssh_h_plus import load as load_fssh_recipe
from libra_py.packages.pyscf.adapter import LibraESAdapter
from libra_py.packages.pyscf.implementations.casscf import CASSCF


NTRAJ = 4
NSTATES = 2
NSTEPS = 25
DT_FS = 0.5
INITIAL_BOND_ANG = 0.7746
OUTPUT_PREFIX = "hehp_casscf_fssh"

HE_MASS_AU = 4.002602 * units.amu
H_MASS_AU = 1.007825 * units.amu


class _Result:
    """Small attribute container matching what Libra inspects."""


def _identity_cmatrix(size: int) -> CMATRIX:
    mat = CMATRIX(size, size)
    for i in range(size):
        mat.set(i, i, 1.0 + 0.0j)
    return mat


def _traj_index(full_id: Any) -> int:
    try:
        return int(Cpp2Py(full_id)[-1])
    except Exception:
        if hasattr(full_id, "__getitem__"):
            return int(full_id[-1])
        return 0


class FirstStepSafeAdapter(LibraESAdapter):
    """Adapter variant that returns an identity time overlap on the first call.

    The current ``LibraESAdapter.compute_model`` unconditionally asks the
    strategy for a time overlap matrix. For step 0 there is no previous frame
    yet, so the example handles that warmup call locally and then delegates all
    later calls back to the base adapter.
    """

    def __init__(self, strategy: CASSCF, atom_labels: list[str], nstates: int) -> None:
        super().__init__(strategy, atom_labels, nstates)
        self._has_previous_frame = False

    def compute_model(self, q: Any, params: dict[str, Any], full_id: Any) -> Any:
        if self._has_previous_frame:
            return super().compute_model(q, params, full_id)

        nstates = self._nstates
        natoms = self._natoms
        ndof = 3 * natoms

        geom = self._libra_q_to_geometry(q, _traj_index(full_id), self._atom_labels)
        self._strategy.set_geom_and_run_hf(geom)

        energies = np.array(
            [self._strategy.compute_energy(root) for root in range(nstates)]
        )
        grads = np.stack(
            [self._strategy.compute_gradient(root) for root in range(nstates)]
        )

        result = _Result()

        ham_adi = CMATRIX(nstates, nstates)
        for state in range(nstates):
            ham_adi.set(state, state, complex(energies[state], 0.0))
        result.ham_adi = ham_adi

        d1ham_adi = []
        for dof in range(ndof):
            atom_idx, xyz_idx = divmod(dof, 3)
            grad_block = CMATRIX(nstates, nstates)
            for state in range(nstates):
                grad_block.set(
                    state,
                    state,
                    complex(grads[state, atom_idx, xyz_idx], 0.0),
                )
            d1ham_adi.append(grad_block)
        result.d1ham_adi = d1ham_adi

        result.dc1_adi = [CMATRIX(nstates, nstates) for _ in range(ndof)]
        result.time_overlap_adi = _identity_cmatrix(nstates)
        result.ovlp_adi = _identity_cmatrix(nstates)

        self._has_previous_frame = True
        return result


def build_adapters(ntraj: int, atom_labels: list[str], nstates: int) -> list[FirstStepSafeAdapter]:
    adapters: list[FirstStepSafeAdapter] = []
    for _ in range(ntraj):
        strategy = CASSCF(
            norbcas=2,
            nelecas=2,
            nroots=nstates,
            basis="sto-3g",
            charge=1,
        )
        adapters.append(FirstStepSafeAdapter(strategy, atom_labels, nstates))
    return adapters


def build_compute_model(adapters: list[FirstStepSafeAdapter]):
    def compute_model(q: Any, params: dict[str, Any], full_id: Any) -> Any:
        return adapters[_traj_index(full_id)].compute_model(q, params, full_id)

    return compute_model


def make_dyn_params(ntraj: int, nsteps: int, dt_au: float, prefix: str) -> dict[str, Any]:
    dyn_params: dict[str, Any] = {}
    load_fssh_recipe(dyn_params)
    dyn_params.update(
        {
            "rep_tdse": 1,
            "rep_sh": 1,
            "rep_force": 1,
            "force_method": 1,
            "ham_update_method": 2,
            "ham_transform_method": 0,
            "time_overlap_method": 1,
            "nac_update_method": 2,
            "nac_algo": 0,
            "hvib_update_method": 1,
            # The adapter currently provides energies, gradients, and
            # time-overlap NAC information, but not derivative-coupling vectors.
            "hop_acceptance_algo": 21,
            "momenta_rescaling_algo": 210,
            "tsh_method": 0,
            "ntraj": ntraj,
            "isNBRA": 0,
            "is_nbra": 0,
            "dt": dt_au,
            "nsteps": nsteps,
            "nprint": 1,
            "prefix": prefix,
            "prefix2": f"{prefix}_txt2",
            "hdf5_output_level": -1,
            "mem_output_level": -1,
            "txt_output_level": 0,
            "txt2_output_level": -1,
        }
    )
    return dyn_params


def make_init_nucl(initial_bond_ang: float) -> dict[str, Any]:
    bond_bohr = initial_bond_ang * units.Angst

    q0 = [
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        bond_bohr,
    ]
    p0 = [
        0.0,
        0.0,
        -0.002,
        0.0,
        0.0,
        0.002,
    ]
    masses = [
        HE_MASS_AU,
        HE_MASS_AU,
        HE_MASS_AU,
        H_MASS_AU,
        H_MASS_AU,
        H_MASS_AU,
    ]

    return {
        "init_type": 0,
        "ndof": len(q0),
        "q": q0,
        "p": p0,
        "mass": masses,
    }


def make_init_elec(nstates: int, ntraj: int, initial_state: int = 1) -> dict[str, Any]:
    populations = [0.0] * nstates
    populations[initial_state] = 1.0
    return {
        "init_type": 3,
        "ndia": nstates,
        "nadi": nstates,
        "nstates": nstates,
        "istates": populations,
        "rep": 1,
        "ntraj": ntraj,
    }


def run_hehp_fssh(
    ntraj: int = NTRAJ,
    nstates: int = NSTATES,
    nsteps: int = NSTEPS,
    dt_fs: float = DT_FS,
    initial_bond_ang: float = INITIAL_BOND_ANG,
    prefix: str = OUTPUT_PREFIX,
) -> Any:
    atom_labels = ["He", "H"]
    ndof = 3 * len(atom_labels)

    adapters = build_adapters(ntraj, atom_labels, nstates)
    compute_model = build_compute_model(adapters)
    model_params = {"model": 0, "model0": 0, "nstates": nstates}

    dyn_params = make_dyn_params(
        ntraj=ntraj,
        nsteps=nsteps,
        dt_au=dt_fs * units.fs2au,
        prefix=prefix,
    )
    init_nucl = make_init_nucl(initial_bond_ang)
    init_elec = make_init_elec(nstates, ntraj)

    rnd = Random()
    dyn_var = dyn_variables(nstates, nstates, ndof, ntraj)
    dyn_var.init_nuclear_dyn_var(init_nucl, rnd)
    dyn_var.init_amplitudes(init_elec, rnd)
    dyn_var.init_density_matrix(init_elec)
    dyn_var.init_auxiliary_variables(init_elec, rnd)

    ham = nHamiltonian(nstates, nstates, ndof)
    ham.add_new_children(nstates, nstates, ndof, ntraj)
    ham.init_all(2, 1)

    model_params["timestep"] = 0
    update_Hamiltonian_variables(dyn_params, dyn_var, ham, ham, compute_model, model_params, 0)
    update_Hamiltonian_variables(dyn_params, dyn_var, ham, ham, compute_model, model_params, 1)

    dyn_var.update_basis_transform(ham)
    dyn_var.update_amplitudes({"rep_tdse": init_elec["rep"]}, ham)
    dyn_var.update_density_matrix(dyn_params, ham, 1)
    dyn_var.init_active_states(init_elec, rnd)

    return run_dynamics(dyn_var, dyn_params, ham, compute_model, model_params, rnd)


def main() -> None:
    run_hehp_fssh()
    print(
        f"Completed HeH+ CASSCF FSSH example. Outputs were written under {OUTPUT_PREFIX}/"
    )


if __name__ == "__main__":
    main()
