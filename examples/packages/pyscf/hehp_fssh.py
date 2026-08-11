"""
Minimal on-the-fly Libra FSSH example for HeH+ with a PySCF CASSCF backend.

This example follows the same CASSCF instantiation pattern as
``implementations/test_casscf.py``:

    CASSCF(norbcas=2, nelecas=2, nroots=NSTATES, basis="sto-3g", charge=1)

The important design point for a trajectory swarm is that ``pyscf_compute_adi``
creates and keeps one stateful ``CASSCF`` object per trajectory. The backend
caches the previous SCF/CASSCF state to build time overlaps, so sharing one
instance across trajectories would mix their electronic histories.
"""

from __future__ import annotations

from typing import Any

from liblibra_core import (
    Random,
    dyn_variables,
    nHamiltonian,
    update_Hamiltonian_variables,
)

import libra_py.units as units
from libra_py.dynamics.tsh.compute import run_dynamics
from libra_py.dynamics.tsh.recipes.fssh_h_plus import load as load_fssh_recipe
from libra_py.packages.pyscf.implementations.casscf import CASSCF
from libra_py.packages.pyscf.methods import pyscf_compute_adi


NTRAJ = 2
NSTATES = 3
NSTEPS = 2
DT_FS = 0.5
INITIAL_BOND_ANG = 0.7746
OUTPUT_PREFIX = "hehp_sa3_casscf_sto3g_fssh"

HE_MASS_AU = 4.002602 * units.amu
H_MASS_AU = 1.007825 * units.amu


def build_strategy_factory(nstates: int):
    def factory() -> CASSCF:
        return CASSCF(
            norbcas=2,
            nelecas=2,
            nroots=nstates,
            basis="sto-3g",
            charge=1,
            unit="Bohr",
        )

    return factory


def make_dyn_params(
    ntraj: int,
    nsteps: int,
    dt_au: float,
    prefix: str,
) -> dict[str, Any]:
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
            # The callback writes time_overlap_adi, so do not overwrite it
            # from basis_transform in the C++ dynamics layer.
            "time_overlap_method": 0,
            "nac_update_method": 2,
            "nac_algo": 0,
            "hvib_update_method": 1,
            "hop_acceptance_algo": 10,
            "momenta_rescaling_algo": 100,
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
    gradient: str = "active",
    initial_state: int = 1,
) -> Any:
    if gradient not in {"active", "all", "none"}:
        raise ValueError("gradient must be one of 'active', 'all', or 'none'.")
    if not 0 <= initial_state < nstates:
        raise ValueError(f"initial_state must be in [0, {nstates}).")

    atom_labels = ["He", "H"]
    ndof = 3 * len(atom_labels)
    dt_au = dt_fs * units.fs2au

    compute_model = pyscf_compute_adi
    model_params = {
        "model": 0,
        "model0": 0,
        "nstates": nstates,
        "atom_labels": atom_labels,
        "strategy_factory": build_strategy_factory(nstates),
        "gradient": gradient,
        "time_overlap": True,
        "dt": dt_au,
        "act_state": {itraj: initial_state for itraj in range(ntraj)},
    }

    dyn_params = make_dyn_params(
        ntraj=ntraj,
        nsteps=nsteps,
        dt_au=dt_au,
        prefix=prefix,
    )
    init_nucl = make_init_nucl(initial_bond_ang)
    init_elec = make_init_elec(nstates, ntraj, initial_state=initial_state)

    rnd = Random()
    dyn_var = dyn_variables(nstates, nstates, ndof, ntraj)
    dyn_var.init_nuclear_dyn_var(init_nucl, rnd)
    dyn_var.init_amplitudes(init_elec, rnd)
    dyn_var.init_density_matrix(init_elec)
    dyn_var.init_auxiliary_variables(init_elec, rnd)

    ham = nHamiltonian(nstates, nstates, ndof)
    ham.add_new_children(nstates, nstates, ndof, ntraj)
    ham.init_all(2, 1)

    warmup_params = dict(model_params)
    warmup_params["timestep"] = 0
    update_Hamiltonian_variables(
        dyn_params,
        dyn_var,
        ham,
        ham,
        compute_model,
        warmup_params,
        0,
    )
    update_Hamiltonian_variables(
        dyn_params,
        dyn_var,
        ham,
        ham,
        compute_model,
        warmup_params,
        1,
    )

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
