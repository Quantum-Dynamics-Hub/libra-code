from __future__ import annotations

from importlib import import_module
from types import SimpleNamespace

import numpy as np

from liblibra_core import CMATRIX, CMATRIXList, Random
from libra_py import units
from libra_py.dynamics.tsh import compute as tsh_dynamics

from libra_py.packages.pyscf.implementations.casscf import CASSCF
from libra_py.packages.pyscf.interfaces import MolecularGeometry


# ============================================================
# 0. User choices
# ============================================================

ATOM_LABELS = ["H", "H"]
NSTATES = 2
NTRAJ = 100
DT = 41.0  # a.u.

INITIAL_STATE = 1
NAM_METHOD = "fssh"
OUTPUT_PREFIX = "h2_nbra_preset_bond"


# ============================================================
# 1. Preset nuclear trajectory: a time series of bond lengths
# ============================================================

def make_bond_length_series() -> list[float]:
    """Return the prescribed H-H bond lengths in Angstrom."""
    nsteps = 40
    center = 0.80
    amplitude = 0.08
    period = 20.0
    return [
        center + amplitude * np.sin(2.0 * np.pi * step / period)
        for step in range(nsteps + 1)
    ]


def make_diatomic_geometry(bond_length_angstrom: float) -> MolecularGeometry:
    half = 0.5 * bond_length_angstrom
    return MolecularGeometry(
        atom_labels=ATOM_LABELS,
        coords_angstrom=np.array(
            [
                [0.0, 0.0, -half],
                [0.0, 0.0, half],
            ],
            dtype=float,
        ),
    )


# ============================================================
# 2. Generate the NBRA electronic data along the preset path
# ============================================================

def identity_cmatrix(nstates: int) -> CMATRIX:
    mat = CMATRIX(nstates, nstates)
    for i in range(nstates):
        mat.set(i, i, 1.0 + 0.0j)
    return mat


def zero_cmatrix_list(nitems: int, nstates: int) -> CMATRIXList:
    mats = CMATRIXList()
    for _ in range(nitems):
        mats.append(CMATRIX(nstates, nstates))
    return mats


def build_hvib_from_time_overlap(energies: list[float], st: np.ndarray) -> CMATRIX:
    hvib = CMATRIX(NSTATES, NSTATES)
    for i, energy in enumerate(energies):
        hvib.set(i, i, energy + 0.0j)

    for i in range(NSTATES):
        for j in range(i + 1, NSTATES):
            dij = (float(st[i, j]) - float(st[j, i])) / (2.0 * DT)
            hvib.set(i, j, -1j * dij)
            hvib.set(j, i, +1j * dij)
    return hvib


def build_ham_from_energies(energies: list[float]) -> CMATRIX:
    ham = CMATRIX(NSTATES, NSTATES)
    for i, energy in enumerate(energies):
        ham.set(i, i, energy + 0.0j)
    return ham


def precompute_hvib_series(bond_lengths_angstrom: list[float]) -> list[SimpleNamespace]:
    """Run PySCF along the preset bond-length path and cache Libra-ready data."""
    engine = CASSCF(
        norbcas=2,
        nelecas=2,
        nroots=NSTATES,
        basis="sto-3g",
        charge=0,
        ntraj=1,
    )

    series = []
    for step, bond_length in enumerate(bond_lengths_angstrom):
        geom = make_diatomic_geometry(bond_length)
        engine.set_geom_and_run_hf(geom, traj_id=0)

        energies = [
            engine.compute_energy(root, traj_id=0)
            for root in range(NSTATES)
        ]

        if step == 0:
            st = np.eye(NSTATES)
        else:
            st = engine.time_overlap_matrix(NSTATES, traj_id=0)

        item = SimpleNamespace()
        item.ham_adi = build_ham_from_energies(energies)
        item.hvib_adi = build_hvib_from_time_overlap(energies, st)
        item.nac_adi = CMATRIX(NSTATES, NSTATES)
        item.basis_transform = identity_cmatrix(NSTATES)
        item.time_overlap_adi = CMATRIX(NSTATES, NSTATES)
        item.d1ham_adi = zero_cmatrix_list(1, NSTATES)

        for i in range(NSTATES):
            for j in range(NSTATES):
                item.time_overlap_adi.set(i, j, float(st[i, j]) + 0.0j)

        series.append(item)

    return series


# ============================================================
# 3. NBRA compute_model callback: read precomputed Hvib(t)
# ============================================================

def nbra_compute_adi(q, params, full_id):
    """Libra callback for NBRA: ignore q and return precomputed step data."""
    series = params["hvib_series"]
    timestep = int(params.get("timestep", 0)) % len(series)
    data = series[timestep]

    obj = SimpleNamespace()
    obj.ham_adi = CMATRIX(data.ham_adi)
    obj.hvib_adi = CMATRIX(data.hvib_adi)
    obj.nac_adi = CMATRIX(data.nac_adi)
    obj.basis_transform = CMATRIX(data.basis_transform)
    obj.time_overlap_adi = CMATRIX(data.time_overlap_adi)
    obj.d1ham_adi = CMATRIXList()
    for mat in data.d1ham_adi:
        obj.d1ham_adi.append(CMATRIX(mat))
    return obj


# ============================================================
# 4. Libra NBRA input
# ============================================================

def load_recipe(dyn_params, method_name: str) -> None:
    recipe_map = {
        "fssh": "fssh_v_plus",
        "fssh2": "fssh2_v_plus",
        "gfsh": "gfsh_v_plus",
        "dish_edc_gfsh": "dish_edc_gfsh_v_plus",
    }
    recipe_name = recipe_map[method_name]
    recipe = import_module(f"libra_py.dynamics.tsh.recipes.{recipe_name}")
    recipe.load(dyn_params)


def make_init_nucl(first_bond_length_angstrom: float) -> dict:
    """NBRA needs a nuclear container, but q/p are not the source of Hvib."""
    return {
        "ndof": 1,
        "q": [first_bond_length_angstrom * units.Angst],
        "p": [0.0],
        "mass": [1.0],
        "init_type": 0,
    }


def make_init_elec() -> dict:
    istates = [0.0 for _ in range(NSTATES)]
    istates[INITIAL_STATE] = 1.0
    return {
        "ndia": NSTATES,
        "nadi": NSTATES,
        "rep": 1,
        "init_type": 1,
        "istate": INITIAL_STATE,
        "istates": istates,
    }


def make_dyn_params(nsteps: int, nfiles: int) -> dict:
    dyn_params = {
        "ntraj": NTRAJ,
        "nsteps": nsteps,
        "nstates": NSTATES,
        "dt": DT,
        "isNBRA": 1,
        "is_nbra": 1,
        "icond": 0,
        "nfiles": nfiles,
        "which_adi_states": list(range(NSTATES)),
        "which_dia_states": list(range(NSTATES)),
        "num_electronic_substeps": 1,
        "mem_output_level": 3,
        "prefix": OUTPUT_PREFIX,
        "prefix2": f"{OUTPUT_PREFIX}_txt2",
    }
    load_recipe(dyn_params, NAM_METHOD)

    # NBRA consumes precomputed adiabatic Hvib data. Do not recompute ES data,
    # forces, NACs, or Hvib from propagated nuclear coordinates.
    dyn_params.update(
        {
            "ham_update_method": 2,
            "ham_transform_method": 0,
            "time_overlap_method": 0,
            "nac_update_method": 0,
            "hvib_update_method": 0,
            "force_method": 0,
            "hop_acceptance_algo": 0,
            "momenta_rescaling_algo": 0,
            "ensemble": 0,
        }
    )
    return dyn_params


def main():
    bond_lengths = make_bond_length_series()
    hvib_series = precompute_hvib_series(bond_lengths)

    init_nucl = make_init_nucl(bond_lengths[0])
    init_elec = make_init_elec()
    dyn_params = make_dyn_params(
        nsteps=len(hvib_series) - 1,
        nfiles=len(hvib_series),
    )
    model_params = {
        "model": 0,
        "model0": 0,
        "nstates": NSTATES,
        "hvib_series": hvib_series,
    }

    rnd = Random()
    return tsh_dynamics.generic_recipe(
        dyn_params,
        nbra_compute_adi,
        model_params,
        init_elec,
        init_nucl,
        rnd,
    )


if __name__ == "__main__":
    main()
