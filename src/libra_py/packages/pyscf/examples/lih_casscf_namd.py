"""
Non-NBRA LiH NAMD example built on the formal PySCF interface + NAMD adapter.

This example keeps the setup under ``pyscf/examples`` but uses the newer
architecture:

* ``CASSCF`` implements ``ElectronicStructureStrategy``
* ``NAMDRunner`` from ``libra_py.namd_adapter`` bridges the strategy into Libra

All inputs are defined directly in Python when instantiating the relevant
objects and parameter dictionaries. This keeps the example close to the
``tutorial.ipynb`` style while using the newer adapter stack.

Examples
--------

List supported aliases:

    python src/libra_py/packages/pyscf/examples/lih_casscf_namd.py --list-methods

Dry-run the default setup without launching dynamics:

    python src/libra_py/packages/pyscf/examples/lih_casscf_namd.py --dry-run

Run FSSH2:

    python src/libra_py/packages/pyscf/examples/lih_casscf_namd.py --method fssh2

Run DISH + EDC + GFSH:

    python src/libra_py/packages/pyscf/examples/lih_casscf_namd.py --method dish_edc_gfsh
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Any

if __name__ == "__main__" and __package__ is None:
    file_path = Path(__file__).resolve()
    for parent in file_path.parents:
        if parent.name == "src":
            sys.path.insert(0, str(parent))
            break
    else:
        raise RuntimeError("Could not locate src/ directory on path for libra_py import")

from libra_py.packages.pyscf.implementations.casscf import CASSCF
from libra_py.packages.pyscf.interfaces import ElectronicStructureStrategy
from libra_py.namd_adapter import NAMDRunner


ATOM_LABELS = ["Li", "H"]

METHOD_RECIPES = {
    "fssh": "fssh_v_plus",
    "fssh2": "fssh2_v_plus",
    "gfsh": "gfsh_v_plus",
    "dish_edc_fssh": "dish_edc_fssh_v_plus",
    "dish_edc_fssh2": "dish_edc_fssh2_v_plus",
    "dish_edc_gfsh": "dish_edc_gfsh_v_plus",
}
DEFAULT_RECIPE = "fssh2_v_plus"

ES_PARAMS = {
    "norbcas": 5,
    "nelecas": 2,
    "nroots": 4,
    "basis": "sto-3g",
    "charge": 0,
}

DYN_PARAMS = {
    "ntraj": 4,
    "nsteps": 5,
    "dt_fs": 1.0,
    "compute_gradients": True,
    "use_nac_vectors": False,
    "seed": 0,
    "model_params": {
        "model": 0,
        "model0": 0,
    },
    "dyn_params_override": {
        "rep_tdse": 1,
        "rep_sh": 1,
        "rep_force": 1,
        "num_electronic_substeps": 1,
        "isNBRA": 0,
        "is_nbra": 0,
        "hdf5_output_level": -1,
        "txt_output_level": -1,
        "txt2_output_level": -1,
        "mem_output_level": 3,
    },
}

NUCLEAR_INIT = {
    "ndof": 6,
    "q": [
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        3.02356158,
    ],
    "p": [
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    ],
    "mass": [
        12760.2210954,
        12760.2210954,
        12760.2210954,
        1837.15267389,
        1837.15267389,
        1837.15267389,
    ],
    "force_constant": [
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    ],
    "q_width": [
        0.25,
        0.25,
        0.25,
        0.25,
        0.25,
        0.25,
    ],
    "p_width": [
        0.1,
        0.1,
        0.1,
        0.1,
        0.1,
        0.1,
    ],
    "init_type": 4,
}

ELECTRONIC_INIT = {
    "verbosity": 2,
    "init_dm_type": 0,
    "ndia": 4,
    "nadi": 4,
    "rep": 1,
    "init_type": 1,
    "istate": 2,
    "istates": [
        0.0,
        0.0,
        1.0,
        0.0,
    ],
}


def sanitize_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]+", "_", value).strip("_")


def build_strategy_factory():
    params = dict(ES_PARAMS)

    def factory() -> ElectronicStructureStrategy:
        return CASSCF(**params)

    return factory


def choose_recipe(
    method: str | None,
    recipe: str | None,
) -> tuple[str, str]:
    if recipe:
        return recipe, sanitize_name(recipe)

    if method:
        if method not in METHOD_RECIPES:
            supported = ", ".join(sorted(METHOD_RECIPES))
            raise ValueError(
                f"Unknown method alias '{method}'. Supported aliases: {supported}"
            )
        return METHOD_RECIPES[method], method

    default_alias = next(
        (name for name, value in METHOD_RECIPES.items() if value == DEFAULT_RECIPE),
        sanitize_name(DEFAULT_RECIPE),
    )
    return DEFAULT_RECIPE, default_alias


def run_example(
    method: str | None = None,
    recipe: str | None = None,
    dry_run: bool = False,
) -> Any:
    chosen_recipe, method_tag = choose_recipe(method, recipe)

    dyn_overrides = dict(DYN_PARAMS["dyn_params_override"])
    dyn_overrides.setdefault("prefix", f"lih_casscf_{method_tag}")
    dyn_overrides.setdefault("prefix2", f"lih_casscf_{method_tag}_txt2")

    if dry_run:
        print("Dry run")
        print(f"recipe={chosen_recipe}")
        print(
            "strategy="
            f"CASSCF(norbcas={ES_PARAMS['norbcas']}, "
            f"nelecas={ES_PARAMS['nelecas']}, "
            f"nroots={ES_PARAMS['nroots']}, "
            f"basis='{ES_PARAMS['basis']}', "
            f"charge={ES_PARAMS.get('charge', 0)})"
        )
        print(
            f"atom_labels={ATOM_LABELS} "
            f"ntraj={DYN_PARAMS['ntraj']} "
            f"nsteps={DYN_PARAMS['nsteps']} "
            f"dt_fs={DYN_PARAMS['dt_fs']}"
        )
        print(f"output_prefix={dyn_overrides['prefix']}")
        return None

    runner = NAMDRunner(
        strategy_factory=build_strategy_factory(),
        atom_labels=ATOM_LABELS,
        compute_gradients=DYN_PARAMS.get("compute_gradients", True),
        use_nac_vectors=DYN_PARAMS.get("use_nac_vectors", False),
    )

    print(f"Running LiH CAS(2e,5o) CASSCF/STO-3G with recipe '{chosen_recipe}'")
    return runner.run(
        recipe=chosen_recipe,
        init_nucl=NUCLEAR_INIT,
        init_elec=ELECTRONIC_INIT,
        nsteps=DYN_PARAMS.get("nsteps", 1000),
        dt=DYN_PARAMS.get("dt_fs", 1.0),
        ntraj=DYN_PARAMS.get("ntraj", 1),
        nstates=ES_PARAMS["nroots"],
        model_params=DYN_PARAMS.get("model_params"),
        dyn_params_override=dyn_overrides,
        seed=DYN_PARAMS.get("seed", 0),
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--method",
        default=None,
        help="Method alias from the in-script METHOD_RECIPES mapping.",
    )
    parser.add_argument(
        "--recipe",
        default=None,
        help="Raw Libra recipe name. Overrides --method.",
    )
    parser.add_argument(
        "--list-methods",
        action="store_true",
        help="Print method aliases and exit.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the resolved setup without launching Libra dynamics.",
    )
    return parser


def main(argv: list[str] | None = None) -> Any:
    args = build_parser().parse_args(argv)

    if args.list_methods:
        for alias, recipe in sorted(METHOD_RECIPES.items()):
            print(f"{alias}: {recipe}")
        return None

    return run_example(
        method=args.method,
        recipe=args.recipe,
        dry_run=args.dry_run,
    )


if __name__ == "__main__":
    main()
