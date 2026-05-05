# *********************************************************************************
# * Copyright (C) 2026 Jieyang Gu <jieyanggu792@gmail.com>
# * This file is distributed under the terms of the GNU General Public License
# * as published by the Free Software Foundation, either version 3 of
# * the License, or (at your option) any later version.
# * See the file LICENSE in the root directory of this distribution
# * or <http://www.gnu.org/licenses/>.
# *
# *********************************************************************************/
"""
Config-driven entry point for the PySCF ↔ Libra adapter stack.

This module is deliberately thin:

* strategy selection and construction live in :mod:`pyscf.factory`
* Libra-specific callback packing lives in :mod:`libra_py.namd_adapter`
* concrete electronic-structure methods stay behind
  :class:`ElectronicStructureStrategy`

Usage
-----

Dry run without ``liblibra_core``::

    python src/libra_py/packages/pyscf/run.py --mode dry-run

Callback smoke test when Libra is installed::

    python src/libra_py/packages/pyscf/run.py --mode callback-demo

Run actual NAMD via ``generic_recipe``::

    python src/libra_py/packages/pyscf/run.py --mode namd --config my_cfg.json
"""

from __future__ import annotations

import argparse
import json
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

from libra_py.namd_adapter import _HAS_LIBRA
from libra_py.packages.pyscf.adapter import LibraESAdapter, NAMDRunner
from libra_py.packages.pyscf.factory import (
    build_strategy,
    build_strategy_factory,
    load_config,
)


_DEFAULT_CFG: dict[str, Any] = {
    "system": {
        "atom_labels": ["Li", "H"],
        "demo_geometries_angstrom": [
            [[0.0, 0.0, 0.0], [0.0, 0.0, 1.60]],
            [[0.0, 0.0, 0.0], [0.0, 0.0, 1.65]],
        ],
    },
    "electronic_structure": {
        "kind": "casscf",
        "params": {
            "norbcas": 2,
            "nelecas": 2,
            "nroots": 2,
            "basis": "sto-3g",
            "charge": 0,
        },
    },
    "dynamics": {
        "recipe": "fssh2_v_plus",
        "ntraj": 1,
        "nsteps": 5,
        "dt_fs": 0.5,
        "compute_gradients": True,
        "use_nac_vectors": False,
        "seed": 0,
        "model_params": {
            "model": 0,
            "model0": 0,
        },
        "dyn_params_override": {
            "prefix": "pyscf_libra_demo",
            "prefix2": "pyscf_libra_demo_txt2",
            "mem_output_level": -1,
            "txt_output_level": -1,
            "txt2_output_level": -1,
            "hdf5_output_level": -1,
        },
    },
    "initial_conditions": {
        "nuclear": {
            "init_type": 0,
            "ndof": 6,
            "q": [0.0, 0.0, 0.0, 0.0, 0.0, 3.02356158],
            "p": [0.0, 0.0, -0.002, 0.0, 0.0, 0.002],
            "mass": [
                12789.3935578,
                12789.3935578,
                12789.3935578,
                1837.15267389,
                1837.15267389,
                1837.15267389,
            ],
        },
        "electronic": {
            "init_type": 3,
            "ndia": 2,
            "nadi": 2,
            "nstates": 2,
            "istates": [0.0, 1.0],
            "rep": 1,
            "ntraj": 1,
        },
    },
}


def get_config(config_path: str | None) -> dict[str, Any]:
    """Load user config or fall back to the built-in example."""
    if config_path is None:
        return json.loads(json.dumps(_DEFAULT_CFG))
    return load_config(config_path)


def describe_strategy(cfg: dict[str, Any]) -> str:
    """Return a short human-readable summary of the configured strategy."""
    strategy = build_strategy(cfg["electronic_structure"])
    atom_labels = cfg["system"]["atom_labels"]
    return (
        f"strategy={strategy.__class__.__name__} "
        f"nstates={strategy.nstates} "
        f"natoms={len(atom_labels)} "
        f"has_nac_vectors={strategy.has_nac_vectors}"
    )


def dry_run(cfg: dict[str, Any]) -> None:
    """Instantiate the configured strategy/runner without touching Libra."""
    summary = describe_strategy(cfg)
    dynamics_cfg = cfg.get("dynamics", {})
    print("Dry run successful")
    print(summary)
    print(
        "recipe="
        f"{dynamics_cfg.get('recipe', 'fssh2_v_plus')} "
        f"ntraj={dynamics_cfg.get('ntraj', 1)} "
        f"compute_gradients={dynamics_cfg.get('compute_gradients', True)}"
    )


def _make_demo_q(atom_labels: list[str], coords_angstrom: list[list[float]]):
    """Construct a single-trajectory Libra MATRIX from Cartesian coordinates."""
    if not _HAS_LIBRA:
        raise RuntimeError("liblibra_core is required for callback-demo mode")

    from liblibra_core import MATRIX  # type: ignore[import-untyped]
    from libra_py import units

    q = MATRIX(3 * len(atom_labels), 1)
    flat_coords = []
    for xyz in coords_angstrom:
        if len(xyz) != 3:
            raise ValueError("Each coordinate must have exactly three components.")
        flat_coords.extend(float(value) * units.Angst for value in xyz)

    for idx, value in enumerate(flat_coords):
        q.set(idx, 0, value)
    return q


def callback_demo(cfg: dict[str, Any]) -> None:
    """Run two callback evaluations to verify first-step overlap handling."""
    if not _HAS_LIBRA:
        raise RuntimeError("liblibra_core is required for callback-demo mode")

    system_cfg = cfg["system"]
    atom_labels = system_cfg["atom_labels"]
    geometries = system_cfg.get("demo_geometries_angstrom")
    if geometries is None or len(geometries) < 2:
        raise ValueError(
            "callback-demo mode requires system.demo_geometries_angstrom "
            "with at least two geometries."
        )

    strategy = build_strategy(cfg["electronic_structure"])
    adapter = LibraESAdapter(strategy, atom_labels, nstates=strategy.nstates)

    first_q = _make_demo_q(atom_labels, geometries[0])
    second_q = _make_demo_q(atom_labels, geometries[1])

    first = adapter.compute_model(first_q, {}, [0, 0])
    second = adapter.compute_model(second_q, {}, [0, 0])

    print("Callback demo successful")
    print("first-call time-overlap diagonal:")
    for state in range(strategy.nstates):
        print(f"  state {state}: {first.time_overlap_adi.get(state, state).real:.10f}")
    print("second-call ham_adi diagonal (Hartree):")
    for state in range(strategy.nstates):
        print(f"  state {state}: {second.ham_adi.get(state, state).real:.10f}")


def run_namd(cfg: dict[str, Any]) -> Any:
    """Run config-driven NAMD using the shared :class:`NAMDRunner`."""
    if not _HAS_LIBRA:
        raise RuntimeError("liblibra_core is required for namd mode")

    system_cfg = cfg["system"]
    dynamics_cfg = cfg.get("dynamics", {})
    init_cfg = cfg.get("initial_conditions", {})

    runner = NAMDRunner(
        strategy_factory=build_strategy_factory(cfg["electronic_structure"]),
        atom_labels=system_cfg["atom_labels"],
        compute_gradients=dynamics_cfg.get("compute_gradients", True),
        use_nac_vectors=dynamics_cfg.get("use_nac_vectors", False),
    )

    result = runner.run(
        recipe=dynamics_cfg.get("recipe", "fssh2_v_plus"),
        init_nucl=init_cfg["nuclear"],
        init_elec=init_cfg["electronic"],
        nsteps=dynamics_cfg.get("nsteps", 1000),
        dt=dynamics_cfg.get("dt_fs", 1.0),
        ntraj=dynamics_cfg.get("ntraj", 1),
        nstates=dynamics_cfg.get("nstates"),
        model_params=dynamics_cfg.get("model_params"),
        dyn_params_override=dynamics_cfg.get("dyn_params_override"),
        seed=dynamics_cfg.get("seed", 0),
    )
    print("NAMD run completed")
    print(f"result_type={type(result).__name__}")
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=str,
        default=None,
        help="Path to a JSON configuration file. Defaults to the built-in example.",
    )
    parser.add_argument(
        "--mode",
        choices=("dry-run", "callback-demo", "namd"),
        default="dry-run",
        help="Execution mode. 'dry-run' does not require liblibra_core.",
    )
    return parser


def main(argv: list[str] | None = None) -> Any:
    args = build_parser().parse_args(argv)
    cfg = get_config(args.config)

    if args.mode == "dry-run":
        dry_run(cfg)
        return None
    if args.mode == "callback-demo":
        callback_demo(cfg)
        return None
    return run_namd(cfg)


if __name__ == "__main__":
    main()
