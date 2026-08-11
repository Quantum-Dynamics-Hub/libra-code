# *********************************************************************************
# * Copyright (C) 2026 Jieyang Gu <jieyanggu792@gmail.com>
# * This file is distributed under the terms of the GNU General Public License
# * as published by the Free Software Foundation, either version 3 of
# * the License, or (at your option) any later version.
# * See the file LICENSE in the root directory of this distribution
# * or <http://www.gnu.org/licenses/>.
# *
# *********************************************************************************/
"""Two-phase LiF CASSCF workflow demo using the ES interface and adapter.

Phase 1:
    Build a simple precomputed trajectory by sweeping the Li-F bond length over a
    range. At each point, request adiabatic energies and the time-overlap matrix
    between adjacent geometries through the adapter.

Phase 2:
    Propagate a simple electronic amplitude over the precomputed trajectory using
    the stored energies and overlaps as a lightweight surrogate for an NBRA
    dynamics step.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import numpy as np


def _prepend_repo_root() -> None:
    file_path = Path(__file__).resolve()
    for parent in file_path.parents:
        if (parent / "src" / "libra_py" / "__init__.py").is_file():
            repo_root = str(parent / "src")
            if repo_root not in sys.path:
                sys.path.insert(0, repo_root)
            return


if __name__ == "__main__" and __package__ is None:
    _prepend_repo_root()

from libra_py.packages.pyscf.adapter import LibraESAdapter
from libra_py.packages.pyscf.implementations.casscf import CASSCF
from libra_py.packages.pyscf.interfaces import ES_Request


class FakeNham:
    def __init__(self) -> None:
        self.values: dict[str, object] = {}

    def set_ham_adi_by_val(self, value):
        self.values["ham"] = value

    def set_ham_dia_by_val(self, value):
        self.values["soc"] = value

    def set_d1ham_adi_by_val(self, value):
        self.values["d1"] = value

    def set_dc1_adi_by_val(self, value):
        self.values["dc1"] = value

    def set_time_overlap_adi_by_val(self, value):
        self.values["time"] = value


def bond_length_vs_time(t: float) -> float:
    """Simple toy trajectory: the bond length increases linearly with time."""
    return 6.0 + 4.0 * t


def phase1_tabulate(adapter: LibraESAdapter, request: ES_Request, times: np.ndarray) -> list[dict[str, Any]]:
    table: list[dict[str, Any]] = []
    previous_strategy = None

    for _, t in enumerate(times):
        distance = bond_length_vs_time(float(t))
        q = np.array([0.0, 0.0, 0.0, 0.0, 0.0, distance], dtype=np.float64)
        nham = FakeNham()
        params = {
            "nham": nham,
            "request": request,
            "previous_strategy": previous_strategy,
        }

        result = adapter.strategy_compute_adi(q, params)
        assert result.H_el is not None
        assert result.time_overlap is not None

        entry = {
            "time": float(t),
            "distance": float(distance),
            "energies": np.asarray(result.H_el, dtype=np.float64),
            "overlap": np.asarray(result.time_overlap, dtype=np.float64),
        }
        table.append(entry)

        previous_strategy = adapter.backend.copy()

    return table


def phase2_propagate(table: list[dict[str, Any]], dt: float = 0.1) -> list[np.ndarray]:
    coeffs = np.zeros(len(table[0]["energies"]), dtype=np.complex128)
    coeffs[0] = 1.0 + 0.0j

    populations: list[np.ndarray] = []
    for idx, entry in enumerate(table):
        if idx == 0:
            populations.append(np.abs(coeffs) ** 2)
            continue

        overlap = np.asarray(table[idx]["overlap"], dtype=np.complex128)
        prev_energies = np.asarray(table[idx - 1]["energies"], dtype=np.float64)
        curr_energies = np.asarray(entry["energies"], dtype=np.float64)
        mean_energies = 0.5 * (prev_energies + curr_energies)
        phase = np.exp(-1j * mean_energies * dt)
        coeffs = phase * (overlap @ coeffs)
        populations.append(np.abs(coeffs) ** 2)

    return populations


def main() -> None:
    nstates = 2
    times = np.linspace(0.0, 1.0, 5)
    basis_dict = {"Li": "sto-3g", "F": "6-311+g*"}
    cas_list = [4, 7, 11, 14, 17]

    backend = CASSCF(
        norbcas=5,
        nelecas=2,
        nroots=nstates,
        basis=basis_dict,
        unit="Bohr",
        charge=0,
        cas_list=cas_list,
    )
    request = ES_Request(
        n_singlets=nstates,
        n_triplets=0,
        H_soc=False,
        gradient_state=None,
        hessian_state=None,
        nacv=True,
        time_overlap=True,
    )
    adapter = LibraESAdapter(
        backend=backend,
        atom_labels=["Li", "F"],
        request_template=request,
    )

    print("Phase 1: tabulating precomputed trajectory")
    table = phase1_tabulate(adapter, request, times)
    for entry in table:
        print(
            f"t={entry['time']:.2f}  bond={entry['distance']:.3f} bohr  "
            f"energies={entry['energies']}  overlap_diag={np.diag(entry['overlap'])}"
        )

    print("\nPhase 2: propagating amplitudes over the precomputed trajectory")
    populations = phase2_propagate(table, dt=0.1)
    for idx, pop in enumerate(populations):
        print(f"step {idx}: populations={pop}")


if __name__ == "__main__":
    main()
