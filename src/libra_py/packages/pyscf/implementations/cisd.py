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
.. module:: pyscf.implementations.cisd
   :platform: Unix, Windows
   :synopsis: PySCF CISD implementation for Libra ES strategy.
.. moduleauthor::
       Jieyang Gu <jieyanggu792@gmail.com>

"""
from __future__ import annotations

from dataclasses import dataclass
from functools import reduce
from typing import Any, Optional, Sequence

import numpy as np
from pyscf import ci, gto, scf

from libra_py.packages.pyscf.interfaces import ElectronicStructureStrategy, MolecularGeometry


@dataclass
class CISDTrajectoryState:
    mol: Optional[Any] = None
    mf: Optional[Any] = None
    ci: Optional[Any] = None
    prev_mol: Optional[Any] = None
    prev_mf: Optional[Any] = None
    prev_ci: Optional[Any] = None


class CISD(ElectronicStructureStrategy):
    """PySCF-based CISD backend for the universal ES interface."""

    def __init__(
        self,
        nroots: int = 1,
        basis: str = "sto-3g",
        unit: str = "Angstrom",
        charge: int = 0,
        use_prev_ci: bool = False,
        ntraj: int = 1,
    ) -> None:
        super().__init__(
            nroots=nroots,
            basis=basis,
            unit=unit,
            charge=int(charge),
        )

        if ntraj <= 0:
            raise ValueError(f"ntraj must be positive, got {ntraj}")
        self._nroots: int = nroots
        self._use_prev_ci: bool = use_prev_ci
        self._traj_states: list[Optional[CISDTrajectoryState]] = [None] * ntraj

    def _get_traj_state(self, traj_id: int) -> CISDTrajectoryState:
        state = self._traj_states[traj_id]
        if state is None:
            state = CISDTrajectoryState()
            self._traj_states[traj_id] = state
        return state

    def set_geom(self, geom: MolecularGeometry, traj_id: int = 0) -> None:
        if traj_id == 0:
            super().set_geom(geom)
            nroots = self.get_nroots()
            self._energies = [None] * nroots
            self._gradients = [None] * nroots

        state = self._get_traj_state(traj_id)
        state.ci = None

        state.mol = gto.M(
            atom=";".join(
                f"{label} {coord[0]} {coord[1]} {coord[2]}"
                for label, coord in zip(geom.atom_labels, geom.coords_angstrom)
            ),
            basis=self._basis,
            unit=self._unit,
            charge=self._charge,
            spin=0,
        )

    def save_cache(self, traj_id: int = 0) -> None:
        state = self._get_traj_state(traj_id)
        state.prev_mol = state.mol
        state.prev_mf = state.mf
        state.prev_ci = state.ci

    def run_hf(self, traj_id: int = 0) -> None:
        state = self._get_traj_state(traj_id)

        if state.mol is None:
            raise ValueError(f"Geometry for trajectory {traj_id} has not been set.")

        dm0 = None
        if state.prev_mf is not None:
            try:
                dm0 = state.prev_mf.make_rdm1()
            except Exception:
                dm0 = None

        mf = scf.RHF(state.mol)
        mf.verbose = 0
        if dm0 is not None:
            mf.kernel(dm0=dm0)
        else:
            mf.kernel()

        state.mf = mf

    def _run_cisd(self, traj_id: int = 0) -> Any:
        state = self._get_traj_state(traj_id)

        if state.ci is None:
            cisd = ci.cisd.CISD(state.mf)
            cisd.nroots = self._nroots
            cisd.verbose = 0

            ci0 = None
            if self._use_prev_ci and state.prev_ci is not None:
                ci0 = getattr(state.prev_ci, "ci", None)

            if ci0 is not None:
                cisd.kernel(ci0=ci0)
            else:
                cisd.kernel()

            state.ci = cisd

        return state.ci

    def _compute_energy(self, root: int, traj_id: int = 0) -> float:
        cisd = self._run_cisd(traj_id=traj_id)

        e_tot: Optional[Sequence[float]] = getattr(cisd, "e_tot", None)
        if isinstance(e_tot, (list, tuple, np.ndarray)):
            return float(np.asarray(e_tot)[root])
        return float(e_tot)

    def compute_energies(self) -> None:
        cisd = self._run_cisd(traj_id=0)
        e_tot: Optional[Sequence[float]] = getattr(cisd, "e_tot", None)
        if isinstance(e_tot, (list, tuple, np.ndarray)):
            self._energies = [float(e) for e in np.asarray(e_tot)[:self._nroots]]
            return

        self._energies[0] = float(cisd.e_tot)

    def compute_energy(self, root: int, traj_id: int = 0) -> float:
        if traj_id != 0:
            return self._compute_energy(root, traj_id=traj_id)
        return self.get_energy(root)

    def _compute_gradient(self, root: int, traj_id: int = 0) -> np.ndarray:
        cisd = self._run_cisd(traj_id=traj_id)

        ci_vectors = cisd.ci
        if isinstance(ci_vectors, (list, tuple)):
            return np.asarray(cisd.nuc_grad_method().kernel(state=root))
        return np.asarray(cisd.nuc_grad_method().kernel())

    def compute_gradient(self, root: int, traj_id: int = 0) -> None:
        if traj_id == 0 and self._energies[root] is None:
            self.compute_energies()
        gradient = self._compute_gradient(root, traj_id=traj_id)
        if traj_id == 0:
            self._gradients[root] = gradient

    def time_overlap_matrix(self, nroots: int, traj_id: int = 0) -> np.ndarray:
        state = self._get_traj_state(traj_id)
        nroots = int(nroots)

        if state.prev_ci is None or state.prev_mf is None or state.prev_mol is None:
            raise ValueError("Previous CISD state must exist for time-overlap computation.")

        ao_overlap = gto.intor_cross("int1e_ovlp", state.prev_mol, state.mol)
        s12_mo = reduce(np.dot, (state.prev_mf.mo_coeff.T, ao_overlap, state.mf.mo_coeff))

        prev_ci_roots = state.prev_ci.ci
        curr_ci_roots = state.ci.ci
        if isinstance(prev_ci_roots, (list, tuple)):
            prev_ci_list = [np.asarray(v) for v in prev_ci_roots]
        else:
            prev_ci_list = [np.asarray(prev_ci_roots)]

        if isinstance(curr_ci_roots, (list, tuple)):
            curr_ci_list = [np.asarray(v) for v in curr_ci_roots]
        else:
            curr_ci_list = [np.asarray(curr_ci_roots)]

        nmo = state.ci.nmo
        nelec = state.mol.nelectron // 2

        overlap = np.zeros((nroots, nroots), dtype=float)
        for i in range(nroots):
            for j in range(nroots):
                overlap[i, j] = ci.cisd.overlap(
                    prev_ci_list[i],
                    curr_ci_list[j],
                    nmo,
                    nelec,
                    s12_mo,
                )

        overlap = np.asarray(np.real_if_close(overlap))

        for i in range(nroots):
            if overlap[i, i] < 0:
                overlap[i, :] = -overlap[i, :]

        return overlap
