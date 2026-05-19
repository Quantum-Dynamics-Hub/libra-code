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
    ao_overlap: Optional[np.ndarray] = None
    prev_mol: Optional[Any] = None
    prev_mf: Optional[Any] = None
    prev_ci: Optional[Any] = None


class CISD(ElectronicStructureStrategy):
    """PySCF-based CISD backend for the universal ES interface."""

    def __init__(
        self,
        mol: Optional[Any] = None,
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
            ntraj=ntraj,
        )

        self._use_prev_ci: bool = use_prev_ci

    def _get_traj_state(self, traj_id: int) -> CISDTrajectoryState:
        state = self._traj_states[traj_id]
        if state is None:
            state = CISDTrajectoryState()
            self._traj_states[traj_id] = state
        return state

    def set_geom(self, geom: MolecularGeometry, traj_id: int = 0) -> None:
        state = self._get_traj_state(traj_id)
        state.ci = None
        state.ao_overlap = None

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

        if state.prev_mol is not None:
            state.ao_overlap = gto.intor_cross("int1e_ovlp", state.prev_mol, state.mol)

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

    def compute_energy(self, root: int, traj_id: int = 0) -> float:
        state = self._get_traj_state(traj_id)

        if state.mf is None:
            raise ValueError("HF must be run before computing CISD energies.")
        if root < 0 or root >= self._nroots:
            raise IndexError(f"Requested root {root}, but only {self._nroots} roots are available.")

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

        e_tot: Optional[Sequence[float]] = getattr(state.ci, "e_tot", None)
        if isinstance(e_tot, (list, tuple, np.ndarray)):
            return float(np.asarray(e_tot)[root])

        if root != 0:
            raise IndexError("Only root 0 is available for this CISD calculation.")
        return float(e_tot)

    def compute_gradient(self, root: int, traj_id: int = 0) -> np.ndarray:
        state = self._get_traj_state(traj_id)
        if state.ci is None:
            self.compute_energy(root=0, traj_id=traj_id)

        ci_vectors = state.ci.ci
        if isinstance(ci_vectors, (list, tuple)):
            ci_roots = [np.asarray(vec) for vec in ci_vectors]
        else:
            ci_roots = [np.asarray(ci_vectors)]
        if root < 0 or root >= len(ci_roots):
            raise IndexError(f"Requested root {root}, but only {len(ci_roots)} roots are available.")

        if len(ci_roots) == 1:
            if root != 0:
                raise IndexError("Only root 0 is available for a single-state CISD calculation.")
            return np.asarray(state.ci.nuc_grad_method().kernel())
        return np.asarray(state.ci.nuc_grad_method().kernel(state=root))

    def time_overlap_matrix(self, nroots: int, traj_id: int = 0) -> np.ndarray:
        state = self._get_traj_state(traj_id)

        if nroots <= 0:
            raise ValueError(f"nroots must be positive, got {nroots}")
        if nroots > self._nroots:
            raise ValueError(f"Requested {nroots} roots, but CISD is configured for {self._nroots} roots")
        if state.ci is None:
            raise ValueError("CISD must be run before computing time-overlap matrix.")

        nroots = int(nroots)

        if state.prev_ci is None or state.prev_mf is None or state.prev_mol is None:
            raise ValueError("Previous CISD state must exist for time-overlap computation.")
        if state.mf is None or state.mol is None:
            raise ValueError("Current CISD state must exist for time-overlap computation.")
        if state.ao_overlap is None:
            raise ValueError("Cached AO overlap between consecutive geometries is not available.")

        s12_mo = reduce(np.dot, (state.prev_mf.mo_coeff.T, state.ao_overlap, state.mf.mo_coeff))

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

        if len(prev_ci_list) < nroots or len(curr_ci_list) < nroots:
            raise ValueError(
                f"Requested {nroots} roots, but only {len(prev_ci_list)} previous and {len(curr_ci_list)} current roots are available."
            )

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

        return np.asarray(np.real_if_close(overlap))
