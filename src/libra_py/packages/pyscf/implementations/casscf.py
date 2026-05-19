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
.. module:: pyscf.implementations.casscf
   :platform: Unix, Windows
   :synopsis: PySCF CASSCF implementation for Libra ES strategy.
.. moduleauthor::
       Jieyang Gu <jieyanggu792@gmail.com>

"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Any, List, Optional, Tuple, Union, Sequence
import numpy as np
from pyscf import fci, gto, mcscf, scf
from libra_py.packages.pyscf.interfaces import ElectronicStructureStrategy, MolecularGeometry


@dataclass
class CASSCFTrajectoryState:
    mol: Optional[Any] = None
    mf: Optional[Any] = None
    mc: Optional[Any] = None
    ao_overlap: Optional[np.ndarray] = None
    prev_mol: Optional[Any] = None
    prev_mf: Optional[Any] = None
    prev_mc: Optional[Any] = None

class CASSCF(ElectronicStructureStrategy):
    """PySCF-based CASSCF backend for the universal ES interface."""

    # 1) when a geom is set the HF is run; the MF is written to the member attribute
    #    `_mf`
    # 2) when the CASSCF energy is requested for a certain root, the CASSCF is run
    #    (number of roots requested = self._nroots). The resulting MC object is
    #    stored in the member attribute `mc`, which also contains energies of other
    #    states.

    def __init__(
        self,
        mol: Optional[Any] = None,
        norbcas: int = 0,
        nelecas: int = 0,
        nroots: int = 1,
        basis: str = "sto-3g",
        unit: str = "Angstrom",
        charge: int = 0,
        cas_list: Optional[List[int]] = None,
        ntraj: int = 1,
        use_prev_ci: bool = False,
    ) -> None:
        super().__init__(
            nroots=nroots,
            basis=basis,
            unit=unit,
            charge=int(charge),
            ntraj=ntraj,
        )

        self._norbcas: int = norbcas
        self._nelecas: Union[int, Tuple[int, int]] = nelecas
        self._cas_list: Optional[List[int]] = cas_list
        self._use_prev_ci: bool = use_prev_ci

    def _get_traj_state(self, traj_id: int) -> CASSCFTrajectoryState:
        state = self._traj_states[traj_id]
        if state is None:
            state = CASSCFTrajectoryState()
            self._traj_states[traj_id] = state
        return state

    def set_geom(self, geom: MolecularGeometry, traj_id: int = 0) -> None:
        state = self._get_traj_state(traj_id)

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

        state.mc = None
        state.ao_overlap = None

        if state.prev_mol is not None:
            state.ao_overlap = gto.intor_cross("int1e_ovlp", state.prev_mol, state.mol)

    def save_cache(self, traj_id: int = 0) -> None:
        # Cache the current state directly (None is fine for first geometry).
        state = self._get_traj_state(traj_id)
        state.prev_mol = state.mol
        state.prev_mc = state.mc
        state.prev_mf = state.mf

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
            raise ValueError("HF must be run before computing CASSCF energies.")
        if root < 0 or root >= self._nroots:
            raise IndexError(f"Requested root {root}, but only {self._nroots} roots are available.")

        if state.mc is None:
            mc = mcscf.CASSCF(state.mf, self._norbcas, self._nelecas)
            mc.fcisolver = fci.direct_spin0.FCI(state.mol)
            mc.fcisolver.nroots = self._nroots

            if self._nroots > 1:
                mc = mc.state_average_([1.0 / self._nroots] * self._nroots)

            mo_coeff = state.mf.mo_coeff
            if self._cas_list is not None:
                mo_coeff = mcscf.sort_mo(mc, mo_coeff, self._cas_list)

            ci0 = None
            if self._use_prev_ci and state.prev_mc is not None:
                ci0 = getattr(state.prev_mc, "ci", None)

            if ci0 is not None:
                mc.kernel(mo_coeff, ci0=ci0)
            else:
                mc.kernel(mo_coeff)

            state.mc = mc

        e_states: Optional[Sequence[float]] = getattr(state.mc, "e_states", None)
        if e_states is not None:
            if root < 0 or root >= len(e_states):
                raise IndexError(f"Requested root {root}, but only {len(e_states)} roots are available.")
            return float(np.asarray(e_states)[root])

        if root != 0:
            raise IndexError("Only root 0 is available for a single-state CASSCF calculation.")
        return float(state.mc.e_tot)

    def compute_gradient(self, root: int, traj_id: int = 0) -> np.ndarray:
        state = self._get_traj_state(traj_id)
        if state.mc is None:
            raise ValueError("CASSCF must be run before computing gradients.")
        e_states: Optional[Sequence[float]] = getattr(state.mc, "e_states", None)
        if e_states is not None:
            if root < 0 or root >= len(e_states):
                raise IndexError(f"Requested root {root}, but only {len(e_states)} roots are available.")
            return np.asarray(state.mc.nuc_grad_method(state=root).kernel())
        if root != 0:
            raise IndexError("Only root 0 is available for a single-state CASSCF calculation.")
        return np.asarray(state.mc.nuc_grad_method().kernel())

    def time_overlap_matrix(self, nroots: int, traj_id: int = 0) -> np.ndarray:
        state = self._get_traj_state(traj_id)
        if state.mc is None:
            raise ValueError("CASSCF must be run before computing time-overlap matrix.")
        if state.prev_mf is None or state.prev_mol is None:
            raise ValueError("Previous and current HF/molecule states are required for time-overlap computation.")
        if state.ao_overlap is None:
            raise ValueError("Cached AO overlap between consecutive geometries is not available.")
        if nroots <= 0:
            raise ValueError(f"nroots must be positive, got {nroots}")

        nroots = int(nroots)
        if state.mf is None or state.mol is None:
            raise ValueError("Current HF/molecule state is required for time-overlap computation.")

        prev_casci = mcscf.CASCI(state.prev_mf, self._norbcas, self._nelecas)
        prev_casci.fcisolver = fci.direct_spin0.FCISolver(state.prev_mol)
        prev_casci.fcisolver.nroots = nroots
        h1prev, _ = prev_casci.get_h1eff(prev_casci.mo_coeff)
        h2prev = prev_casci.get_h2cas(prev_casci.mo_coeff)
        _, prev_roots = prev_casci.fcisolver.kernel(h1prev, h2prev, prev_casci.ncas, prev_casci.nelecas, nroots=nroots)

        curr_casci = mcscf.CASCI(state.mf, self._norbcas, self._nelecas)
        curr_casci.fcisolver = fci.direct_spin0.FCISolver(state.mol)
        curr_casci.fcisolver.nroots = nroots
        h1curr, _ = curr_casci.get_h1eff(curr_casci.mo_coeff)
        h2curr = curr_casci.get_h2cas(curr_casci.mo_coeff)
        _, curr_roots = curr_casci.fcisolver.kernel(h1curr, h2curr, curr_casci.ncas, curr_casci.nelecas, nroots=nroots)

        if not isinstance(prev_roots, (list, tuple)):
            prev_roots = [prev_roots]
        if not isinstance(curr_roots, (list, tuple)):
            curr_roots = [curr_roots]
        if len(prev_roots) < nroots or len(curr_roots) < nroots:
            raise ValueError(
                f"Requested {nroots} roots, but only {len(prev_roots)} previous and {len(curr_roots)} current roots are available."
            )

        prev_roots = [np.asarray(v) for v in prev_roots[:nroots]]
        curr_roots = [np.asarray(v) for v in curr_roots[:nroots]]

        mo_prev_act = np.asarray(prev_casci.mo_coeff)[:, prev_casci.ncore : prev_casci.ncore + prev_casci.ncas]
        mo_curr_act = np.asarray(curr_casci.mo_coeff)[:, curr_casci.ncore : curr_casci.ncore + curr_casci.ncas]
        s12_mo = mo_prev_act.T @ state.ao_overlap @ mo_curr_act

        overlap = np.zeros((nroots, nroots), dtype=float)
        nelecas = prev_casci.nelecas
        for i in range(nroots):
            for j in range(nroots):
                overlap[i, j] = fci.addons.overlap(
                    prev_roots[i],
                    curr_roots[j],
                    prev_casci.ncas,
                    nelecas,
                    s=s12_mo,
                )

        overlap = np.asarray(np.real_if_close(overlap))

        for i in range(nroots):
            if overlap[i, i] < 0:
                overlap[i, :] = -overlap[i, :]

        return overlap
    
    def compute_nac_vectors(self, use_etfs: bool = True, traj_id: int = 0) -> np.ndarray:
        state = self._get_traj_state(traj_id)
        if state.mc is None:
            raise ValueError("CASSCF must be run before computing NAC vectors.")
        if state.mol is None:
            raise ValueError("Current molecule is required before computing NAC vectors.")

        nstates = self._nroots
        natm = int(state.mol.natm)
        
        mc_nacs = state.mc.nac_method()
        nacv = np.zeros((nstates, nstates, natm, 3), dtype=np.float64)

        for ket in range(nstates):
            for bra in range(nstates):
                if bra == ket:
                    continue
                nacv[bra, ket] = np.asarray(
                    mc_nacs.kernel(state=(ket, bra), use_etfs=use_etfs, mult_ediff=False),
                    dtype=np.float64,
                )

        return nacv
