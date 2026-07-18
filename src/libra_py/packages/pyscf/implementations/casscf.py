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
import copy as pycopy
from dataclasses import dataclass
from typing import Any, List, Optional, Tuple, Union
import numpy as np
from pyscf import fci, gto, mcscf, scf

def _compute_ao_overlap(prev_mol: Any, curr_mol: Any) -> np.ndarray:
    return gto.intor_cross("int1e_ovlp", prev_mol, curr_mol)

from libra_py.packages.pyscf.interfaces import ES_Strategy, ES_Request, MolecularGeometry

BOHR_TO_ANG = 0.529177210903

@dataclass
class CASSCF_States:
    mol: Optional[Any] = None
    mf: Optional[Any] = None
    mc: Optional[Any] = None

class CASSCF(ES_Strategy):
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
        nelecas:int = 0, 
        nroots: int = 1,
        basis: str = "sto-3g",
        unit: str = "Bohr",
        charge: int = 0,
        cas_list: Optional[List[int]] = None,
    ) -> None:
        #setting up the initial state 
        self._mol: Optional[Any] = mol
        self._norbcas: int = norbcas
        self._nelecas: Union[int, Tuple[int, int]] = nelecas  
        self._nroots: int = nroots    #can only be >1 
        self._basis: str = basis
        self._unit: str = unit  #default to Bohr, but can be set to Angstrom
        self._charge: int = int(charge)
        self._cas_list: Optional[List[int]] = cas_list
        self._mf: Optional[Any] = None
        self._mc: Optional[Any] = None
        self._geom: Optional[MolecularGeometry] = None
        self._request: Optional[ES_Request] = None
        self._ao_overlap: Optional[np.ndarray] = None  # AO overlap between consecutive geoms for time-overlap computation
        self._state: CASSCF_States | None = None
        self._previous_state: CASSCF_States | None = None

    @staticmethod
    def _as_ci_vector_seq(ci_data: Any) -> Optional[tuple[np.ndarray, ...]]:
        if ci_data is None:
            return None
        if isinstance(ci_data, (list, tuple)):
            return tuple(np.asarray(vec) for vec in ci_data)
        return (np.asarray(ci_data),)

    @staticmethod
    def _normalize_mc(mc: Any) -> None:
        """Promote single-root mc result to list layout so all callers see a uniform interface."""
        if not hasattr(mc, "e_states"):
            mc.e_states = [float(mc.e_tot)]
        if not isinstance(mc.ci, list):
            mc.ci = [mc.ci]

    @staticmethod
    def _coerce_ci_roots(ci_data: Any) -> list[np.ndarray]:
        if ci_data is None:
            return []
        if isinstance(ci_data, (list, tuple)):
            return [np.asarray(vec) for vec in ci_data]
        return [np.asarray(ci_data)]

    @staticmethod
    def _phase_align_ci_roots(
        prev_roots: list[np.ndarray], curr_roots: list[np.ndarray]
    ) -> tuple[list[np.ndarray], list[np.ndarray]]:
        prev = [np.asarray(vec) for vec in prev_roots]
        curr = [np.asarray(vec) for vec in curr_roots]
        for idx in range(min(len(prev), len(curr))):
            if prev[idx].size == 0 or curr[idx].size == 0:
                continue
            overlap0 = np.vdot(prev[idx], curr[idx])
            if abs(overlap0) > 1e-14:
                phase = overlap0 / abs(overlap0)
                curr[idx] = curr[idx] * np.conjugate(phase)
        return prev, curr

    def _active_mo_coeff(self, mc_obj: Any) -> np.ndarray:
        mo_coeff = np.asarray(mc_obj.mo_coeff)
        ncore = int(getattr(mc_obj, "ncore", 0))
        ncas = int(getattr(mc_obj, "ncas", self._norbcas))
        return mo_coeff[:, ncore:ncore + ncas]

    def set_geom(self, geom: MolecularGeometry) -> None:
        self._geom = geom
        self._previous_state = self.get_state()

        self._mc = None
        self._ao_overlap = None

        charge: int = self._charge
        coords = getattr(geom, "coords_bohr", None)
        if coords is None:
            coords = np.asarray(getattr(geom, "coords_angstrom"), dtype=float) / BOHR_TO_ANG
        else:
            coords = np.asarray(coords, dtype=float)

        self._mol = gto.M(
            atom=";".join(
                f"{label} {coord[0]} {coord[1]} {coord[2]}"
                for label, coord in zip(geom.atom_labels, coords)
            ),
            basis=self._basis,
            unit=self._unit,
            charge=charge,
            spin=0,
        )

        prev_state = self.get_previous_state()
        if prev_state is not None and prev_state.mol is not None:
            self._ao_overlap = _compute_ao_overlap(prev_state.mol, self._mol)

        if prev_state is not None and prev_state.mf is not None:
            self._mf = scf.RHF(self._mol)
            try:
                self._mf.init_guess_by_mo(prev_state.mf.mo_coeff)
            except Exception:
                pass
            self._mf = self._mf.run(verbose=0)
        else:
            self._mf = scf.RHF(self._mol).run(verbose=0)

        self._state = CASSCF_States(mol=self._mol, mf=self._mf, mc=None)

    def get_geom(self) -> MolecularGeometry:
        if self._geom is None:
            raise ValueError("Geometry has not been set.")
        return self._geom

    def get_state(self) -> Optional[CASSCF_States]:
        return self._state

    def get_previous_state(self) -> Optional[CASSCF_States]:
        return self._previous_state

    def copy(self) -> "CASSCF":
        """Return an independent snapshot of the full strategy state."""
        return pycopy.deepcopy(self)

    def _ensure_mc(self, previous: Optional[ES_Strategy] = None) -> None:
        
        
        if self._mf is None:
            raise ValueError("HF must be run before computing CASSCF energies.")

        if self._mc is None:
            prev_state = self.get_previous_state()
            if prev_state is None and previous is not None and isinstance(previous, CASSCF):
                prev_state = previous.get_state()
            prev_mc = prev_state.mc if prev_state is not None else None
            self._mc = mcscf.CASSCF(self._mf, self._norbcas, self._nelecas)
            self._mc.fcisolver = fci.direct_spin0.FCI(self._mol)
            self._mc.fcisolver.nroots = self._nroots
            if self._nroots > 1:
                self._mc = self._mc.state_average_([1.0 / self._nroots] * self._nroots)  # equal weights as default; required for gradients in pyscf

            if prev_mc is not None:
                mo_coeff = getattr(prev_mc, "mo_coeff", None)
            else:
                mo_coeff = self._mf.mo_coeff
                if self._cas_list is not None:
                    mo_coeff = mcscf.sort_mo(self._mc, mo_coeff, self._cas_list)
            try:
                ci0 = prev_mc.ci if prev_mc is not None else None
            except AttributeError:
                ci0 = None

            self._mc.kernel(mo_coeff=mo_coeff, ci0=ci0)
            self._normalize_mc(self._mc)
            state = self.get_state()
            if state is not None:
                state.mc = self._mc

    def compute_H_el(self, previous: Optional[ES_Strategy] = None) -> np.ndarray:
        self._ensure_mc(previous)
        n_total = self._request.n_singlets if self._request is not None else self._nroots
        e_states: List[float] = self._mc.e_states
        if len(e_states) < n_total:
            raise IndexError(f"Requested {n_total} states, but only {len(e_states)} are available.")
        return np.asarray(e_states[:n_total], dtype=np.float64)

    def compute_energy(self, root: int) -> float:
        energies = self.compute_H_el()
        return float(energies[root])

    def compute_gradient(self, root: int = 0) -> np.ndarray:
        self._ensure_mc()
        e_states: List[float] = self._mc.e_states
        if root < 0 or root >= len(e_states):
            raise IndexError(f"Requested root {root}, but only {len(e_states)} roots are available.")
        return np.asarray(self._mc.nuc_grad_method(state=root).kernel())

    # time overlap
    def _require_time_overlap_state(self) -> None:
        if self._mol is None or self._mf is None:
            raise ValueError("Molecule and HF state are required for time-overlap computation.")
        if self._mc is None:
            raise ValueError("CASSCF must be run before computing time-overlap matrix.")
        if getattr(self._mc, "ci", None) is None:
            raise ValueError("CI vectors are required for time-overlap computation.")

    def _time_overlap_nroots(self, right: "CASSCF") -> int:
        if right._request is not None:
            return int(right._request.n_singlets)
        if self._request is not None:
            return int(self._request.n_singlets)
        return min(int(self._nroots), int(right._nroots))

    def compute_ao_overlap(self, right: CASSCF) -> np.ndarray:
        if not isinstance(right, CASSCF):
            raise TypeError(f"right must be CASSCF, got {type(right).__name__}")
        if self._mol is None or right._mol is None:
            raise ValueError("Both CASSCF objects need molecule states before AO overlap.")
        return _compute_ao_overlap(self._mol, right._mol)

    def _compute_active_mo_overlap(self, right: "CASSCF", ao_overlap: np.ndarray) -> np.ndarray:
        mo_left_act = self._active_mo_coeff(self._mc)
        mo_right_act = right._active_mo_coeff(right._mc)
        return mo_left_act.T @ ao_overlap @ mo_right_act

    def _ci_roots(self) -> list[np.ndarray]:
        roots = self._coerce_ci_roots(self._mc.ci)
        if not roots:
            raise ValueError("CI vectors are required for time-overlap computation.")
        return roots

    def _compute_ci_overlap_matrix(
        self,
        right: "CASSCF",
        nroots: int,
        active_mo_overlap: np.ndarray,
    ) -> np.ndarray:
        left_roots = self._ci_roots()
        right_roots = right._ci_roots()
        if len(left_roots) < nroots or len(right_roots) < nroots:
            raise ValueError(
                f"Requested {nroots} roots, but only {len(left_roots)} left "
                f"and {len(right_roots)} right roots are available."
            )

        overlap = np.zeros((nroots, nroots), dtype=float)
        for i in range(nroots):
            for j in range(nroots):
                overlap[i, j] = fci.addons.overlap(
                    left_roots[i],
                    right_roots[j],
                    self._mc.ncas,
                    self._mc.nelecas,
                    s=active_mo_overlap,
                )
        return np.asarray(np.real_if_close(overlap))

    @staticmethod
    def _align_time_overlap_phases(overlap: np.ndarray) -> np.ndarray:
        overlap = np.array(overlap, copy=True)
        for j in range(overlap.shape[1]):
            if overlap[j, j] < 0:
                overlap[:, j] = -overlap[:, j]
        return overlap

    def compute_time_overlap(self, right: CASSCF) -> np.ndarray:
        if not isinstance(right, CASSCF):
            raise TypeError(f"right must be CASSCF, got {type(right).__name__}")
        if self._mf is None or right._mf is None:
            raise ValueError("HF states are required for time-overlap computation.")
        if self._mol is None or right._mol is None:
            raise ValueError("Molecule states are required for time-overlap computation.")

        nroots = self._time_overlap_nroots(right)
        if nroots <= 0:
            raise ValueError(f"nroots must be positive, got {nroots}")

        ao_overlap = self.compute_ao_overlap(right)
        prev_casci = mcscf.CASCI(self._mf, self._norbcas, self._nelecas)
        prev_casci.fcisolver = fci.direct_spin0.FCISolver(self._mol)
        prev_casci.fcisolver.nroots = nroots
        h1prev, _ = prev_casci.get_h1eff(prev_casci.mo_coeff)
        h2prev = prev_casci.get_h2cas(prev_casci.mo_coeff)
        _, prev_roots = prev_casci.fcisolver.kernel(
            h1prev,
            h2prev,
            prev_casci.ncas,
            prev_casci.nelecas,
            nroots=nroots,
        )

        curr_casci = mcscf.CASCI(right._mf, right._norbcas, right._nelecas)
        curr_casci.fcisolver = fci.direct_spin0.FCISolver(right._mol)
        curr_casci.fcisolver.nroots = nroots
        h1curr, _ = curr_casci.get_h1eff(curr_casci.mo_coeff)
        h2curr = curr_casci.get_h2cas(curr_casci.mo_coeff)
        _, curr_roots = curr_casci.fcisolver.kernel(
            h1curr,
            h2curr,
            curr_casci.ncas,
            curr_casci.nelecas,
            nroots=nroots,
        )

        prev_roots = [np.asarray(vec) for vec in prev_roots[:nroots]]
        curr_roots = [np.asarray(vec) for vec in curr_roots[:nroots]]

        prev_act = np.asarray(prev_casci.mo_coeff)[:, : prev_casci.ncas]
        curr_act = np.asarray(curr_casci.mo_coeff)[:, : curr_casci.ncas]
        s12_mo = prev_act.T @ ao_overlap @ curr_act

        overlap = np.zeros((nroots, nroots), dtype=float)
        for i in range(nroots):
            for j in range(nroots):
                overlap[i, j] = fci.addons.overlap(
                    prev_roots[i],
                    curr_roots[j],
                    prev_casci.ncas,
                    prev_casci.nelecas,
                    s=s12_mo,
                )

        overlap = np.asarray(np.real_if_close(overlap))
        for j in range(nroots):
            if overlap[j, j] < 0:
                overlap[:, j] = -overlap[:, j]
        return overlap

    def time_overlap_matrix(self, nroots: int) -> np.ndarray:
        """Legacy same-object API: compare cached previous state to current state."""
        if nroots <= 0:
            raise ValueError(f"nroots must be positive, got {nroots}")
        prev_state = self.get_previous_state()
        if prev_state is None or prev_state.mc is None:
            raise ValueError("Previous CASSCF state must exist for time-overlap computation.")

        left = self.copy()
        left._mol = prev_state.mol
        left._mf = prev_state.mf
        left._mc = prev_state.mc
        left._state = CASSCF_States(mol=prev_state.mol, mf=prev_state.mf, mc=prev_state.mc)
        left._previous_state = None
        left._nroots = int(nroots)

        right = self.copy()
        right._nroots = int(nroots)
        return left.compute_time_overlap(right)

    def compute_nac_vectors(self) -> np.ndarray:
        use_etfs = True
        if self._mc is None:
            raise ValueError("CASSCF must be run before computing NAC vectors.")
        nstates = self._nroots
        natm = int(self._mol.natm)
        
        mc_nacs = self._mc.nac_method()
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

    
