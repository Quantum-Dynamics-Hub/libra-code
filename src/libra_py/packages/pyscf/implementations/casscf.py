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
from typing import Any, Optional, Tuple, Union, Sequence
import numpy as np
from pyscf import fci, gto, mcscf, scf
from ..interfaces import ElectronicStructureStrategy, MolecularGeometry


@dataclass
class Cache:
    prev_mol: Optional[Any] = None
    prev_mc: Optional[Any] = None
    prev_mf: Optional[Any] = None
    prev_ci_vecs: Optional[Sequence[np.ndarray]] = None

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
        nelecas:int = 0, 
        nroots: int = 1,
        basis: str = "sto-3g",
        unit: str = "Angstrom",
        charge: int = 0,
    ) -> None:
        #setting up the initial state 
        self._mol: Optional[Any] = mol
        self._norbcas: int = norbcas
        self._nelecas: Union[int, Tuple[int, int]] = nelecas  
        self._nroots: int = nroots
        self._basis: str = basis
        self._unit: str = unit  # default to Angstrom 
        self._charge: int = int(charge)
        self._mf: Optional[Any] = None
        self._mc: Optional[Any] = None
        self._ci_vecs: Optional[Sequence[np.ndarray]] = None  # list of CI vectors as numpy arrays, one per root
        self._ao_overlap: Optional[np.ndarray] = None  # AO overlap between consecutive geoms for time-overlap computation
        self._cache: Cache = Cache()

    @staticmethod
    def _as_ci_vector_list(ci_data: Any) -> Optional[list[np.ndarray]]:
        if ci_data is None:
            return None
        if isinstance(ci_data, (list, tuple)):
            return [np.asarray(vec).copy() for vec in ci_data]
        return [np.asarray(ci_data).copy()]

    def set_geom_and_run_hf(self, geom: MolecularGeometry) -> None:
        # Cache the current state directly (None is fine for first geometry).
        self._cache.prev_mol = self._mol
        self._cache.prev_mc = self._mc
        self._cache.prev_mf = self._mf
        self._cache.prev_ci_vecs = self._ci_vecs

        # Clear the current state.

        self._mc = None
        self._ci_vecs = None
        self._ao_overlap = None

        # Set up the new molecule and run HF.

        charge: int = self._charge
        self._mol = gto.M(
            atom=";".join(
                f"{label} {coord[0]} {coord[1]} {coord[2]}"
                for label, coord in zip(geom.atom_labels, geom.coords_angstrom)
            ),
            basis=self._basis,
            unit=self._unit,
            charge=charge,
            spin=0,
        )
        #compute the AO overlap between new geom and the previous one 
        if self._cache.prev_mol is not None:
            self._ao_overlap = gto.intor_cross("int1e_ovlp", self._cache.prev_mol, self._mol)

        self._mf = scf.RHF(self._mol).run(verbose=0)

    def compute_energy(self, root: int) -> float:
        if self._mf is None:
            raise ValueError("HF must be run before computing CASSCF energies.")

        if self._mc is None:
            self._mc = mcscf.CASSCF(self._mf, self._norbcas, self._nelecas)
            self._mc.fcisolver = fci.direct_spin0.FCI(self._mol)
            self._mc.fcisolver.nroots = self._nroots
            if self._nroots > 1:
                self._mc = self._mc.state_average_([1.0 / self._nroots] * self._nroots)  # equal weights as default; required for gradients in pyscf

            if self._cache.prev_mc is not None:
                mo_coeff = getattr(self._cache.prev_mc, "mo_coeff", None)
            else:
                mo_coeff = self._mf.mo_coeff

            if mo_coeff is not None:
                self._mc.kernel(mo_coeff=mo_coeff)
            else:
                self._mc.kernel()

        if self._ci_vecs is None:
            self._ci_vecs = self._as_ci_vector_list(getattr(self._mc, "ci", None))

        e_states: Optional[Sequence[float]] = getattr(self._mc, "e_states", None)
        if e_states is not None:
            if root < 0 or root >= len(e_states):
                raise IndexError(f"Requested root {root}, but only {len(e_states)} roots are available.")
            return float(np.asarray(e_states)[root])
        if root != 0:
            raise IndexError("Only root 0 is available for a single-state CASSCF calculation.")
        return float(self._mc.e_tot)

    def compute_gradient(self, root: int) -> np.ndarray:
        if self._mc is None:
            raise ValueError("CASSCF must be run before computing gradients.")
        e_states: Optional[Sequence[float]] = getattr(self._mc, "e_states", None)
        if e_states is not None:
            if root < 0 or root >= len(e_states):
                raise IndexError(f"Requested root {root}, but only {len(e_states)} roots are available.")
            return np.asarray(self._mc.nuc_grad_method(state=root).kernel())
        if root != 0:
            raise IndexError("Only root 0 is available for a single-state CASSCF calculation.")
        return np.asarray(self._mc.nuc_grad_method().kernel())

    def time_overlap_matrix(self, nroots: int) -> np.ndarray:
        if self._mc is None:
            raise ValueError("CASSCF must be run before computing time-overlap matrix.")
        if self._cache is None or self._cache.prev_mf is None or self._cache.prev_mol is None:
            raise ValueError("Previous and current HF/molecule states are required for time-overlap computation.")
        if self._ao_overlap is None:
            raise ValueError("Cached AO overlap between consecutive geometries is not available.")
        if nroots <= 0:
            raise ValueError(f"nroots must be positive, got {nroots}")

        nroots = int(nroots)
        if self._mf is None or self._mol is None:
            raise ValueError("Current HF/molecule state is required for time-overlap computation.")

        prev_casci = mcscf.CASCI(self._cache.prev_mf, self._norbcas, self._nelecas)
        prev_casci.fcisolver = fci.direct_spin0.FCISolver(self._cache.prev_mol)
        prev_casci.fcisolver.nroots = nroots
        h1prev, _ = prev_casci.get_h1eff(prev_casci.mo_coeff)
        h2prev = prev_casci.get_h2cas(prev_casci.mo_coeff)
        _, prev_roots = prev_casci.fcisolver.kernel(h1prev, h2prev, prev_casci.ncas, prev_casci.nelecas, nroots=nroots)

        curr_casci = mcscf.CASCI(self._mf, self._norbcas, self._nelecas)
        curr_casci.fcisolver = fci.direct_spin0.FCISolver(self._mol)
        curr_casci.fcisolver.nroots = nroots
        h1curr, _ = curr_casci.get_h1eff(curr_casci.mo_coeff)
        h2curr = curr_casci.get_h2cas(curr_casci.mo_coeff)
        _, curr_roots = curr_casci.fcisolver.kernel(h1curr, h2curr, curr_casci.ncas, curr_casci.nelecas, nroots=nroots)

        if len(prev_roots) < nroots or len(curr_roots) < nroots:
            raise ValueError(
                f"Requested {nroots} roots, but only {len(prev_roots)} previous and {len(curr_roots)} current roots are available."
            )

        prev_roots = [np.asarray(v) for v in prev_roots[:nroots]]
        curr_roots = [np.asarray(v) for v in curr_roots[:nroots]]

        mo_prev_act = np.asarray(prev_casci.mo_coeff)[:, : prev_casci.ncas]
        mo_curr_act = np.asarray(curr_casci.mo_coeff)[:, : curr_casci.ncas]
        s12_mo = mo_prev_act.T @ self._ao_overlap @ mo_curr_act

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
       
