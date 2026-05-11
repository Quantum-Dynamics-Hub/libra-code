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
from pyscf import ci

from libra_py.packages.pyscf.interfaces import ElectronicStructureStrategy, MolecularGeometry
#from libra_py.packages.pyscf.utils import run_rhf_for_geometry

@dataclass
class CISDCache:
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
    ) -> None:
        super().__init__(mol=mol, nroots=nroots, basis=basis, unit=unit, charge=int(charge))
        self._ci: Optional[Any] = None
        self._cache: CISDCache = CISDCache()
        self._use_prev_ci: bool = use_prev_ci

    def save_cache(self) -> None:
        self._cache.prev_mol = self._mol
        self._cache.prev_mf = self._mf
        self._cache.prev_ci = self._ci

    def run_hf(self) -> None:
        self._ci = None
        self._ao_overlap = None

        dm0 = None
        if self._cache.prev_mf is not None:
            try:
                dm0 = self._cache.prev_mf.make_rdm1()
            except Exception:
                dm0 = None

        self._mol, self._mf, self._ao_overlap = run_rhf_for_geometry(
            self._geom,
            basis=self._basis,
            unit=self._unit,
            charge=self._charge,
            prev_mol=self._cache.prev_mol,
            overlap_integral="int1e_ovlp",
            spin=0,
            verbose=0,
            dm0=dm0,
        )


    def compute_energy(self, root: int) -> float:
        if self._mf is None:
           raise ValueError("HF must be run before computing CISD energies.")
        if root < 0 or root >= self._nroots:
           raise IndexError(f"Requested root {root}, but only {self._nroots} roots are available.")

        if self._ci is None:
           self._ci = ci.cisd.CISD(self._mf)
           self._ci.nroots = self._nroots
           self._ci.verbose = 0
           # Use CI coefficients from previous step if requested
           ci0 = None
           if self._use_prev_ci and self._cache.prev_ci is not None:
               ci0 = getattr(self._cache.prev_ci, "ci", None)
           if ci0 is not None:
               self._ci.kernel(ci0=ci0)
           else:
               self._ci.kernel()

        e_tot: Optional[Sequence[float]] = getattr(self._ci, "e_tot", None)
        if isinstance(e_tot, (list, tuple, np.ndarray)):
           return float(np.asarray(e_tot)[root])

        if root != 0:
           raise IndexError("Only root 0 is available for this CISD calculation.")
        return float(e_tot)

    def compute_gradient(self, root: int) -> np.ndarray:
        if self._ci is None:
           self.compute_energy(root=0)
        # normalize CI roots into a list
        ci_vectors = self._ci.ci
        if isinstance(ci_vectors, (list, tuple)):
           ci_roots = [np.asarray(vec) for vec in ci_vectors]
        else:
           ci_roots = [np.asarray(ci_vectors)]
        if root < 0 or root >= len(ci_roots):
           raise IndexError(f"Requested root {root}, but only {len(ci_roots)} roots are available.")

        if len(ci_roots) == 1:
           if root != 0:
               raise IndexError("Only root 0 is available for a single-state CISD calculation.")
           return np.asarray(self._ci.nuc_grad_method().kernel())
        return np.asarray(self._ci.nuc_grad_method().kernel(state=root))

    def time_overlap_matrix(self, nroots: int) -> np.ndarray:
        if nroots <= 0:
           raise ValueError(f"nroots must be positive, got {nroots}")
        if nroots > self._nroots:
           raise ValueError(f"Requested {nroots} roots, but CISD is configured for {self._nroots} roots")

        if self._ci is None:
           raise ValueError("CISD must be run before computing time-overlap matrix.")

        nroots = int(nroots)

        if self._cache.prev_ci is None or self._cache.prev_mf is None or self._cache.prev_mol is None:
           raise ValueError("Previous CISD state must exist for time-overlap computation.")
        if self._mf is None or self._mol is None:
           raise ValueError("Current CISD state must exist for time-overlap computation.")
        if self._ao_overlap is None:
           raise ValueError("Cached AO overlap between consecutive geometries is not available.")

        s12_mo = reduce(np.dot, (self._cache.prev_mf.mo_coeff.T, self._ao_overlap, self._mf.mo_coeff))

        prev_ci_roots = self._cache.prev_ci.ci
        curr_ci_roots = self._ci.ci
        if isinstance(prev_ci_roots, (list, tuple)):
           prev_ci_list = [np.asarray(v) for v in prev_ci_roots]
        else:
           prev_ci_list = [np.asarray(prev_ci_roots)]

        if isinstance(curr_ci_roots, (list, tuple)):
           curr_ci_list = [np.asarray(v) for v in curr_ci_roots]
        else:
           curr_ci_list = [np.asarray(curr_ci_roots)]

        nmo = self._ci.nmo
        nelec = self._mol.nelectron // 2

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
