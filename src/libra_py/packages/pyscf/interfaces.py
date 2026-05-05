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
.. module:: pyscf.interfaces
   :platform: Unix, Windows
   :synopsis: Core interface definitions for electronic structure strategies.

   The interface is backend-agnostic: implementations may wrap PySCF, DFTB+,
   CP2K, or any other quantum-chemistry code.  All Libra-specific types
   (CMATRIX, MATRIX, etc.) live exclusively in the adapter layer so that
   strategy implementations never depend on liblibra_core.

.. moduleauthor::
       Jieyang Gu <jieyanggu792@gmail.com>

"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Any, Optional

import numpy as np


@dataclass
class MolecularGeometry:
    """Nuclear geometry in Angstrom."""
    atom_labels: list[str]
    coords_angstrom: np.ndarray  # shape (natoms, 3)


class ElectronicStructureStrategy(ABC):
    """Base interface for electronic structure backends.

    Implementations are allowed to store any backend-specific state internally.

    The interface is intentionally minimal: it only describes the operations
    required by a consumer (e.g. NAMD adapter) and does not prescribe how a
    backend achieves those operations.

    **Required** abstract methods/properties must be implemented by every
    concrete strategy.  **Optional** methods have default implementations
    that signal "not available"; the adapter queries capability flags before
    calling them.
    """

    # ------------------------------------------------------------------
    #  Metadata (required)
    # ------------------------------------------------------------------

    def __init__(
        self,
        mol: Optional[Any] = None,
        nroots: int = 1,
        basis: str = "sto-3g",
        unit: str = "Angstrom",
        charge: int = 0,
    ) -> None:
        self._mol = mol
        self._nroots = nroots
        self._basis = basis
        self._unit = unit
        self._charge = charge
        self._mf = None
        self._ao_overlap = None
        self._geom = None

    @property
    def nstates(self) -> int:
        return self._nroots

    # ------------------------------------------------------------------
    #  Core computation (required)
    # ------------------------------------------------------------------

    def save_cache(self) -> None:
        pass

    def set_geom(self, geom: MolecularGeometry) -> None:
        self._geom = geom

    @abstractmethod
    def run_hf(self) -> None:
        pass

    def set_geom_and_run_hf(self, geom: MolecularGeometry) -> None:
        self.save_cache()
        self.set_geom(geom)
        self.run_hf()

    @abstractmethod
    def compute_energy(self, root: int) -> float:
        """Return the total energy (Hartree) for *root*."""

    @abstractmethod
    def compute_gradient(self, root: int) -> np.ndarray:
        """Return the nuclear gradient for *root*.

        Returns
        -------
        np.ndarray
            Shape ``(natoms, 3)`` in Hartree/Bohr.
        """


    def time_overlap_matrix(self, nroots: int) -> np.ndarray:
        """Return the time-overlap matrix ``<psi_i(t)|psi_j(t+dt)>``.

        Returns
        -------
        np.ndarray
            Shape ``(nroots, nroots)``.
        """
        raise NotImplementedError(
            "This backend does not provide explicit NAC vectors; "
            "use time-overlap-based NACs instead."
        )

    def compute_nac_vectors(self, **kwargs: Any) -> np.ndarray:
        """Return all NAC vectors ``d_{ij}`` between states.

        Returns
        -------
        np.ndarray
            Shape ``(nstates, nstates, natoms, 3)`` in 1/Bohr.
            Only off-diagonal elements are meaningful.

        Raises
        ------
        NotImplementedError
            If the backend does not support explicit NAC vectors.
        """
        raise NotImplementedError(
            "This backend does not provide explicit NAC vectors; "
            "use time-overlap-based NACs instead."
        )
