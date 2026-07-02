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
from turtle import left, right

import numpy as np


@dataclass
class MolecularGeometry:
    """Nuclear geometry in Angstrom."""
    atom_labels: list[str]
    coords_angstrom: np.ndarray  # shape (natoms, 3)


class ElectronicStructureStrategy(ABC):
    """Base interface for electronic structure backends.

    Implementations are allowed to store backend-specific state for one
    computed electronic-structure snapshot internally.

    **Required** abstract methods must be implemented by every
    concrete strategy.  
    """

    # ------------------------------------------------------------------
    #  Metadata (required)
    # ------------------------------------------------------------------
    def __init__(
        self,
        nroots: int = 1,
        basis: str = "sto-3g",
        unit: str = "Bohr",
        charge: int = 0,
    ) -> None:
        if nroots <= 0:
            raise ValueError(f"nroots must be positive, got {nroots}")

        self._energies = [None]*nroots
        self._gradients = [None]*nroots
        self._basis = basis
        self._unit = unit
        self._charge = int(charge)
        self._geom: MolecularGeometry | None = None
        self._cache = None  # for derived class specific in ram or disk caching  

    #Public API

    def _set_geom(self, geom: MolecularGeometry) -> None:
        self._geom = geom
        self._energies = [None]*len(self._energies)
        self._gradients = [None]*len(self._gradients)
        pass
    
    def get_geometry(self) -> MolecularGeometry:
        if self._geom is None:
            raise ValueError("Geometry has not been set.")
        return self._geom

    def get_nroots(self) -> int:
        """Return the number of roots (electronic states) computed by this strategy."""
        return len(self._energies)

    def get_energy(self, root: int) -> float: # if energy is None, call the compute_energies method 
        """Return the energy of one state in Hartree.

        If the energies have not been cached yet, compute and cache them first.
        """
        if self._energies[root] is None:
            self._compute_energies()
        return float(self._energies[root])

    def get_gradient(self, root: int) -> np.ndarray:
        """Return the nuclear gradient for *root*.

        If the gradient has not been cached yet, compute and cache it first.

        Returns
        -------
        np.ndarray
            Shape ``(natoms, 3)`` in Hartree/Bohr.
        """
        if self._energies[root] is None:
            self._compute_energies()
        if self._gradients[root] is None:
            self._compute_gradient(root)
        return self._gradients[root]

    def time_overlap(self, right: ElectronicStructureStrategy) -> np.ndarray: #error handling wrapper for _compute_time_overlap
        """Return the state-overlap matrix ``<psi_i(left)|psi_j(right)>``.

        Returns
        -------
        np.ndarray
            Shape ``(nroots, nroots)``.
        """
        if type(self) is not type(right):
                raise TypeError(
                    f"time overlap requires same strategy type, got "
                    f"{type(self).__name__} and {type(right).__name__}"
                )
        return self._compute_time_overlap(right)
    
    def get_nac_vectors(self) -> np.ndarray:
        """Return all NAC vectors ``d_{ij}`` between states.

        Returns
        -------
        np.ndarray
            Shape ``(nroots, nroots, natoms, 3)`` in 1/Bohr.
            Only off-diagonal elements are meaningful.
        """
        if self._energies[0] is None:
            self._compute_energies()
        return self._compute_nac_vectors()

    # ------------------------------------------------------------------
    #  Core computation
    # ------------------------------------------------------------------

    @abstractmethod
    def _compute_energies(self) -> None: 
        """Compute and cache the total energies (Hartree) for all roots."""
        

    def _compute_gradient(self, root: int) -> None:
        """Compute and cache the nuclear gradient for *root*.

        The computed gradient must be stored in ``self._gradients[root]`` with
        shape ``(natoms, 3)`` in Hartree/Bohr.
        """
        if self._energies[root] is None:
            raise ValueError(f"Energy for root {root} has not been computed yet.")

        raise NotImplementedError(
            "This backend does not provide nuclear gradients."
        )

    def _compute_time_overlap(self, other: "ElectronicStructureStrategy") -> np.ndarray:
        """Compute state time-overlap matrix with another snapshot.

        Returns
        -------
        np.ndarray
            Shape ``(nroots, nroots)``. Element ``S[i, j]`` is
            ``<psi_i(self) | psi_j(other)>``, where ``self`` is the left/previous
            electronic-structure snapshot and ``other`` is the right/current
            snapshot. Dimensionless.

        Notes
        -----
        Implementations should return the matrix and should not mutate base-class
        caches directly. Phase/state tracking, if required, should be handled by
        the backend or adapter before returning.
        """
        raise NotImplementedError(
            "This backend does not implement time-overlap-based NACs."
        )


    def _compute_nac_vectors(self) -> np.ndarray:
        """Return all NAC vectors ``d_{ij}`` between states.

        Returns
        -------
        np.ndarray
            Shape ``(nroots, nroots, natoms, 3)`` in 1/Bohr.
            Only off-diagonal elements are meaningful.

        Raises
        ------
        NotImplementedError
            If the backend does not support explicit NAC vectors.
        """
        if self._energies[0] is None:
            raise ValueError("Energies must be computed before NAC vectors.")   
        raise NotImplementedError(
            "This backend does not provide explicit NAC vectors; "
            "use time-overlap-based NACs instead."
        )
