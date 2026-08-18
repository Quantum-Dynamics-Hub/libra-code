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
   CP2K, or other quantum-chemistry codes.

   One ES_Strategy instance represents one electronic-structure calculation
   at one nuclear geometry. Consecutive calculations are represented by
   separate ES_Strategy instances.

.. moduleauthor::
       Jieyang Gu <jieyanggu792@gmail.com>
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Literal, Sequence

import numpy as np

# ---------------------------------------------------------------------------
# Data containers
# ---------------------------------------------------------------------------

@dataclass
class MolecularGeometry:
    """Molecular geometry in atomic units."""

    atom_labels: tuple[str, ...]
    coords_bohr: np.ndarray  # shape: (natoms, 3)


@dataclass
class ES_Request:
    """Quantities requested at one geometry snapshot."""

    n_singlets: int = 1
    n_triplets: int = 0

    H_soc: bool = False

    # None  -> no gradient
    # int   -> gradient for one state
    # "all" -> gradients for all states
    gradient_state: int | Literal["all"] | None = None

    # None -> no Hessian
    # int  -> Hessian for one state
    hessian_state: int | None = None

    nacv: bool = False
    time_overlap: bool = False

    @property
    def n_total(self) -> int:
        """Total number of electronic states.

        Currently each triplet manifold is counted as one state.
        """
        return self.n_singlets + self.n_triplets


@dataclass
class ES_Result:
    """Computed quantities at one geometry snapshot.

    Units:
        energies: Hartree
        gradients: Hartree / Bohr
        Hessians: Hartree / Bohr^2
        NAC vectors: Bohr^-1
    """
    H_el: np.ndarray | None = None  # shape: (n_total, n_total)

    H_soc: np.ndarray | None = None  # shape: (n_total, n_total)

    gradients: list[np.ndarray | None] | None = None  # one entry per electronic state, each non-None entry has shape (natoms, 3)

    hessians: list[np.ndarray | None] | None = None  # one entry per electronic state, each non-None entry has shape (3*natoms, 3*natoms)

    nac_vectors: np.ndarray | None = None  # shape: (n_total, n_total, natoms, 3)

    time_overlap: np.ndarray | None = None  # shape: (n_total, n_total)


# ---------------------------------------------------------------------------
# Backend ABC
# ---------------------------------------------------------------------------

class ES_Strategy(ABC):
    """Abstract electronic-structure strategy.

    One ES_Strategy instance represents one geometry and its associated
    electronic-structure calculation.
    """

    @abstractmethod
    def snapshot_state(self) -> None: # save the state of the current calculation to previous, overwriting whatever was there before. 
        """Save the current state of the calculation to a previous-state snapshot."""
        raise NotImplementedError   

    @abstractmethod
    def get_state(self) -> object: #object is a dataclass or dict containing the state, or a name or a path to a folder for disk based cache.
        """Return the current state of the calculation."""
        raise NotImplementedError

    @abstractmethod
    def get_previous_state(self) -> object | None: #object is a dataclass or dict containing the state, or a name or a path to a folder for disk based cache.
        """Return the cached previous-state snapshot, path, or handle if present."""
        raise NotImplementedError

    # -----------------------------------------------------------------------
    # Geometry
    # -----------------------------------------------------------------------

    @abstractmethod
    def set_geom(self, geom: MolecularGeometry) -> None:
        """Set the nuclear geometry for this calculation."""
        raise NotImplementedError

    @abstractmethod
    def get_geom(self) -> MolecularGeometry:
        """Return the geometry represented by this strategy."""
        raise NotImplementedError

    # -----------------------------------------------------------------------
    # Required electronic-structure calculation
    # -----------------------------------------------------------------------

    @abstractmethod
    def compute_H_el( self ) -> np.ndarray:
        """Compute adiabatic electronic energies."""
        raise NotImplementedError

    # -----------------------------------------------------------------------
    # Optional electronic-structure quantities
    # -----------------------------------------------------------------------

    def compute_H_soc(self) -> np.ndarray:
        """Compute the spin-orbit Hamiltonian.

        Returns
        -------
        np.ndarray
            Shape ``(n_total, n_total)`` in Hartree.
        """
        raise NotImplementedError(
            f"{type(self).__name__}: "
            "override compute_H_soc or do not request H_soc."
        )

    def compute_gradient(self, roots: Sequence[int]) -> list[np.ndarray]:
        """Compute nuclear gradients for the given electronic states.

        Backends receive the whole set at once so they can build whatever the
        states share -- integrals, response-equation operators -- a single time.

        Returns
        -------
        list[np.ndarray]
            One entry per element of ``roots``, in that order, each of shape
            ``(natoms, 3)`` in Hartree/Bohr.
        """
        raise NotImplementedError(
            f"{type(self).__name__}: "
            "override compute_gradient or do not request gradients."
        )

    def compute_hessian(self, root: int = 0) -> np.ndarray:
        """Compute the nuclear Hessian for one electronic state.

        Returns
        -------
        np.ndarray
            Shape ``(3N, 3N)`` in Hartree/Bohr^2.
        """
        raise NotImplementedError(
            f"{type(self).__name__}: "
            "override compute_hessian or do not request Hessians."
        )

    def compute_nac_vectors(self) -> np.ndarray:
        """Compute nonadiabatic coupling vectors.

        Returns
        -------
        np.ndarray
            Shape ``(n_total, n_total, natoms, 3)`` in Bohr^-1.
        """
        raise NotImplementedError(
            f"{type(self).__name__}: "
            "override compute_nac_vectors or do not request NAC vectors."
        )

    
    def compute_time_overlap( self, state1: object, state2: object ) -> np.ndarray:
        """Compute the time-overlap matrix between two states.
        np.ndarray
            Shape ``(n_total, n_total)``.
        """
        raise NotImplementedError(
            f"{type(self).__name__}: "
            "override compute_time_overlap or do not request time overlaps."
        )

    # -----------------------------------------------------------------------
    # Sequencing contract
    # -----------------------------------------------------------------------

    def compute_result( self, geom: MolecularGeometry, request: ES_Request ) -> ES_Result:

        result = ES_Result()


        """Compute all requested quantities at one geometry and enforce a sequence of calculations."""
        if request.n_triplets != 0:
            raise NotImplementedError( "Triplet states are not supported yet." )

        n_total = request.n_total
        natoms = len(geom.coords_bohr)
            
        self.set_geom(geom)

        result.H_el = np.asarray( self.compute_H_el(), dtype=np.float64 )

        if request.H_soc:
            result.H_soc = np.asarray( self.compute_H_soc() )

        if request.gradient_state is not None:
            if request.gradient_state == "all":
                roots = list(range(n_total))
            else:
                roots = [int(request.gradient_state)]
            for root in roots:
                if not 0 <= root < n_total:
                    raise ValueError(
                        f"gradient_state={root} is outside the valid range [0, {n_total})."
                    )

            # Scatter into a slot-per-state layout; the backend only ever sees
            # the roots that were asked for, in the order it was asked for them.
            result.gradients = [None] * n_total
            for root, grad in zip(roots, self.compute_gradient(roots)):
                result.gradients[root] = grad

        if request.hessian_state is not None:

            root = request.hessian_state

            result.hessians = [None] * n_total

            hessian = np.asarray( self.compute_hessian(root), dtype=np.float64 )
            result.hessians[root] = hessian

        if request.nacv is True:

            result.nac_vectors = np.asarray( self.compute_nac_vectors(), dtype=np.float64 )

        if request.time_overlap is True:
            if self.get_previous_state() is  not None: # if none, then this is the first geometry, and there is no previous state to compute time-overlap with.
                state1 = self.get_state()
                state2 = self.get_previous_state()
                result.time_overlap = np.asarray(
                    self.compute_time_overlap(state1, state2),
                    dtype=np.float64,
                )   
            else:
                result.time_overlap = None 


        self.snapshot_state() #copy state to previous after the calculation is done, so that the next geometry can use it for time-overlap and NACV calculations.

        return result



# Backward-compatible alias used by older PySCF package imports.
ElectronicStructureStrategy = ES_Strategy
