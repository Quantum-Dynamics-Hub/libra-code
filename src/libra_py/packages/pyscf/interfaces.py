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

import copy as pycopy

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

        The default implementation is the numerical fallback at the bottom of
        this module: it central-differences the wavefunction overlap with
        respect to every nuclear coordinate, so a backend that can only produce
        ``compute_time_overlap`` gets NAC vectors for free.  Backends with an
        analytic coupling (CASSCF, TDDFT response) should override this -- the
        fallback costs ``6 * natoms`` extra electronic-structure calculations
        per geometry.

        Returns
        -------
        np.ndarray
            Shape ``(n_total, n_total, natoms, 3)`` in Bohr^-1.
        """
        return numerical_nac_vectors(self)

    
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



# ---------------------------------------------------------------------------
# Fallbacks
# ---------------------------------------------------------------------------
#
# Generic implementations of the optional quantities, written entirely against
# the ES_Strategy interface.  They let a backend that implements only the
# cheap, always-available pieces still answer the full request; a backend with
# an analytic route overrides the corresponding method and never reaches here.

# Central-difference displacement for the numerical NAC vectors, in Bohr.
#
# Hard-coded rather than exposed as a parameter: the useful window is narrow
# and roughly the same for every backend.  The finite-difference error is the
# sum of a truncation term ~ h^2 * (third derivative) and a noise term
# ~ eps_S / (2h), where eps_S is how reproducibly the backend converges its
# wavefunction overlaps (~1e-7..1e-6 for a default-threshold CASSCF/CASCI).
# At h = 1e-3 Bohr both sit near 1e-4 Bohr^-1, which is well below the
# ~1e-1..1e+1 Bohr^-1 scale of a coupling that matters.
NUMERICAL_NACV_STEP_BOHR = 1.0e-3


def _clone_strategy(strategy: ES_Strategy) -> ES_Strategy:
    """Return an independent copy of ``strategy``.

    The displaced-geometry calculations must not touch the caller's strategy:
    it is mid-``compute_result`` at the reference geometry, and set_geom would
    overwrite both its current and its previous state.  Same copy/clone/deepcopy
    ladder the Libra adapter uses when it hands one strategy template to many
    trajectories.
    """
    for name in ("copy", "clone"):
        cloner = getattr(strategy, name, None)
        if callable(cloner):
            return cloner()
    return pycopy.deepcopy(strategy)


def numerical_nac_vectors(
    strategy: ES_Strategy,
    step: float = NUMERICAL_NACV_STEP_BOHR,
) -> np.ndarray:
    """NAC vectors by central-differencing the wavefunction overlap.

    The coupling is the derivative of an overlap::

        d_ij^(A,x) = < psi_i(R) | d/dR_Ax | psi_j(R) >
                   ~ [ <psi_i(R)|psi_j(R+h)> - <psi_i(R)|psi_j(R-h)> ] / (2h)

    with ``h = step`` along one Cartesian coordinate at a time, so the whole
    thing is built out of ``compute_time_overlap`` calls: the bra stays pinned
    at the reference geometry and only the ket moves.  That matches the
    argument order used by ``compute_result``, namely
    ``compute_time_overlap(state1, state2) -> <state2_i | state1_j>``, where
    state2 is the earlier -- here, the reference -- snapshot.

    Two properties of the exact d_ij are imposed rather than hoped for: it is
    antisymmetric for real wavefunctions, and its diagonal vanishes.  Averaging
    d against -d^T also cancels the leading symmetric part of the numerical
    error, so the symmetrization buys a little accuracy and is not only
    cosmetic.

    Relies on the backend's own overlap phase convention to keep the displaced
    roots gauge-aligned with the reference ones (CASSCF, for one, fixes the
    sign of each column so the diagonal comes out positive).  Without that the
    two displacements could disagree on the sign of a root and the difference
    would be meaningless.

    Requires a converged reference calculation to already be in place -- i.e.
    this runs after compute_H_el, which is what compute_result guarantees.

    Returns
    -------
    np.ndarray
        Shape ``(n_total, n_total, natoms, 3)`` in Bohr^-1.
    """
    geom = strategy.get_geom()
    coords = np.asarray(geom.coords_bohr, dtype=np.float64)
    natoms = coords.shape[0]

    if natoms == 0:
        raise ValueError("Numerical NAC vectors need at least one atom.")

    # Frozen copy of the reference wavefunction: it is the bra of every overlap
    # below, and the worker would otherwise mutate it out from under us.
    reference_state = pycopy.deepcopy(strategy.get_state())
    if reference_state is None:
        raise ValueError(
            f"{type(strategy).__name__}: numerical NAC vectors need a converged "
            "reference calculation; call compute_H_el before compute_nac_vectors."
        )

    # One worker for all 6*natoms displaced calculations.  Reusing it keeps the
    # backend's own warm-start chain alive -- each displaced SCF/CASSCF starts
    # from the previous one -- which is both faster and steadier in the root
    # ordering across the scan.
    worker = _clone_strategy(strategy)

    nacv = None

    for atom in range(natoms):
        for xyz in range(3):

            overlaps = []

            for direction in (+1.0, -1.0):

                displaced_coords = coords.copy()
                displaced_coords[atom, xyz] += direction * step

                worker.set_geom(
                    MolecularGeometry(
                        atom_labels=geom.atom_labels,
                        coords_bohr=displaced_coords,
                    )
                )
                worker.compute_H_el()

                # <reference_i | displaced_j>
                overlaps.append(
                    np.asarray(
                        worker.compute_time_overlap(worker.get_state(), reference_state),
                        dtype=np.float64,
                    )
                )

            derivative = (overlaps[0] - overlaps[1]) / (2.0 * step)

            if nacv is None:
                n_total = derivative.shape[0]
                nacv = np.zeros((n_total, n_total, natoms, 3), dtype=np.float64)

            nacv[:, :, atom, xyz] = derivative

    # Enforce antisymmetry and a zero diagonal.
    nacv = 0.5 * (nacv - np.swapaxes(nacv, 0, 1))
    for state in range(nacv.shape[0]):
        nacv[state, state] = 0.0

    return nacv


# Backward-compatible alias used by older PySCF package imports.
ElectronicStructureStrategy = ES_Strategy
