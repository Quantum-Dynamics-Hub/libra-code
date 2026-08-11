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
.. module:: lvc_strategy
   :platform: Unix, Windows
   :synopsis: Linear vibronic coupling (LVC) implementation of ES_Strategy.

   Implements the ES_Strategy interface for systems described by a linear
   (or quadratic) vibronic coupling Hamiltonian parametrised from a reference
   electronic-structure calculation.

   Supported quantities
   --------------------
   - ``compute_H_el``        adiabatic energies (always available)
   - ``compute_H_soc``       spin-orbit Hamiltonian (requires SOC block in template)
   - ``compute_gradient``    nuclear gradient for one or all states
   - ``compute_nac_vectors`` non-adiabatic coupling vectors
   - ``compute_time_overlap`` wavefunction overlap between consecutive steps

   Not implemented
   ---------------
   - ``compute_hessian``     LVC has no analytic second derivatives; the base-class
                             ``NotImplementedError`` is intentionally preserved.

.. moduleauthor:: Jieyang Gu <jieyanggu792@gmail.com>
"""

from __future__ import annotations

import copy
import datetime
import os
from dataclasses import dataclass, field
from typing import Literal

import numpy as np

# ---------------------------------------------------------------------------
# Base-class interface
# ---------------------------------------------------------------------------
from pyscf.interfaces import (
    ES_Request,
    ES_Result,
    ES_Strategy,
    MolecularGeometry,
)

# ---------------------------------------------------------------------------
# LVC parametrisation helpers (internal – not re-exported)
# ---------------------------------------------------------------------------
from constants import U_TO_AMU
from kabsch import kabsch_w as _kabsch, kabsch_w_with_deriv as _kabsch_with_deriv
from utils import expand_path, phase_correction, readfile


# ===========================================================================
# State snapshot
# ===========================================================================

@dataclass
class LVCState:
    """Minimal data needed to restart or form a time-overlap from one step.

    Attributes
    ----------
    coords_bohr:
        Nuclear coordinates at this step, shape ``(natoms, 3)``.
    U:
        Adiabatic-to-diabatic transformation matrix at this step,
        shape ``(nmstates, nmstates)``.  The time-overlap between two
        consecutive steps is ``U_prev.T @ U_curr``.
    H_el:
        Diagonal adiabatic energies at this step, shape ``(nmstates,)``.
    """

    coords_bohr: np.ndarray
    U: np.ndarray
    H_el: np.ndarray


# ===========================================================================
# LVC parametrisation  (reads template + V0 file)
# ===========================================================================

@dataclass
class _LVCParams:
    """Internal container for all LVC model parameters read from the template."""

    # geometry reference
    ref_coords: np.ndarray          # (natoms, 3)  Bohr
    masses: np.ndarray              # (natoms,)    AMU
    Msa: np.ndarray                 # (3*natoms,)  sqrt(m * U_TO_AMU)
    Om: np.ndarray                  # (nmodes,)    harmonic frequencies
    Km: np.ndarray                  # (nmodes, 3*natoms)  mass-weighted normal modes

    # per-multiplicity diabatic blocks
    h: dict       # {im: (n,n)}          on-diagonal (epsilon + eta)
    H_i: dict     # {im: (3N, n, n)}     linear coupling (kappa + lambda)
    G: dict       # {im: (n,n,3N,3N)}    quadratic coupling (gamma)
    gammas: bool  # whether gamma terms are nonzero

    # constant operators (all in nmstates basis)
    soc: np.ndarray       # (nmstates, nmstates)  real or complex
    dipole: np.ndarray    # (3, nmstates, nmstates)
    lambda_soc: np.ndarray | None  # (nmstates, nmstates, 3N) or None

    # state bookkeeping
    states: list[int]
    nmstates: int


# ===========================================================================
# Main class
# ===========================================================================

class LVCStrategy(ES_Strategy):
    """LVC implementation of the ES_Strategy interface.

    One instance represents **all** consecutive geometry steps of one
    trajectory.  Call ``compute_result`` at each step; the class manages
    the previous-step snapshot internally.

    Parameters
    ----------
    template_filename:
        Path to the LVC template file (default ``"LVC.template"``).
    do_kabsch:
        Whether to Kabsch-align the input geometry to the reference before
        projecting onto normal modes (default ``False``).
    diagonalize:
        Diagonalise the diabatic Hamiltonian to obtain adiabatic states
        (default ``True``).  Set to ``False`` to stay in the diabatic basis.
    """

    def __init__(
        self,
        template_filename: str = "LVC.template",
        do_kabsch: bool = False,
        diagonalize: bool = True,
    ) -> None:
        self._params: _LVCParams = _read_template(template_filename)
        self._do_kabsch = do_kabsch
        self._diagonalize = diagonalize

        # Current geometry and cached computation outputs
        self._geom: MolecularGeometry | None = None
        self._U: np.ndarray | None = None          # set after each _run()
        self._cache: _RunCache | None = None

        # State snapshots
        self._current_state: LVCState | None = None
        self._previous_state: LVCState | None = None

    # -----------------------------------------------------------------------
    # ES_Strategy: state management
    # -----------------------------------------------------------------------

    def snapshot_state(self) -> None:
        if self._current_state is not None:
            self._previous_state = copy.deepcopy(self._current_state)

    def get_state(self) -> LVCState | None:
        return self._current_state

    def get_previous_state(self) -> LVCState | None:
        return self._previous_state

    # -----------------------------------------------------------------------
    # ES_Strategy: geometry
    # -----------------------------------------------------------------------

    def set_geom(self, geom: MolecularGeometry) -> None:
        self._geom = geom
        self._cache = None  # invalidate on geometry change

    def get_geom(self) -> MolecularGeometry:
        if self._geom is None:
            raise RuntimeError("No geometry has been set yet.")
        return self._geom

    # -----------------------------------------------------------------------
    # ES_Strategy: required quantity
    # -----------------------------------------------------------------------

    def compute_H_el(self) -> np.ndarray:
        """Adiabatic energies, shape ``(n_total,)`` in Hartree."""
        return self._run().H_el

    # -----------------------------------------------------------------------
    # ES_Strategy: optional quantities
    # -----------------------------------------------------------------------

    def compute_H_soc(self) -> np.ndarray:
        """Spin-orbit Hamiltonian, shape ``(n_total, n_total)``."""
        cache = self._run()
        if cache.H_soc is None:
            raise NotImplementedError(
                "LVCStrategy: SOC matrix was not computed. "
                "Add a 'SOC' block to the LVC template."
            )
        return cache.H_soc

    def compute_gradient(self, root: int = 0) -> np.ndarray:
        """Nuclear gradient for state ``root``, shape ``(natoms, 3)`` in Ha/Bohr."""
        cache = self._run()
        if cache.gradients is None or cache.gradients[root] is None:
            raise NotImplementedError(
                f"LVCStrategy: gradient for root {root} was not computed."
            )
        return cache.gradients[root]

    def compute_nac_vectors(self) -> np.ndarray:
        """NAC vectors, shape ``(n_total, n_total, natoms, 3)`` in Bohr^{-1}."""
        cache = self._run()
        if cache.nac_vectors is None:
            raise NotImplementedError(
                "LVCStrategy: NAC vectors were not computed."
            )
        return cache.nac_vectors

    def compute_time_overlap(
        self,
        state1: LVCState,
        state2: LVCState,
    ) -> np.ndarray:
        """Wavefunction time-overlap between two consecutive steps.

        The overlap matrix is defined as

        .. math::

            S_{ij}(t, t{+}\\Delta t)
            = \\langle \\psi_i(t) \\mid \\psi_j(t{+}\\Delta t) \\rangle
            = U_{\\text{prev}}^\\top \\, U_{\\text{curr}}

        where ``state2`` is the *previous* step and ``state1`` is the *current*
        step, matching the calling convention in ``ES_Strategy.compute_result``.

        Parameters
        ----------
        state1:
            Current-step ``LVCState`` (returned by ``get_state()``).
        state2:
            Previous-step ``LVCState`` (returned by ``get_previous_state()``).

        Returns
        -------
        np.ndarray
            Shape ``(n_total, n_total)``, real-valued.
        """
        # Mimics the exact expression used internally:
        #   overlap = Uold.T @ self._U
        return state2.U.T @ state1.U

    # -----------------------------------------------------------------------
    # Internal: run LVC once per geometry, cache outputs
    # -----------------------------------------------------------------------

    def _run(self) -> "_RunCache":
        """Execute the LVC Hamiltonian at the current geometry (cached)."""
        if self._cache is not None:
            return self._cache

        if self._geom is None:
            raise RuntimeError("Call set_geom() before any compute_* method.")

        p = self._params
        coords = np.array(self._geom.coords_bohr, dtype=float)
        natom = coords.shape[0]
        r3N = 3 * natom
        nmstates = p.nmstates
        states = p.states
        req_nmstates = nmstates  # expose all states

        # ------------------------------------------------------------------
        # Optional Kabsch alignment
        # ------------------------------------------------------------------
        coords_ref_basis = coords
        Trot = np.eye(3)
        if self._do_kabsch:
            Trot, com_ref, com_coords = _kabsch(p.ref_coords, 