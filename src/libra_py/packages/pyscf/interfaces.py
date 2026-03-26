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
   :synopsis: Core interface definitions for PySCF-backed electronic structure strategies.
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
    atom_labels: list[str]
    coords_angstrom: np.ndarray


class ElectronicStructureStrategy(ABC):
    """Base interface for electronic structure backends.

    Implementations are allowed to store any backend-specific state internally.

    The interface is intentionally minimal: it only describes the operations
    required by a consumer (e.g. NAMD) and does not prescribe how a backend
    achieves those operations.
    """

    @abstractmethod
    def set_geom_and_run_hf(self, geom: MolecularGeometry) -> None: 
        """Set the current geometry and run HF on it (storing results)."""

    @abstractmethod
    def compute_energy(self, root: int) -> float:     #singlets only
        """Compute the energy for a given root."""

    @abstractmethod
    def compute_gradient(self, root: int) -> np.ndarray:
        """Compute the nuclear gradient for a given root."""

    @abstractmethod
    def time_overlap_matrix(self, nroots: int) -> np.ndarray:
        """Compute the time-overlap matrix for all roots."""
