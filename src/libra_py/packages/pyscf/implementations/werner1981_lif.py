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
.. module:: pyscf.implementations.werner1981_lif
   :platform: Unix, Windows
   :synopsis: Reference LiF benchmark backend from Werner and Meyer (1981).
.. moduleauthor::
       Jieyang Gu <jieyanggu792@gmail.com>

This implementation provides a lightweight two-state reference model for the
LiF avoided crossing using the tabulated adiabatic energies from Table II and
nonadiabatic couplings from Table V of:

    H.-J. Werner and W. Meyer, J. Chem. Phys. 74, 5802 (1981)

Important: the paper tabulates the bond length in bohr (a.u.). Incoming
geometries are expressed in Angstrom via ``MolecularGeometry`` and are
converted internally before interpolation.
"""

from __future__ import annotations
from typing import Any
from typing import Any

import numpy as np

from ..interfaces import ElectronicStructureStrategy, MolecularGeometry


class Werner1981LiF(ElectronicStructureStrategy):
    """Two-state LiF benchmark strategy derived from Werner--Meyer 1981."""

    _BOHR_TO_ANG = 0.529177210903

    _R_BOHR = np.array(
        [
            2.0,
            2.2,
            2.4,
            2.6,
            2.8,
            3.0,
            3.2,
            3.5,
            4.0,
            5.0,
            6.0,
            7.0,
            8.0,
            9.0,
            10.0,
            11.0,
            12.0,
            13.0,
            13.3,
            13.4,
            14.0,
            15.0,
            20.0,
        ],
        dtype=float,
    )

    _E1_HARTREE = np.array(
        [
            -106.840315,
            -106.953123,
            -107.015849,
            -107.048481,
            -107.063090,
            -107.066935,
            -107.064418,
            -107.054238,
            -107.030834,
            -106.985211,
            -106.950284,
            -106.924411,
            -106.904927,
            -106.889927,
            -106.878116,
            -106.868592,
            -106.860753,
            -106.854347,
            -106.853027,
            -106.852783,
            -106.852415,
            -106.852386,
            -106.852417,
        ],
        dtype=float,
    )

    _E2_HARTREE = np.array(
        [
            -106.566486,
            -106.669433,
            -106.731402,
            -106.768743,
            -106.791343,
            -106.804870,
            -106.813689,
            -106.821097,
            -106.827183,
            -106.834217,
            -106.840516,
            -106.845424,
            -106.848548,
            -106.850308,
            -106.851233,
            -106.851707,
            -106.851940,
            -106.851882,
            -106.851451,
            -106.851129,
            -106.848270,
            -106.843501,
            -106.826751,
        ],
        dtype=float,
    )

    _NAC_R_BOHR = np.array(
        [
            9.025,
            10.025,
            11.025,
            12.025,
            12.525,
            13.025,
            13.125,
            13.225,
            13.275,
            13.375,
            13.425,
            13.525,
            13.725,
            14.025,
            14.525,
            16.025,
        ],
        dtype=float,
    )

    _NAC12_INV_BOHR = np.array(
        [
            -0.005,
            -0.003,
            0.011,
            0.065,
            0.183,
            0.872,
            1.265,
            1.694,
            1.829,
            1.745,
            1.557,
            1.118,
            0.535,
            0.219,
            0.078,
            0.012,
        ],
        dtype=float,
    )

    _ENERGIES_HARTREE = np.vstack((_E1_HARTREE, _E2_HARTREE))
    _DEDR_HARTREE_PER_BOHR = np.vstack(
        (
            np.gradient(_E1_HARTREE, _R_BOHR, edge_order=2),
            np.gradient(_E2_HARTREE, _R_BOHR, edge_order=2),
        )
    )

    def __init__(self, mol: Any = None, nroots: int = 2, basis: str = "sto-3g", unit: str = "Angstrom", charge: int = 0) -> None:
        super().__init__(mol=mol, nroots=nroots, basis=basis, unit=unit, charge=charge)
        self._distance_bohr: float | None = None
        self._bond_unit: np.ndarray | None = None
        self._prev_distance_bohr: float | None = None

    @property
    def nstates(self) -> int:
        return 2

    @property
    def has_nac_vectors(self) -> bool:
        return True

    def _validate_geometry(self, geom: MolecularGeometry) -> tuple[float, np.ndarray]:
        coords = np.asarray(geom.coords_angstrom, dtype=float)
        if coords.shape != (2, 3):
            raise ValueError(
                "Werner1981LiF expects a diatomic geometry with shape (2, 3)."
            )

        labels = list(geom.atom_labels)
        if len(labels) != 2 or set(labels) != {"Li", "F"}:
            raise ValueError(
                "Werner1981LiF is parameterized only for LiF geometries."
            )

        li_idx = labels.index("Li")
        f_idx = labels.index("F")

        # Define the bond vector pointing from F to Li
        # Hardcode Z-axis only logic: require X and Y coordinates to be near zero
        for idx in (li_idx, f_idx):
            if abs(coords[idx, 0]) > 1e-6 or abs(coords[idx, 1]) > 1e-6:
                raise ValueError("Werner1981LiF hardcoded mode requires atoms to lie exactly on the Z-axis (x=0, y=0).")

        # Evaluate distance based solely on the Z coordinate
        bond_len_ang = abs(coords[li_idx, 2] - coords[f_idx, 2])
        if bond_len_ang <= 0.0:
            raise ValueError("LiF bond length must be positive.")

        bond_len_bohr = float(bond_len_ang) / self._BOHR_TO_ANG
        # We don't bother returning a proper 3D bond_unit anymore, since NAC is hardcoded Z
        return bond_len_bohr, np.array([0.0, 0.0, 1.0])

    @classmethod
    def _interp_energy(cls, root: int, r_bohr: float) -> float:
        r_min = float(cls._R_BOHR[0])
        r_max = float(cls._R_BOHR[-1])
        if r_bohr < r_min or r_bohr > r_max:
            raise ValueError(
                f"Bond length {r_bohr:.6f} bohr is outside the Table II range "
                f"[{r_min:.3f}, {r_max:.3f}] bohr."
            )
        return float(np.interp(r_bohr, cls._R_BOHR, cls._ENERGIES_HARTREE[root]))

    @classmethod
    def _interp_dedr(cls, root: int, r_bohr: float) -> float:
        r_min = float(cls._R_BOHR[0])
        r_max = float(cls._R_BOHR[-1])
        if r_bohr < r_min or r_bohr > r_max:
            raise ValueError(
                f"Bond length {r_bohr:.6f} bohr is outside the Table II range "
                f"[{r_min:.3f}, {r_max:.3f}] bohr."
            )
        return float(
            np.interp(r_bohr, cls._R_BOHR, cls._DEDR_HARTREE_PER_BOHR[root])
        )

    @classmethod
    def _interp_nac_scalar(cls, r_bohr: float) -> float:
        return float(
            np.interp(
                r_bohr,
                cls._NAC_R_BOHR,
                cls._NAC12_INV_BOHR,
                left=0.0,
                right=0.0,
            )
        )

    def _require_geometry(self) -> tuple[float, np.ndarray]:
        if self._distance_bohr is None or self._bond_unit is None:
            raise ValueError("Geometry must be set before requesting LiF data.")
        return self._distance_bohr, self._bond_unit

    def save_cache(self) -> None:
        self._prev_distance_bohr = self._distance_bohr

    def run_hf(self) -> None:
        if self._geom is None:
            raise ValueError("Geometry must be set before running.")
        self._distance_bohr, self._bond_unit = self._validate_geometry(self._geom)

    def compute_energy(self, root: int) -> float:
        if root < 0 or root >= self.nstates:
            raise IndexError(
                f"Requested root {root}, but only {self.nstates} roots are available."
            )
        r_bohr, _ = self._require_geometry()
        return self._interp_energy(root, r_bohr)

    def compute_nac_vectors(self, **kwargs: Any) -> np.ndarray:
        r_bohr, bond_unit = self._require_geometry()
        nac12 = self._interp_nac_scalar(r_bohr)

        nac = np.zeros((self.nstates, self.nstates, 2, 3), dtype=float)
        labels = list(self._geom.atom_labels)
        li_idx = labels.index("Li")
        f_idx = labels.index("F")

        # Hardcode the NAC to be strictly along the Z-axis (index 2)
        nac[0, 1, li_idx, 2] = nac12
        nac[0, 1, f_idx, 2] = -nac12
        nac[1, 0, li_idx, 2] = -nac12
        nac[1, 0, f_idx, 2] = nac12
        return nac

