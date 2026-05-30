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
from typing import Any, List, Optional, Tuple, Union, Sequence
import numpy as np
from pyscf import fci, gto, mcscf, scf
from libra_py.packages.pyscf.interfaces import ElectronicStructureStrategy, MolecularGeometry


class CASSCF(ElectronicStructureStrategy):
    """PySCF-based CASSCF backend for the universal ES interface."""

    def __init__(
        self,
        norbcas: int = 0,
        nelecas: int = 0,
        nroots: int = 1,
        basis: str = "sto-3g",
        unit: str = "Angstrom",
        charge: int = 0,
        cas_list: Optional[List[int]] = None,
    ) -> None:
        super().__init__(
            nroots=nroots,
            basis=basis,
            unit=unit,
            charge=int(charge),
        )

        self._norbcas: int = norbcas
        self._nelecas: Union[int, Tuple[int, int]] = nelecas
        self._cas_list: Optional[List[int]] = cas_list

        self.mol: Optional[Any] = None
        self.mf: Optional[Any] = None
        self.mc: Optional[Any] = None

    def set_geom(self, geom: MolecularGeometry) -> None:
        super().set_geom(geom)
        nroots = self.get_nroots()
        self._energies = [None] * nroots
        self._gradients = [None] * nroots

        self.mol = gto.M(
            atom=";".join(
                f"{label} {coord[0]} {coord[1]} {coord[2]}"
                for label, coord in zip(geom.atom_labels, geom.coords_angstrom)
            ),
            basis=self._basis,
            unit=self._unit,
            charge=self._charge,
            spin=0,
        )

        self.mf = None
        self.mc = None

    def run_hf(self) -> None:
        mf = scf.RHF(self.mol)
        mf.verbose = 0
        mf.kernel()

        self.mf = mf

    def compute_energies(self) -> None:
        nroots = self.get_nroots()
        if self.mc is None:
            mc = mcscf.CASSCF(self.mf, self._norbcas, self._nelecas)
            mc.fcisolver = fci.direct_spin0.FCI(self.mol)
            mc.fcisolver.nroots = nroots

            if nroots > 1:
                mc = mc.state_average_([1.0 / nroots] * nroots) # hardcoded equal weights for now;

            mo_coeff = self.mf.mo_coeff
            if self._cas_list is not None:
                mo_coeff = mcscf.sort_mo(mc, mo_coeff, self._cas_list)

            mc.kernel(mo_coeff)

            self.mc = mc

        e_states: Optional[Sequence[float]] = getattr(self.mc, "e_states", None)
        if e_states is not None:
            e_arr = np.asarray(e_states, dtype=float)
            self._energies = [float(e) for e in e_arr[:nroots]]
            return

        self._energies[0] = float(self.mc.e_tot)

    def compute_energy(self, root: int) -> float:
        return self.get_energy(root)

    def compute_gradient(self, root: int) -> None:
        if self._energies[root] is None:
            self.compute_energies()
        e_states: Optional[Sequence[float]] = getattr(self.mc, "e_states", None)
        if e_states is not None:
            self._gradients[root] = np.asarray(self.mc.nuc_grad_method(state=root).kernel())
            return
        self._gradients[root] = np.asarray(self.mc.nuc_grad_method().kernel())

    def _compute_time_overlap(self, right: "CASSCF") -> np.ndarray:
        nroots = self.get_nroots()

        left_roots = getattr(self.mc, "ci", None)
        right_roots = getattr(right.mc, "ci", None)

        if not isinstance(left_roots, (list, tuple)):
            left_roots = [left_roots]
        if not isinstance(right_roots, (list, tuple)):
            right_roots = [right_roots]

        left_roots = [np.asarray(v) for v in left_roots[:nroots]]
        right_roots = [np.asarray(v) for v in right_roots[:nroots]]

        ao_overlap = gto.intor_cross("int1e_ovlp", self.mol, right.mol)
        mo_left_act = np.asarray(self.mc.mo_coeff)[:, self.mc.ncore : self.mc.ncore + self.mc.ncas]
        mo_right_act = np.asarray(right.mc.mo_coeff)[:, right.mc.ncore : right.mc.ncore + right.mc.ncas]
        s12_mo = mo_left_act.T @ ao_overlap @ mo_right_act

        overlap = np.zeros((nroots, nroots), dtype=float)
        nelecas = self.mc.nelecas
        for i in range(nroots):
            for j in range(nroots):
                overlap[i, j] = fci.addons.overlap(
                    left_roots[i],
                    right_roots[j],
                    self.mc.ncas,
                    nelecas,
                    s=s12_mo,
                )

        overlap = np.asarray(np.real_if_close(overlap))

        for i in range(nroots):
            if overlap[i, i] < 0:
                overlap[i, :] = -overlap[i, :]

        return overlap

    def _compute_nac_vectors(self, use_etfs: bool = True) -> np.ndarray:
        nstates = self.get_nroots()
        natm = int(self.mol.natm)

        mc_nacs = self.mc.nac_method()
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
