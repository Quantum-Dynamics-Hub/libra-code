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
.. module:: pyscf.implementations.tddft
   :platform: Unix, Windows
   :synopsis: PySCF TDDFT/TDA implementation for Libra ES strategy.

   Class names are deliberately TDDFT-prefixed (``TDDFT``, ``TDDFT_States``) so
   callers can tell they are using linear-response TDDFT rather than CASSCF
   (``CASSCF``, ``CASSCF_States``).

.. moduleauthor::
       Jieyang Gu <jieyanggu792@gmail.com>

"""

from __future__ import annotations

import copy as pycopy
from dataclasses import dataclass
from typing import Any, List, Optional, Sequence

import numpy as np
from pyscf import dft, gto, scf, tdscf

from libra_py.packages.pyscf.interfaces import ES_Strategy, ES_Request, MolecularGeometry

# Electron masses per amu; Libra nuclear dynamics use atomic units throughout.
AMU = 1822.888486209


def _compute_ao_overlap(prev_mol: Any, curr_mol: Any) -> np.ndarray:
    return gto.intor_cross("int1e_ovlp", prev_mol, curr_mol)


@dataclass
class TDDFT_States:
    """Cached electronic-structure snapshot for one TDDFT geometry.

    Distinct from ``CASSCF_States`` (which stores ``mc``): here the excited-state
    object is a PySCF ``tdscf`` / TDA instance, plus the CIS amplitudes used for
    time overlaps.
    """

    mol: Optional[Any] = None
    mf: Optional[Any] = None
    td: Optional[Any] = None
    mo_coeff: Optional[np.ndarray] = None
    dm: Optional[np.ndarray] = None
    amplitudes: Optional[List[np.ndarray]] = None
    energies: Optional[np.ndarray] = None


class TDDFT(ES_Strategy):
    """PySCF-based TDDFT/TDA backend for the universal ES interface.

    States are ordered ``[S0, S1, ..., S_nexc]`` so ``nstates = nexc + 1``.
    Unlike ``CASSCF``, the ground state comes from DFT/SCF and the excited
    manifold from linear-response TDA/TDDFT.

    Notes
    -----
    * Mainline PySCF has no general analytic derivative couplings for TDDFT,
      so ``compute_nac_vectors`` is not implemented; use time overlaps.
    * State overlaps use the factorized (neglect-of-orbital-relaxation) CIS form.
    """

    def __init__(
        self,
        atom_labels: Optional[Sequence[str]] = None,
        symbols: Optional[Sequence[str]] = None,
        nexc: int = 1,
        basis: str = "sto-3g",
        xc: str = "pbe0",
        unit: str = "Bohr",
        charge: int = 0,
        spin: int = 0,
        grid_level: int = 3,
        conv_tol: float = 1e-9,
        use_tda: bool = True,
        phase_tol: float = 0.5,
        mol: Optional[Any] = None,
    ) -> None:
        labels = atom_labels if atom_labels is not None else symbols
        if labels is None and mol is None:
            raise ValueError("TDDFT requires atom_labels/symbols or a prebuilt mol.")

        self._atom_labels: tuple[str, ...] = (
            tuple(labels) if labels is not None else tuple(mol.atom_symbol(i) for i in range(mol.natm))
        )
        self._nexc: int = int(nexc)
        self._nstates: int = self._nexc + 1
        self._basis: str = basis
        self._xc: str = xc
        self._unit: str = unit
        self._charge: int = int(charge)
        self._spin: int = int(spin)
        self._grid_level: int = int(grid_level)
        self._conv_tol: float = float(conv_tol)
        self._use_tda: bool = bool(use_tda)
        self._phase_tol: float = float(phase_tol)

        self._mol: Optional[Any] = mol
        self._mf: Optional[Any] = None
        self._td: Optional[Any] = None
        self._geom: Optional[MolecularGeometry] = None
        self._request: Optional[ES_Request] = None
        self._ao_overlap: Optional[np.ndarray] = None
        self._state: TDDFT_States | None = None
        self._previous_state: TDDFT_States | None = None
        self._occ: Optional[np.ndarray] = None
        self._vir: Optional[np.ndarray] = None

    # ------------------------------------------------------------------ aliases
    # Notebook / older code often used these PySCFSource attribute names.

    @property
    def atom_labels(self) -> tuple[str, ...]:
        return self._atom_labels

    @property
    def symbols(self) -> list[str]:
        return list(self._atom_labels)

    @property
    def nexc(self) -> int:
        return self._nexc

    @property
    def nstates(self) -> int:
        return self._nstates

    @property
    def natoms(self) -> int:
        return len(self._atom_labels)

    @property
    def basis(self) -> str:
        return self._basis

    @property
    def xc(self) -> str:
        return self._xc

    @property
    def grid_level(self) -> int:
        return self._grid_level

    @property
    def conv_tol(self) -> float:
        return self._conv_tol

    def init_kwargs(self) -> dict[str, Any]:
        """Keyword arguments needed to construct an empty twin of this strategy."""
        return {
            "atom_labels": self._atom_labels,
            "nexc": self._nexc,
            "basis": self._basis,
            "xc": self._xc,
            "unit": self._unit,
            "charge": self._charge,
            "spin": self._spin,
            "grid_level": self._grid_level,
            "conv_tol": self._conv_tol,
            "use_tda": self._use_tda,
            "phase_tol": self._phase_tol,
        }

    def clone_empty(self) -> "TDDFT":
        """Return a fresh TDDFT with the same settings and no cached ES state."""
        return TDDFT(**self.init_kwargs())

    def reset(self) -> None:
        """Clear cached electronic-structure state (keep constructor settings)."""
        self._geom = None
        self._state = None
        self._previous_state = None
        self._mol = None
        self._mf = None
        self._td = None
        self._ao_overlap = None
        self._request = None
        self._occ = None
        self._vir = None

    def build_mol(self, coords_bohr: np.ndarray):
        """Build a PySCF ``Mole`` at ``coords_bohr`` without running SCF."""
        coords = np.asarray(coords_bohr, dtype=float)
        mol = gto.Mole()
        mol.atom = [
            [label, tuple(float(x) for x in xyz)]
            for label, xyz in zip(self._atom_labels, coords)
        ]
        mol.basis = self._basis
        mol.unit = self._unit
        mol.charge = self._charge
        mol.spin = self._spin
        mol.verbose = 0
        mol.build()
        return mol

    def masses_au(self, coords_bohr: np.ndarray) -> list[float]:
        """Length-3N nuclear masses in a.u., Libra flat dof ordering."""
        masses = self.build_mol(coords_bohr).atom_mass_list()
        return [float(mass) * AMU for mass in masses for _ in range(3)]

    # ------------------------------------------------------------------ geometry

    def set_geom(self, geom: MolecularGeometry) -> None:
        self._geom = geom
        self._previous_state = self.get_state()
        self._td = None
        self._ao_overlap = None

        coords = np.asarray(geom.coords_bohr, dtype=float)
        labels = tuple(geom.atom_labels)
        self._atom_labels = labels

        self._mol = gto.M(
            atom=";".join(
                f"{label} {coord[0]} {coord[1]} {coord[2]}"
                for label, coord in zip(labels, coords)
            ),
            basis=self._basis,
            unit=self._unit,
            charge=self._charge,
            spin=self._spin,
            verbose=0,
        )

        prev_state = self.get_previous_state()
        if prev_state is not None and prev_state.mol is not None:
            self._ao_overlap = _compute_ao_overlap(prev_state.mol, self._mol)

        self._mf = dft.RKS(self._mol, xc=self._xc)
        self._mf.grids.level = self._grid_level
        self._mf.conv_tol = self._conv_tol

        dm0 = None
        if prev_state is not None and prev_state.mol is not None and prev_state.dm is not None:
            try:
                dm0 = scf.addons.project_dm_nr2nr(prev_state.mol, prev_state.dm, self._mol)
            except Exception:
                dm0 = None

        self._mf.kernel(dm0=dm0)
        if not self._mf.converged:
            print("  WARNING: TDDFT SCF not converged")

        nocc = int((self._mf.mo_occ > 0).sum())
        nmo = self._mf.mo_coeff.shape[1]
        self._occ = np.arange(nocc)
        self._vir = np.arange(nocc, nmo)

        self._state = TDDFT_States(
            mol=self._mol,
            mf=self._mf,
            td=None,
            mo_coeff=np.asarray(self._mf.mo_coeff),
            dm=np.asarray(self._mf.make_rdm1()),
            amplitudes=None,
            energies=None,
        )

    def get_geom(self) -> MolecularGeometry:
        if self._geom is None:
            raise ValueError("Geometry has not been set.")
        return self._geom

    def get_state(self) -> Optional[TDDFT_States]:
        return self._state

    def get_previous_state(self) -> Optional[TDDFT_States]:
        return self._previous_state

    def copy(self) -> "TDDFT":
        """Return an independent snapshot of the full strategy state."""
        return pycopy.deepcopy(self)

    def compute_result(
        self,
        geom: MolecularGeometry,
        request: ES_Request,
        result,
        previous: Optional[ES_Strategy] = None,
    ) -> None:
        """Store the active request, then run the shared sequencing contract."""
        self._request = request
        super().compute_result(geom, request, result, previous=previous)

    # ------------------------------------------------------------------ TDA/TDDFT

    def _ensure_td(self, previous: Optional[ES_Strategy] = None) -> None:
        if self._mf is None:
            raise ValueError("SCF must be run before computing TDDFT energies.")

        # When this instance was built empty (factory path), dens-project from the
        # previous strategy if set_geom had no local previous state.
        if previous is not None and isinstance(previous, TDDFT):
            prev_state = previous.get_state()
            if (
                self.get_previous_state() is None
                and prev_state is not None
                and prev_state.mol is not None
                and prev_state.dm is not None
                and self._mol is not None
            ):
                try:
                    dm0 = scf.addons.project_dm_nr2nr(prev_state.mol, prev_state.dm, self._mol)
                    self._mf.kernel(dm0=dm0)
                    if self._state is not None:
                        self._state.mf = self._mf
                        self._state.mo_coeff = np.asarray(self._mf.mo_coeff)
                        self._state.dm = np.asarray(self._mf.make_rdm1())
                except Exception:
                    pass

        if self._td is None:
            if self._use_tda:
                self._td = tdscf.TDA(self._mf)
            else:
                self._td = tdscf.TDDFT(self._mf)
            self._td.nstates = self._nexc
            self._td.kernel()
            if not all(np.atleast_1d(self._td.converged)):
                label = "TDA" if self._use_tda else "TDDFT"
                print(f"  WARNING: {label} not converged for all roots")

            state = self.get_state()
            if state is not None:
                state.td = self._td
                state.amplitudes = self._amplitudes(self._td)
                energies = np.concatenate(
                    ([float(self._mf.e_tot)], float(self._mf.e_tot) + np.asarray(self._td.e, dtype=float))
                )
                state.energies = energies
                state.mo_coeff = np.asarray(self._mf.mo_coeff)
                state.dm = np.asarray(self._mf.make_rdm1())

    @staticmethod
    def _amplitudes(td: Any) -> List[np.ndarray]:
        """Normalized TDA/TDDFT X amplitudes as (nocc, nvir) arrays."""
        X: List[np.ndarray] = []
        for n in range(len(td.e)):
            x = np.asarray(td.xy[n][0])
            norm = np.linalg.norm(x)
            X.append(x / norm if norm > 0 else x)
        return X

    def compute_H_el(self, previous: Optional[ES_Strategy] = None) -> np.ndarray:
        self._ensure_td(previous)
        n_total = self._request.n_singlets if self._request is not None else self._nstates
        energies = self._state.energies if self._state is not None else None
        if energies is None:
            raise RuntimeError("TDDFT energies are unavailable after _ensure_td.")
        if len(energies) < n_total:
            raise IndexError(
                f"Requested {n_total} TDDFT states, but only {len(energies)} are available."
            )
        return np.asarray(energies[:n_total], dtype=np.float64)

    def compute_energy(self, root: int) -> float:
        return float(self.compute_H_el()[root])

    def compute_gradient(self, root: int = 0) -> np.ndarray:
        self._ensure_td()
        if root < 0 or root >= self._nstates:
            raise IndexError(
                f"Requested TDDFT root {root}, but only {self._nstates} states are available."
            )
        if root == 0:
            return np.asarray(self._mf.nuc_grad_method().kernel(), dtype=np.float64)

        tdg = self._td.nuc_grad_method()
        # PySCF excited-state gradient indexing is 1-based over the LR roots.
        return np.asarray(tdg.kernel(state=root), dtype=np.float64)

    # ------------------------------------------------------------------ overlaps

    def compute_ao_overlap(self, right: "TDDFT") -> np.ndarray:
        if not isinstance(right, TDDFT):
            raise TypeError(f"right must be TDDFT, got {type(right).__name__}")
        if self._mol is None or right._mol is None:
            raise ValueError("Both TDDFT objects need molecule states before AO overlap.")
        return _compute_ao_overlap(self._mol, right._mol)

    def _time_overlap_nstates(self, right: "TDDFT") -> int:
        if right._request is not None:
            return int(right._request.n_singlets)
        if self._request is not None:
            return int(self._request.n_singlets)
        return min(int(self._nstates), int(right._nstates))

    def _factorized_cis_overlap(
        self,
        left: TDDFT_States,
        right: TDDFT_States,
        ao_overlap: np.ndarray,
        nstates: int,
    ) -> np.ndarray:
        if left.mo_coeff is None or right.mo_coeff is None:
            raise ValueError("MO coefficients are required for TDDFT time-overlap.")
        if left.amplitudes is None or right.amplitudes is None:
            raise ValueError("TDDFT amplitudes are required for time-overlap.")

        nexc = nstates - 1
        if len(left.amplitudes) < nexc or len(right.amplitudes) < nexc:
            raise ValueError(
                f"Requested {nexc} excited amplitudes, but only "
                f"{len(left.amplitudes)} left and {len(right.amplitudes)} right are available."
            )

        s_mo = left.mo_coeff.T @ ao_overlap @ right.mo_coeff
        nocc = left.amplitudes[0].shape[0]
        nmo = s_mo.shape[0]
        occ = np.arange(nocc)
        vir = np.arange(nocc, nmo)

        s_oo = s_mo[np.ix_(occ, occ)]
        s_vv = s_mo[np.ix_(vir, vir)]
        s_ov = s_mo[np.ix_(occ, vir)]
        s_vo = s_mo[np.ix_(vir, occ)]

        st = np.zeros((nstates, nstates), dtype=float)
        st[0, 0] = float(np.linalg.det(s_oo))
        for j in range(1, nstates):
            st[0, j] = float(np.sum(right.amplitudes[j - 1] * s_ov))
            st[j, 0] = float(np.sum(left.amplitudes[j - 1] * s_vo.T))
        for i in range(1, nstates):
            for j in range(1, nstates):
                st[i, j] = float(
                    np.einsum(
                        "ia,jb,ij,ab->",
                        left.amplitudes[i - 1],
                        right.amplitudes[j - 1],
                        s_oo,
                        s_vv,
                    )
                )
        return st

    @staticmethod
    def _align_time_overlap_phases(overlap: np.ndarray, phase_tol: float) -> np.ndarray:
        """Flip columns only when the diagonal sign is well determined."""
        overlap = np.array(overlap, copy=True)
        diag = np.diag(overlap)
        sgn = np.where(np.abs(diag) > phase_tol, np.sign(diag), 1.0)
        sgn[sgn == 0.0] = 1.0
        return overlap * sgn[None, :]

    def compute_time_overlap(self, right: "TDDFT") -> np.ndarray:
        if not isinstance(right, TDDFT):
            raise TypeError(f"right must be TDDFT, got {type(right).__name__}")

        left_state = self.get_state()
        right_state = right.get_state()
        if left_state is None or right_state is None:
            raise ValueError("Both TDDFT objects need completed states for time-overlap.")
        if left_state.td is None:
            self._ensure_td()
            left_state = self.get_state()
        if right_state.td is None:
            right._ensure_td()
            right_state = right.get_state()
        if left_state is None or right_state is None:
            raise ValueError("TDDFT states missing after ensuring LR calculation.")

        nstates = self._time_overlap_nstates(right)
        if nstates <= 0:
            raise ValueError(f"nstates must be positive, got {nstates}")

        ao_overlap = self.compute_ao_overlap(right)
        overlap = self._factorized_cis_overlap(left_state, right_state, ao_overlap, nstates)
        return self._align_time_overlap_phases(overlap, self._phase_tol)

    def compute_nac_vectors(self) -> np.ndarray:
        raise NotImplementedError(
            "TDDFT: analytic NAC vectors are not provided; "
            "request time_overlap instead (nac_update_method from St)."
        )
