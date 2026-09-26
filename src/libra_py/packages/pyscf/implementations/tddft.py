# *********************************************************************************
# * Copyright (C) 2026 Aniket Mandal
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
       Aniket Mandal

"""

from __future__ import annotations

import copy as pycopy
from dataclasses import dataclass
from typing import Any, List, Optional, Sequence

import numpy as np
from pyscf import dft, gto, scf, tdscf

from libra_py.packages.pyscf.interfaces import (
    ES_Request,
    ES_Result,
    ES_Strategy,
    MolecularGeometry,
)

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
    amplitudes: Optional[List[Any]] = None
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
    * State overlaps are determinant overlaps of normalized X-only CIS-like
      pseudo-wavefunctions.  Restricted response roots are expanded as
      spin-adapted singlets; unrestricted roots retain their separate alpha and
      beta response amplitudes.
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
    def nroots(self) -> int:
        """Number of spin-free roots, including the reference state."""
        return self._nstates

    @property
    def spin_multiplicity(self) -> int:
        """Spin multiplicity ``2*S+1`` represented by this strategy."""
        return self._spin + 1

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

        dft_cls = dft.RKS if self._spin == 0 else dft.UKS
        self._mf = dft_cls(self._mol, xc=self._xc)
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

        if self._spin == 0:
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

    def snapshot_state(self) -> None:
        """Save the current electronic structure for the next time overlap."""
        self._previous_state = pycopy.deepcopy(self._state)

    def compute_result(
        self,
        geom: MolecularGeometry,
        request: ES_Request,
    ) -> ES_Result:
        """Store the active request, then use the shared sequencing contract."""
        self._request = request
        return super().compute_result(geom, request)

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
            if self._nexc == 0:
                state = self.get_state()
                if state is not None:
                    state.amplitudes = []
                    state.energies = np.asarray([float(self._mf.e_tot)])
                return
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
    def _amplitudes(td: Any) -> List[Any]:
        """Normalized TDA/TDDFT X amplitudes as (nocc, nvir) arrays."""
        X: List[Any] = []
        for n in range(len(td.e)):
            x = td.xy[n][0]
            if isinstance(x, (tuple, list)):
                blocks = tuple(np.asarray(block) for block in x)
                norm = np.sqrt(sum(np.linalg.norm(block) ** 2 for block in blocks))
                X.append(tuple(block / norm if norm > 0 else block for block in blocks))
            else:
                x = np.asarray(x)
                norm = np.linalg.norm(x)
                X.append(x / norm if norm > 0 else x)
        return X

    def compute_H_el(self) -> np.ndarray:
        self._ensure_td()
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

    def compute_gradient(self, roots: Sequence[int]) -> list[np.ndarray]:
        """Nuclear gradients for ``roots``, returned in request order."""
        self._ensure_td()
        roots = list(roots)
        for root in roots:
            if root < 0 or root >= self._nstates:
                raise IndexError(
                    f"Requested TDDFT root {root}, but only "
                    f"{self._nstates} states are available."
                )

        ground_gradient = None
        excited_gradient = None
        gradients = []
        for root in roots:
            if root == 0:
                if ground_gradient is None:
                    ground_gradient = self._mf.nuc_grad_method()
                gradient = ground_gradient.kernel()
            else:
                if excited_gradient is None:
                    excited_gradient = self._td.nuc_grad_method()
                # PySCF excited-state gradient indexing is 1-based over LR roots.
                gradient = excited_gradient.kernel(state=root)
            gradients.append(np.asarray(gradient, dtype=np.float64))
        return gradients

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

    def _cis_overlap(
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

        if np.asarray(left.mo_coeff).ndim == 3:
            return self._unrestricted_cis_overlap(
                left, right, ao_overlap, nstates
            )
        return self._restricted_singlet_cis_overlap(
            left, right, ao_overlap, nstates
        )

    @staticmethod
    def _restricted_singlet_configurations(amplitude, nmo, nocc):
        """Expand one restricted X-only root in normalized singlet determinants.

        The reference is ``|Phi_0>``.  A normalized restricted response root is

        ``sum_ia X_ia (|Phi_i^a(alpha)> + |Phi_i^a(beta)>)/sqrt(2)``,

        where ``sum_ia |X_ia|^2 = 1`` after :meth:`_amplitudes`.  Occupations
        retain the replaced orbital's list position so the determinant phase is
        consistent with the excitation-operator convention.
        """
        occ = tuple(range(nocc))
        if amplitude is None:
            return [(1.0, occ, occ)]

        amplitude = np.asarray(amplitude)
        expected = (nocc, nmo - nocc)
        if amplitude.shape != expected:
            raise ValueError(
                "Restricted TDDFT amplitude shape does not match the MO "
                f"occupation: expected {expected}, got {amplitude.shape}."
            )

        configs = []
        spin_factor = 1.0 / np.sqrt(2.0)
        for i in range(amplitude.shape[0]):
            for a in range(amplitude.shape[1]):
                coefficient = amplitude[i, a]
                if abs(coefficient) < 1e-14:
                    continue
                excited = list(occ)
                excited[i] = nocc + a
                excited = tuple(excited)
                coefficient = coefficient * spin_factor
                configs.append((coefficient, excited, occ))
                configs.append((coefficient, occ, excited))
        return configs

    @staticmethod
    def _configuration_overlap(configs_i, configs_j, s_mo_alpha, s_mo_beta):
        """Contract two alpha/beta determinant expansions."""
        dtype = np.result_type(s_mo_alpha, s_mo_beta, complex)
        value = np.asarray(0.0, dtype=dtype)[()]
        for coeff_i, occ_ia, occ_ib in configs_i:
            for coeff_j, occ_ja, occ_jb in configs_j:
                value += (
                    np.conjugate(coeff_i) * coeff_j
                    * np.linalg.det(s_mo_alpha[np.ix_(occ_ia, occ_ja)])
                    * np.linalg.det(s_mo_beta[np.ix_(occ_ib, occ_jb)])
                )
        return np.real_if_close(value)

    def _restricted_singlet_cis_overlap(
        self, left, right, ao_overlap, nstates
    ):
        """Determinant overlap of normalized restricted singlet pseudo-states."""
        left_mo = np.asarray(left.mo_coeff)
        right_mo = np.asarray(right.mo_coeff)
        s_mo = left_mo.T.conj() @ ao_overlap @ right_mo
        nocc_left = int(np.count_nonzero(np.asarray(left.mf.mo_occ)))
        nocc_right = int(np.count_nonzero(np.asarray(right.mf.mo_occ)))
        if nocc_left != nocc_right:
            raise ValueError(
                "Restricted TDDFT snapshots have different occupied-orbital counts: "
                f"{nocc_left} and {nocc_right}."
            )

        left_states = [None] + list(left.amplitudes[: nstates - 1])
        right_states = [None] + list(right.amplitudes[: nstates - 1])
        expected_left = (nocc_left, left_mo.shape[1] - nocc_left)
        expected_right = (nocc_right, right_mo.shape[1] - nocc_right)
        for side, states, expected in (
            ("left", left_states, expected_left),
            ("right", right_states, expected_right),
        ):
            for amplitude in states[1:]:
                if np.asarray(amplitude).shape != expected:
                    raise ValueError(
                        f"Restricted TDDFT {side} amplitude shape does not match "
                        f"the MO occupation: expected {expected}, got "
                        f"{np.asarray(amplitude).shape}."
                    )
        occupied = s_mo[:nocc_left, :nocc_right]
        # The determinant-minor identities below are exact when the occupied
        # block is invertible and avoid the quadratic number of nocc-by-nocc
        # determinants in an explicit singles expansion.  Near singularities,
        # use the general determinant contraction instead.
        try:
            condition = np.linalg.cond(occupied)
        except np.linalg.LinAlgError:
            condition = np.inf
        if np.isfinite(condition) and condition < 1.0 / np.sqrt(np.finfo(float).eps):
            return self._restricted_singlet_cis_overlap_minors(
                left_states, right_states, s_mo, nocc_left
            )
        return self._restricted_singlet_cis_overlap_determinants(
            left_states, right_states, left_mo.shape[1], right_mo.shape[1],
            nocc_left, s_mo,
        )

    def _restricted_singlet_cis_overlap_minors(
        self, left_states, right_states, s_mo, nocc
    ):
        """Evaluate restricted singlet overlaps from exact determinant minors."""
        occupied = s_mo[:nocc, :nocc]
        occ_virtual = s_mo[:nocc, nocc:]
        virtual_occ = s_mo[nocc:, :nocc]
        virtual_virtual = s_mo[nocc:, nocc:]
        inverse = np.linalg.inv(occupied)
        determinant = np.linalg.det(occupied)
        right_minor = inverse @ occ_virtual
        left_minor = virtual_occ @ inverse
        double_remainder = virtual_virtual - virtual_occ @ inverse @ occ_virtual
        reference_overlap = determinant * determinant

        dtype = np.result_type(s_mo, *(left_states[1:] + right_states[1:]))
        overlap = np.zeros((len(left_states), len(right_states)), dtype=dtype)
        overlap[0, 0] = reference_overlap
        left_contractions = [None]
        right_contractions = [None]
        for i, amplitude in enumerate(left_states[1:], start=1):
            contraction = np.einsum(
                "ia,ai->", np.conjugate(amplitude), left_minor, optimize=True
            )
            left_contractions.append(contraction)
            overlap[i, 0] = np.sqrt(2.0) * reference_overlap * contraction
        for j, amplitude in enumerate(right_states[1:], start=1):
            contraction = np.einsum(
                "jb,jb->", amplitude, right_minor, optimize=True
            )
            right_contractions.append(contraction)
            overlap[0, j] = np.sqrt(2.0) * reference_overlap * contraction
        for i, amplitude_i in enumerate(left_states[1:], start=1):
            for j, amplitude_j in enumerate(right_states[1:], start=1):
                connected = np.einsum(
                    "ia,jb,ji,ab->",
                    np.conjugate(amplitude_i), amplitude_j, inverse,
                    double_remainder, optimize=True,
                )
                overlap[i, j] = reference_overlap * (
                    2.0 * left_contractions[i] * right_contractions[j] + connected
                )
        return np.real_if_close(overlap)

    def _restricted_singlet_cis_overlap_determinants(
        self, left_states, right_states, nmo_left, nmo_right, nocc, s_mo
    ):
        """General determinant-sum fallback for a singular occupied block."""
        dtype = np.result_type(s_mo, *(left_states[1:] + right_states[1:]))
        overlap = np.zeros((len(left_states), len(right_states)), dtype=dtype)
        left_configs = [
            self._restricted_singlet_configurations(
                amplitude, nmo_left, nocc
            )
            for amplitude in left_states
        ]
        right_configs = [
            self._restricted_singlet_configurations(
                amplitude, nmo_right, nocc
            )
            for amplitude in right_states
        ]
        for i, configs_i in enumerate(left_configs):
            for j, configs_j in enumerate(right_configs):
                overlap[i, j] = self._configuration_overlap(
                    configs_i, configs_j, s_mo, s_mo
                )
        return np.real_if_close(overlap)

    @staticmethod
    def _unrestricted_configurations(amplitude, nmo, nocc):
        """Return determinant coefficients and occupations for one LR state."""
        occ_a = tuple(range(nocc[0]))
        occ_b = tuple(range(nocc[1]))
        if amplitude is None:
            return [(1.0, occ_a, occ_b)]
        configs = []
        for spin, block in enumerate(amplitude):
            for i in range(block.shape[0]):
                for a in range(block.shape[1]):
                    if abs(block[i, a]) < 1e-14:
                        continue
                    occupations = [list(occ_a), list(occ_b)]
                    occupations[spin][i] = nocc[spin] + a
                    configs.append(
                        (block[i, a], tuple(occupations[0]), tuple(occupations[1]))
                    )
        return configs

    def _unrestricted_cis_overlap(self, left, right, ao_overlap, nstates):
        """Determinant overlap for UKS X-only response amplitudes."""
        left_mo = np.asarray(left.mo_coeff)
        right_mo = np.asarray(right.mo_coeff)
        s_mo = tuple(
            left_mo[spin].T.conj() @ ao_overlap @ right_mo[spin]
            for spin in range(2)
        )
        nmo_left = tuple(left_mo[spin].shape[1] for spin in range(2))
        nmo_right = tuple(right_mo[spin].shape[1] for spin in range(2))
        nocc_left = tuple(
            int(np.count_nonzero(left.mf.mo_occ[spin])) for spin in range(2)
        )
        nocc_right = tuple(
            int(np.count_nonzero(right.mf.mo_occ[spin])) for spin in range(2)
        )
        if nocc_left != nocc_right:
            raise ValueError(
                "Unrestricted TDDFT snapshots have different alpha/beta "
                f"occupied-orbital counts: {nocc_left} and {nocc_right}."
            )
        left_states = [None] + list(left.amplitudes[: nstates - 1])
        right_states = [None] + list(right.amplitudes[: nstates - 1])
        amplitude_dtypes = [
            np.asarray(block).dtype
            for amplitude in left_states[1:] + right_states[1:]
            for block in amplitude
        ]
        dtype = np.result_type(*(list(s_mo) + amplitude_dtypes))
        overlap = np.zeros((nstates, nstates), dtype=dtype)
        for i, amplitude_i in enumerate(left_states):
            configs_i = self._unrestricted_configurations(
                amplitude_i, nmo_left, nocc_left
            )
            for j, amplitude_j in enumerate(right_states):
                configs_j = self._unrestricted_configurations(
                    amplitude_j, nmo_right, nocc_right
                )
                overlap[i, j] = self._configuration_overlap(
                    configs_i, configs_j, s_mo[0], s_mo[1]
                )
        return np.real_if_close(overlap)

    @staticmethod
    def _align_time_overlap_phases(overlap: np.ndarray, phase_tol: float) -> np.ndarray:
        """Flip columns only when the diagonal sign is well determined."""
        overlap = np.array(overlap, copy=True)
        diag = np.diag(overlap)
        sgn = np.where(np.abs(diag) > phase_tol, np.sign(diag), 1.0)
        sgn[sgn == 0.0] = 1.0
        return overlap * sgn[None, :]

    def compute_time_overlap(self, state1: object, state2: object) -> np.ndarray:
        """Return ``<Psi(t)|Psi(t+dt)>`` from previous/current snapshots."""
        if not isinstance(state1, TDDFT_States) or not isinstance(state2, TDDFT_States):
            raise TypeError(
                "state1/state2 must be TDDFT_States, got "
                f"{type(state1).__name__} / {type(state2).__name__}."
            )
        current_state = state1
        previous_state = state2
        nstates = int(self._request.n_singlets) if self._request else self._nstates
        if nstates <= 0:
            raise ValueError(f"nstates must be positive, got {nstates}")

        ao_overlap = _compute_ao_overlap(previous_state.mol, current_state.mol)
        overlap = self._cis_overlap(
            previous_state, current_state, ao_overlap, nstates
        )
        return self._align_time_overlap_phases(overlap, self._phase_tol)

    def compute_nac_vectors(self) -> np.ndarray:
        raise NotImplementedError(
            "TDDFT: analytic NAC vectors are not provided; "
            "request time_overlap instead (nac_update_method from St)."
        )
