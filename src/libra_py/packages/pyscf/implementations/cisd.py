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
.. module:: pyscf.implementations.cisd
   :platform: Unix, Windows
   :synopsis: PySCF CISD implementation for the universal ES strategy.
.. moduleauthor::
       Jieyang Gu <jieyanggu792@gmail.com>

"""
from __future__ import annotations
import copy as pycopy
from dataclasses import dataclass
from functools import reduce
from typing import Any, Optional
import numpy as np
from pyscf import ci, gto, scf
from libra_py.packages.pyscf.interfaces import ES_Strategy, ES_Request, MolecularGeometry

# Backward-compatible alias used by older imports.
ElectronicStructureStrategy = ES_Strategy

BOHR_TO_ANG = 0.529177210903


@dataclass
class CISD_States:
    """Snapshot of the CISD calculation at one geometry."""

    mol: Optional[Any] = None
    mf: Optional[Any] = None
    ci: Optional[Any] = None


class CISD(ES_Strategy):
    """PySCF-based CISD backend for the universal ES interface."""

    def __init__(
        self,
        mol: Optional[Any] = None,
        nroots: int = 1,
        basis: str = "sto-3g",
        unit: str = "Angstrom",
        charge: int = 0,
    ) -> None:
        self._mol: Optional[Any] = mol
        self._nroots: int = nroots
        self._basis: str = basis
        self._unit: str = unit
        self._charge: int = int(charge)
        self._geom: Optional[MolecularGeometry] = None
        self._request: Optional[ES_Request] = None
        self._state: Optional[CISD_States] = None
        self._previous_state: Optional[CISD_States] = None

    # --------------------------------------------------------------- settings
    def _n_total(self) -> int:
        """Number of electronic states requested for the current calculation."""
        if self._request is not None:
            return int(self._request.n_singlets)
        return int(self._nroots)

    @staticmethod
    def _as_ci_vector_list(ci_data: Any) -> list[np.ndarray]:
        if ci_data is None:
            return []
        if isinstance(ci_data, (list, tuple)):
            return [np.asarray(vec).copy() for vec in ci_data]
        return [np.asarray(ci_data).copy()]

    def _ci_guess(self) -> Optional[Any]:
        """CI guess from the previous geometry's CISD vectors (hot start)."""
        prev_state = self.get_previous_state()
        if prev_state is None or prev_state.ci is None:
            return None
        ci_guess = self._as_ci_vector_list(getattr(prev_state.ci, "ci", None))
        if not ci_guess:
            return None
        if self._n_total() == 1:
            return ci_guess[0]
        if len(ci_guess) < self._n_total():
            return None
        return ci_guess[: self._n_total()]

    # ---------------------------------------------------------------- geometry
    def set_geom(self, geom: MolecularGeometry) -> None:
        self._geom = geom
        self._previous_state = pycopy.deepcopy(self._state)

        self._ci = None

        coords_bohr = getattr(geom, "coords_bohr", None)
        if coords_bohr is None:
            # Legacy: coords_angstrom holds values in the strategy's configured unit.
            coords = np.asarray(getattr(geom, "coords_angstrom"), dtype=float)
        else:
            coords = np.asarray(coords_bohr, dtype=float)
            if self._unit.lower().startswith("ang"):
                coords = coords / BOHR_TO_ANG  # bohr -> angstrom

        self._mol = gto.M(
            atom=";".join(
                f"{label} {coord[0]} {coord[1]} {coord[2]}"
                for label, coord in zip(geom.atom_labels, coords)
            ),
            basis=self._basis,
            unit=self._unit,
            charge=self._charge,
            spin=0,
        )

        self._mf = scf.RHF(self._mol).run(verbose=0)
        self._state = CISD_States(mol=self._mol, mf=self._mf, ci=None)

    def get_geom(self) -> MolecularGeometry:
        if self._geom is None:
            raise ValueError("Geometry has not been set.")
        return self._geom

    # ------------------------------------------------------------------ state
    def get_state(self) -> Optional[CISD_States]:
        return self._state

    def get_previous_state(self) -> Optional[CISD_States]:
        return self._previous_state

    def snapshot_state(self) -> None:
        """Save the current state of the calculation to a previous-state snapshot."""
        self._previous_state = pycopy.deepcopy(self._state)

    def copy(self) -> "CISD":
        """Return an independent snapshot of the full strategy state."""
        return pycopy.deepcopy(self)

    # --------------------------------------------------------------- required
    def _build_ci(self) -> None:
        if self._mf is None:
            raise ValueError("HF must be run before computing CISD energies.")
        if self._ci is None:
            self._ci = ci.cisd.CISD(self._mf)
        self._ci.nroots = self._n_total()
        self._ci.verbose = 0
        self._ci.kernel(ci0=self._ci_guess())
        state = self.get_state()
        if state is not None:
            state.ci = self._ci

    def compute_H_el(self) -> np.ndarray:
        self._build_ci()
        n_total = self._n_total()
        e_tot = getattr(self._ci, "e_tot", None)
        if isinstance(e_tot, (list, tuple, np.ndarray)):
            energies = np.asarray(e_tot, dtype=np.float64)
        else:
            energies = np.asarray([float(e_tot)], dtype=np.float64)
        if len(energies) < n_total:
            raise IndexError(
                f"Requested {n_total} states, but only {len(energies)} are available."
            )
        return energies[:n_total]

    def compute_energy(self, root: int) -> float:
        """Legacy convenience: energy of a single root."""
        return float(self.compute_H_el()[root])

    def compute_gradient(self, root: int = 0) -> np.ndarray:
        self._build_ci()
        ci_vectors = self._ci.ci
        if isinstance(ci_vectors, (list, tuple)):
            ci_roots = [np.asarray(vec) for vec in ci_vectors]
        else:
            ci_roots = [np.asarray(ci_vectors)]
        if root < 0 or root >= len(ci_roots):
            raise IndexError(
                f"Requested root {root}, but only {len(ci_roots)} roots are available."
            )
        if len(ci_roots) == 1:
            if root != 0:
                raise IndexError("Only root 0 is available for a single-state CISD calculation.")
            return np.asarray(self._ci.nuc_grad_method().kernel())
        return np.asarray(self._ci.nuc_grad_method().kernel(state=root))

    def compute_all_gradients(self) -> list[np.ndarray]:
        """Compute the nuclear gradients for all requested states together."""
        self._build_ci()
        n_total = self._n_total()
        ci_vectors = self._ci.ci
        n_avail = len(ci_vectors) if isinstance(ci_vectors, (list, tuple)) else 1
        if n_total > n_avail:
            raise IndexError(
                f"Requested {n_total} states, but only {n_avail} are available."
            )
        return [self.compute_gradient(root) for root in range(n_total)]

    # ------------------------------------------------------------- time overlap

    def _compute_ao_overlap(
        self,
        prev_state: CISD_States,
        curr_state: CISD_States,
    ) -> np.ndarray:
        """Compute the AO overlap matrix between the previous and current geometries."""
        prev_mol = prev_state.mol
        curr_mol = curr_state.mol
        if prev_mol is None or curr_mol is None:
            raise ValueError("Both previous and current molecule objects are required.")
        return gto.intor_cross("int1e_ovlp", prev_mol, curr_mol)

    def compute_time_overlap(self, state1: object, state2: object) -> np.ndarray:
        """Adiabatic time-overlap between two CISD state snapshots.

        state1: the current-state snapshot (``CISD_States``).
        state2: the previous-state snapshot (``CISD_States``).

        Returns
        -------
        np.ndarray
            Shape ``(n_total, n_total)``.
        """
        if not isinstance(state1, CISD_States) or not isinstance(state2, CISD_States):
            raise TypeError(
                f"state1/state2 must be CISD_States, got "
                f"{type(state1).__name__} / {type(state2).__name__}."
            )
        prev_state = state2  # previous geometry
        curr_state = state1  # current geometry

        if prev_state.mol is None or curr_state.mol is None:
            raise ValueError("Molecule states are required for time-overlap computation.")
        if prev_state.mf is None or curr_state.mf is None:
            raise ValueError("HF states are required for time-overlap computation.")
        if prev_state.ci is None or curr_state.ci is None:
            raise ValueError("Both CISD states must be run before computing time-overlap.")

        nroots = self._n_total()
        if nroots <= 0:
            raise ValueError(f"nroots must be positive, got {nroots}")


        s12_ao = self._compute_ao_overlap(prev_state, curr_state)
        s12_mo = reduce(np.dot, (prev_state.mf.mo_coeff.T, s12_ao, curr_state.mf.mo_coeff))

        prev_ci_list = self._as_ci_vector_list(prev_state.ci.ci) or []
        curr_ci_list = self._as_ci_vector_list(curr_state.ci.ci) or []
        if len(prev_ci_list) < nroots or len(curr_ci_list) < nroots:
            raise ValueError(
                f"Requested {nroots} roots, but only {len(prev_ci_list)} previous "
                f"and {len(curr_ci_list)} current roots are available."
            )

        nmo = curr_state.ci.nmo
        nelec = curr_state.mol.nelectron // 2

        overlap = np.zeros((nroots, nroots), dtype=float)
        for i in range(nroots):
            for j in range(nroots):
                overlap[i, j] = ci.cisd.overlap(
                    prev_ci_list[i],
                    curr_ci_list[j],
                    nmo,
                    nelec,
                    s12_mo,
                )

        overlap = np.asarray(np.real_if_close(overlap))
        for j in range(nroots):
            if overlap[j, j] < 0:
                overlap[:, j] = -overlap[:, j]
        return overlap

    def time_overlap_matrix(self, nroots: int) -> np.ndarray:
        """Legacy same-object API: compare the cached previous state to the current state."""
        if nroots <= 0:
            raise ValueError(f"nroots must be positive, got {nroots}")
        prev_state = self.get_previous_state()
        curr_state = self.get_state()
        if prev_state is None or prev_state.ci is None:
            raise ValueError(
                "Previous CISD state must exist for time-overlap computation. "
                "Run at least two geometries before requesting time overlaps."
            )
        if curr_state is None or curr_state.ci is None:
            raise ValueError("Current CISD state must exist for time-overlap computation.")
        if nroots != self._n_total():
            raise ValueError(
                f"Requested {nroots} roots, but the strategy is configured for "
                f"{self._n_total()} states."
            )
        return self.compute_time_overlap(curr_state, prev_state)


if __name__ == "__main__":
    # Smoke tests
    # Geometry coordinates are given explicitly in Bohr.
    geom1 = MolecularGeometry(
        atom_labels=('He', 'H'),
        coords_bohr=np.array([
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 1.46379],
        ], dtype=np.float64),
    )

    geom2 = MolecularGeometry(
        atom_labels=('He', 'H'),
        coords_bohr=np.array([
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 1.88973],
        ], dtype=np.float64),
    )

    # 1) Demonstrate: set geom1 → run SCF → print the AO (overlap) matrix.
    demo = CISD(
        nroots=3,
        basis="sto-3g",
        charge=1,
        unit="Bohr",
    )
    demo.set_geom(geom1)  # sets the geometry and triggers the HF SCF
    ao_overlap = demo._mol.intor("int1e_ovlp")
    np.set_printoptions(precision=4, suppress=True)
    print("\n[set geom1 + run SCF] natm =", demo._mol.natm,
          "  HF energy = %.8f Ha" % demo._mf.e_tot)
    print("[AO overlap matrix] shape =", ao_overlap.shape)
    print(ao_overlap)
    assert demo._mol.natm == 2
    assert demo._mf.converged

    # 2) Full two-geometry sequencing flow on a fresh strategy.
    # (geom1 and geom2 already defined above — reuse them.)
    cisd = CISD(
        nroots=3,
        basis="sto-3g",
        charge=1,
        unit="Bohr",
    )
    request = ES_Request(
        n_singlets=3,
        gradient_state="all",
        nacv=False,
        time_overlap=True,
    )

    result1 = cisd.compute_result(geom1, request)
    print("\nResult 1 H_el:", result1.H_el)
    print("Result 1 gradients:\n", result1.gradients)
    print("Result 1 time_overlap:", result1.time_overlap)

    result2 = cisd.compute_result(geom2, request)
    print("\nResult 2 H_el:", result2.H_el)
    print("Result 2 gradients:\n", result2.gradients)
    print("Result 2 time_overlap:\n", result2.time_overlap)

    assert result1.time_overlap is None
    assert result2.time_overlap is not None and result2.time_overlap.shape == (3, 3)
    assert isinstance(cisd, ES_Strategy)
    print("\nCISD __main__ smoke test passed.")

    

