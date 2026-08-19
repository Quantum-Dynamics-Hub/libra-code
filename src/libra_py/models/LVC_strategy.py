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
.. module:: models.LVC_strategy
   :platform: Unix, Windows
   :synopsis: Analytic Linear Vibronic Coupling model as an ES_Strategy.

   A dependency-free stand-in for a real electronic-structure backend. It
   provides every quantity the ES interface can be asked for -- energies,
   gradients, NAC vectors, time-overlaps -- from closed-form expressions, so
   the NAMD machinery can be exercised without PySCF, DFTB+, or any external
   code.

   The whole model is three statements:

       1. Each diabatic state is a parabola:   V_i(q) = E_i + sum_n 1/2 k_in (q_n - c_in)^2
       2. The states are coupled:              V_ij(q) = constant, or linear in q
       3. Adiabatic states are the eigenvalues/eigenvectors of that matrix.

   Everything else in this file is bookkeeping around those three lines.

   Usage::

       from libra_py.models.LVC_strategy import LVC, Parabola, Coupling

   It lives under ``libra_py.models`` rather than beside the PySCF backends in
   ``packages/pyscf/implementations/`` for one concrete reason: that package's
   ``__init__`` imports the PySCF-backed classes, so importing anything from it
   requires PySCF. ``libra_py.models.__init__`` imports nothing and
   ``libra_py.packages.pyscf.interfaces`` needs only NumPy, so this module runs
   with no quantum-chemistry code installed -- which is the whole point of it.

.. moduleauthor::
       Jieyang Gu <jieyanggu792@gmail.com>
"""

from __future__ import annotations

import copy as pycopy
import math
from dataclasses import dataclass, field
from typing import Literal, Optional, Sequence

import numpy as np

from libra_py.packages.pyscf.interfaces import (
    ES_Request,
    ES_Strategy,
    MolecularGeometry,
)

__all__ = ["Parabola", "Coupling", "LVC_State", "LVC"]


# =============================================================================
# SECTION 1 -- What you specify when you instantiate the model
# =============================================================================

def _broadcast(value, ndof: int, name: str) -> np.ndarray:
    """Accept a scalar (same for every dof) or one value per dof."""
    array = np.atleast_1d(np.asarray(value, dtype=np.float64))
    if array.size == 1:
        return np.full(ndof, float(array[0]))
    if array.size != ndof:
        raise ValueError(
            f"{name}: expected a scalar or {ndof} values, got {array.size}."
        )
    return array.astype(np.float64, copy=True)


@dataclass(frozen=True)
class Parabola:
    """One diabatic state: an ndof-dimensional harmonic well.

        V(q) = energy_min + sum_n 1/2 * mass_n * omega_n^2 * (q_n - center_n)^2

    Attributes
    ----------
    energy_min : float
        Energy at the bottom of the well [Hartree].
    center : float or sequence of float
        Position of the minimum along each dof [Bohr]. A scalar applies to
        every dof.
    omega : float or sequence of float
        Harmonic frequency along each dof [Hartree]. The curvature (force
        constant) is ``k_n = mass_n * omega_n**2``; ``mass`` is supplied once
        on the model, since it belongs to the nuclei rather than to a state.
    """

    energy_min: float = 0.0
    center: float | Sequence[float] = 0.0
    omega: float | Sequence[float] = 0.01

    def centers(self, ndof: int) -> np.ndarray:
        return _broadcast(self.center, ndof, "Parabola.center")

    def omegas(self, ndof: int) -> np.ndarray:
        return _broadcast(self.omega, ndof, "Parabola.omega")


@dataclass(frozen=True)
class Coupling:
    """Off-diagonal element between two diabatic states.

    Use the constructors rather than the raw fields:

        Coupling.constant(0.005)                  ->  V_ij = 0.005
        Coupling.linear(slopes=0.002)             ->  V_ij = 0.002 * q
        Coupling.linear(slopes=[1e-3, 2e-3])      ->  V_ij = sum_n slope_n * q_n
        Coupling.linear(slopes=1e-3, origin=1.5)  ->  V_ij = 1e-3 * (q - 1.5)
        Coupling.none()                           ->  V_ij = 0  (states never mix)

    ``constant`` gives the textbook avoided crossing; ``linear`` reproduces the
    genuine LVC form, where the coupling vanishes on a seam and the two
    adiabatic surfaces can touch (a conical intersection for ndof >= 2).
    """

    kind: Literal["constant", "linear"] = "constant"
    strength: float = 0.0
    slopes: float | Sequence[float] = 0.0
    origin: float | Sequence[float] = 0.0

    @classmethod
    def none(cls) -> "Coupling":
        return cls(kind="constant", strength=0.0)

    @classmethod
    def constant(cls, strength: float) -> "Coupling":
        return cls(kind="constant", strength=float(strength))

    @classmethod
    def linear(cls, slopes, origin=0.0) -> "Coupling":
        return cls(kind="linear", slopes=slopes, origin=origin)

    # -- value and derivative, the only two things the model asks of a coupling
    def value(self, q: np.ndarray, ndof: int) -> float:
        if self.kind == "constant":
            return float(self.strength)
        slopes = _broadcast(self.slopes, ndof, "Coupling.slopes")
        origin = _broadcast(self.origin, ndof, "Coupling.origin")
        return float(np.dot(slopes, q - origin))

    def gradient(self, ndof: int) -> np.ndarray:
        if self.kind == "constant":
            return np.zeros(ndof)
        return _broadcast(self.slopes, ndof, "Coupling.slopes")


@dataclass
class LVC_State:
    """Snapshot of the model at one geometry (the ES_Strategy 'state')."""

    coords: np.ndarray = field(default_factory=lambda: np.zeros(0))
    H_dia: np.ndarray = field(default_factory=lambda: np.zeros((0, 0)))
    dH_dia: np.ndarray = field(default_factory=lambda: np.zeros((0, 0, 0)))
    energies: np.ndarray = field(default_factory=lambda: np.zeros(0))
    U: np.ndarray = field(default_factory=lambda: np.zeros((0, 0)))


# =============================================================================
# SECTION 2 -- The model
# =============================================================================

class LVC(ES_Strategy):
    """Analytic vibronic-coupling model exposed through the ES interface.

    Parameters
    ----------
    ndof : int
        Number of nuclear degrees of freedom. This is the "1D / 2D / 3D" knob.
    states : sequence of Parabola
        One entry per diabatic state. Two entries is the classic LVC problem.
    coupling : Coupling or dict[(i, j) -> Coupling], optional
        A single Coupling is applied to every off-diagonal pair; a dict sets
        pairs individually (unlisted pairs are zero). Defaults to no coupling.
    mass : float or sequence of float
        Nuclear mass per dof [a.u. of mass]. Enters the curvature as
        ``k_n = mass_n * omega_n**2`` and must match the masses handed to the
        Libra dynamics driver.

    Examples
    --------
    A 1-D two-state avoided crossing -- two displaced wells, constant coupling::

        model = LVC(
            ndof=1,
            states=[
                Parabola(energy_min=0.000, center=-1.0, omega=0.005),
                Parabola(energy_min=0.010, center=+1.0, omega=0.005),
            ],
            coupling=Coupling.constant(0.002),
            mass=2000.0,
        )

    A 2-D conical intersection -- linear coupling along the second mode::

        model = LVC(
            ndof=2,
            states=[
                Parabola(energy_min=0.0, center=(-1.0, 0.0), omega=(0.005, 0.005)),
                Parabola(energy_min=0.0, center=(+1.0, 0.0), omega=(0.005, 0.005)),
            ],
            coupling=Coupling.linear(slopes=(0.0, 0.003)),
            mass=(2000.0, 2000.0),
        )

    Notes
    -----
    The diabatic basis is geometry-independent, so the adiabatic time-overlap
    ``<phi_i(t)|phi_j(t+dt)> = (U_prev.T @ U_curr)[i, j]`` is *exact* rather
    than approximate -- which makes this model a clean reference for testing
    state tracking, phase correction, and NAC-from-time-overlap code paths.
    """

    # ---------------------------------------------------------------- setup
    def __init__(
        self,
        ndof: int,
        states: Sequence[Parabola],
        coupling=None,
        mass: float | Sequence[float] = 1.0,
    ) -> None:
        if int(ndof) < 1:
            raise ValueError(f"ndof must be >= 1, got {ndof}.")
        if len(states) < 1:
            raise ValueError("At least one diabatic state (Parabola) is required.")

        self.ndof = int(ndof)
        self.nstates = len(states)

        self._energy_min = np.array(
            [float(state.energy_min) for state in states], dtype=np.float64
        )
        self._center = np.array(
            [state.centers(self.ndof) for state in states], dtype=np.float64
        )
        self._omega = np.array(
            [state.omegas(self.ndof) for state in states], dtype=np.float64
        )
        self._mass = _broadcast(mass, self.ndof, "mass")

        # k[i, n] = mass[n] * omega[i, n]**2
        self._k = self._mass[None, :] * self._omega**2

        self._coupling = self._normalize_coupling(coupling)

        # ES_Strategy bookkeeping
        self._geom: Optional[MolecularGeometry] = None
        self._coords: Optional[np.ndarray] = None
        self._request: Optional[ES_Request] = None
        self._state: Optional[LVC_State] = None
        self._previous_state: Optional[LVC_State] = None

    def _normalize_coupling(self, coupling) -> dict:
        """Turn any accepted spelling into {(i, j): Coupling} with i < j."""
        pairs = [
            (i, j)
            for i in range(self.nstates)
            for j in range(i + 1, self.nstates)
        ]
        if coupling is None:
            return {}
        if isinstance(coupling, Coupling):
            return {pair: coupling for pair in pairs}
        if isinstance(coupling, dict):
            normalized = {}
            for (i, j), value in coupling.items():
                i, j = int(i), int(j)
                if not (0 <= i < self.nstates and 0 <= j < self.nstates):
                    raise ValueError(
                        f"coupling pair {(i, j)} is outside [0, {self.nstates})."
                    )
                if i == j:
                    raise ValueError("coupling pairs must connect distinct states.")
                normalized[(min(i, j), max(i, j))] = value
            return normalized
        raise TypeError(
            "coupling must be a Coupling, a dict of {(i, j): Coupling}, or None."
        )

    @property
    def natoms(self) -> int:
        """Atoms the dofs are packed into, for the (natoms, 3) ES conventions."""
        return int(math.ceil(self.ndof / 3))

    # ============================================================== the physics
    def _diabatic(self, q: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Diabatic matrix and its nuclear derivatives.

        Returns
        -------
        H : (nstates, nstates)              Hartree
        dH : (ndof, nstates, nstates)       Hartree / Bohr
        """
        H = np.zeros((self.nstates, self.nstates), dtype=np.float64)
        dH = np.zeros((self.ndof, self.nstates, self.nstates), dtype=np.float64)

        # 1. Diagonal: one parabola per diabatic state.
        displacement = q[None, :] - self._center          # (nstates, ndof)
        H[np.diag_indices(self.nstates)] = (
            self._energy_min + 0.5 * np.sum(self._k * displacement**2, axis=1)
        )
        for i in range(self.nstates):
            dH[:, i, i] = self._k[i] * displacement[i]

        # 2. Off-diagonal: the coupling that lets the states mix.
        for (i, j), coupling in self._coupling.items():
            H[i, j] = H[j, i] = coupling.value(q, self.ndof)
            slope = coupling.gradient(self.ndof)
            dH[:, i, j] = slope
            dH[:, j, i] = slope

        return H, dH

    @staticmethod
    def _fix_phase(U: np.ndarray) -> np.ndarray:
        """Make eigenvector signs reproducible.

        ``numpy.linalg.eigh`` fixes each column only up to a sign. A geometry-
        local convention -- largest-magnitude component positive -- is used so
        repeated runs agree, while genuine rotations through a crossing still
        flip the sign. The signs are deliberately *not* aligned against the
        previous step: that is the job of Libra's phase-correction code, which
        this model exists to test.
        """
        U = np.array(U, dtype=np.float64, copy=True)
        for column in range(U.shape[1]):
            pivot = int(np.argmax(np.abs(U[:, column])))
            if U[pivot, column] < 0.0:
                U[:, column] *= -1.0
        return U

    def _solve(self) -> LVC_State:
        """Build (and cache) the adiabatic solution at the current geometry."""
        if self._state is not None:
            return self._state
        if self._coords is None:
            raise ValueError("Geometry has not been set; call set_geom first.")

        H, dH = self._diabatic(self._coords)
        energies, U = np.linalg.eigh(H)          # ascending, orthonormal columns
        self._state = LVC_State(
            coords=self._coords.copy(),
            H_dia=H,
            dH_dia=dH,
            energies=energies,
            U=self._fix_phase(U),
        )
        return self._state

    def _pack(self, per_dof: np.ndarray) -> np.ndarray:
        """Lay ndof values out as (natoms, 3), zero-padding the tail."""
        packed = np.zeros(3 * self.natoms, dtype=np.float64)
        packed[: self.ndof] = per_dof
        return packed.reshape(self.natoms, 3)


# =============================================================================
# SECTION 3 -- ES_Strategy: geometry and state bookkeeping
# =============================================================================

    def set_geom(self, geom: MolecularGeometry) -> None:
        """Set the geometry and invalidate the cached solution.

        Accepts either an (natoms, 3) block -- what the Libra adapter builds --
        or a flat array of dof values. Only the leading ``ndof`` entries are
        used, so a 1-D or 2-D model can ride inside a 3-coordinate atom.
        """
        coords = np.asarray(geom.coords_bohr, dtype=np.float64).reshape(-1)
        if coords.size < self.ndof:
            raise ValueError(
                f"Geometry supplies {coords.size} coordinates, "
                f"but the model needs {self.ndof}."
            )
        self._geom = geom
        self._coords = coords[: self.ndof].copy()
        self._state = None                       # force a recompute

    def get_geom(self) -> MolecularGeometry:
        if self._geom is None:
            raise ValueError("Geometry has not been set.")
        return self._geom

    def get_state(self) -> Optional[LVC_State]:
        return self._state

    def get_previous_state(self) -> Optional[LVC_State]:
        return self._previous_state

    def snapshot_state(self) -> None:
        """Save the current solution as the previous-geometry snapshot."""
        self._previous_state = pycopy.deepcopy(self._state)

    def copy(self) -> "LVC":
        """Independent clone, including the previous-geometry snapshot."""
        return pycopy.deepcopy(self)


# =============================================================================
# SECTION 4 -- ES_Strategy: the requested quantities
# =============================================================================

    def compute_H_el(self) -> np.ndarray:
        """Adiabatic energies, ascending.

        Returns a 1-D array of length nstates -- the convention the other
        implementations and the Libra adapter use, despite the shape comment
        on ``ES_Result.H_el``.
        """
        return self._solve().energies.copy()

    def compute_gradient(self, root: int = 0) -> np.ndarray:
        """Gradient of one adiabatic energy, shape (natoms, 3), Hartree/Bohr.

        Hellmann-Feynman on a real symmetric matrix:
        ``dE_k/dq_n = <k| dH/dq_n |k>``.
        """
        state = self._solve()
        if not 0 <= root < self.nstates:
            raise IndexError(
                f"Requested root {root}, but the model has {self.nstates} states."
            )
        vector = state.U[:, root]
        per_dof = np.einsum("i,nij,j->n", vector, state.dH_dia, vector)
        return self._pack(per_dof)

    def compute_all_gradients(self) -> list[np.ndarray]:
        """Gradients for every adiabatic state; one diagonalization for all."""
        state = self._solve()
        per_state = np.einsum("ik,nij,jk->kn", state.U, state.dH_dia, state.U)
        return [self._pack(per_state[k]) for k in range(self.nstates)]

    def compute_hessian(self, root: int = 0) -> np.ndarray:
        """Diabatic-limit Hessian, shape (3N, 3N), Hartree/Bohr^2.

        Exact only where the states do not mix (the second derivative of the
        adiabatic energy picks up a term from the rotation of ``U``); it is the
        curvature of the corresponding parabola, which is what a test of the
        Hessian plumbing needs.
        """
        self._solve()
        if not 0 <= root < self.nstates:
            raise IndexError(
                f"Requested root {root}, but the model has {self.nstates} states."
            )
        size = 3 * self.natoms
        hessian = np.zeros((size, size), dtype=np.float64)
        hessian[np.diag_indices(self.ndof)] = self._k[root]
        return hessian

    def compute_nac_vectors(self) -> np.ndarray:
        """NAC vectors, shape (nstates, nstates, natoms, 3), Bohr^-1.

        ``d_kl = <k| dH/dq |l> / (E_l - E_k)``, zero on the diagonal and where
        two states are degenerate to within 1e-12 Hartree.
        """
        state = self._solve()
        coupling = np.einsum("ik,nij,jl->kln", state.U, state.dH_dia, state.U)
        gaps = state.energies[None, :] - state.energies[:, None]     # E_l - E_k

        nac = np.zeros(
            (self.nstates, self.nstates, self.natoms, 3), dtype=np.float64
        )
        for k in range(self.nstates):
            for l in range(self.nstates):
                if k == l or abs(gaps[k, l]) < 1e-12:
                    continue
                nac[k, l] = self._pack(coupling[k, l] / gaps[k, l])
        return nac

    def compute_time_overlap(self, state1: object, state2: object) -> np.ndarray:
        """Adiabatic time-overlap matrix, shape (nstates, nstates).

        ``compute_result`` passes ``state1`` = current geometry and ``state2``
        = previous geometry. Libra's convention is
        ``S[i][j] = <phi_i(t)|phi_j(t+dt)>`` -- matching
        ``dyn_ham.cpp``, which builds it as
        ``basis_transform_prev.H() * basis_transform`` -- so this returns
        ``U_previous.T @ U_current``.

        Because the diabatic basis does not depend on geometry, this is exact.
        """
        if state1 is None or state2 is None:
            raise ValueError("Both the current and previous states are required.")
        return np.asarray(state2.U, dtype=np.float64).T @ np.asarray(
            state1.U, dtype=np.float64
        )

    # -- extras that the ABC does not require but are handy in tests ----------
    def compute_basis_transform(self) -> np.ndarray:
        """Diabatic-to-adiabatic transformation ``U`` at the current geometry."""
        return self._solve().U.copy()

    def compute_H_dia(self) -> np.ndarray:
        """Diabatic matrix at the current geometry, shape (nstates, nstates)."""
        return self._solve().H_dia.copy()


# =============================================================================
# SECTION 5 -- Interop: legacy parameter sets and the Libra adapter
# =============================================================================

    @classmethod
    def from_lvc_params(cls, params: dict, ndof: Optional[int] = None) -> "LVC":
        """Build the model from a legacy ``libra_py.models.LVC`` parameter dict.

        The legacy diabatic diagonal

            V_i(q) = Delta_i + sum_n [ 1/2 m_n w_n^2 q_n^2 + sqrt(m_n) d_i[n] q_n ]

        is the same parabola written about the origin instead of about its own
        minimum. Completing the square gives the exact reparametrization

            center_i[n]   = -d_i[n] / ( sqrt(m_n) * w_n^2 )
            energy_min_i  = Delta_i - sum_n d_i[n]^2 / ( 2 w_n^2 )

        and the legacy off-diagonal ``sum_n sqrt(m_n) coup[n] q_n`` becomes a
        linear coupling with slopes ``sqrt(m_n) * coup[n]``.

        So ``get_LVC_set1() / set2() / set3()`` (Fulvene, BMA, MIA) can be fed
        straight in, given a ``mass`` entry -- which those helpers do not
        provide, so supply one alongside them.
        """
        omega = np.asarray(params["omega"], dtype=np.float64)
        ndof = int(ndof) if ndof is not None else omega.size
        omega = omega[:ndof]

        mass = _broadcast(params.get("mass", 1.0), ndof, "mass")
        d1 = np.asarray(params["d1"], dtype=np.float64)[:ndof]
        d2 = np.asarray(params["d2"], dtype=np.float64)[:ndof]
        coup = np.asarray(params["coup"], dtype=np.float64)[:ndof]

        def parabola(delta: float, d: np.ndarray) -> Parabola:
            center = -d / (np.sqrt(mass) * omega**2)
            energy_min = float(delta) - float(np.sum(d**2 / (2.0 * omega**2)))
            return Parabola(
                energy_min=energy_min,
                center=tuple(center),
                omega=tuple(omega),
            )

        return cls(
            ndof=ndof,
            states=[
                parabola(params["Delta1"], d1),
                parabola(params["Delta2"], d2),
            ],
            coupling=Coupling.linear(slopes=tuple(np.sqrt(mass) * coup)),
            mass=tuple(mass),
        )

    def libra_model_params(self, dt: float = 41.0, **extra) -> dict:
        """A ready ``model_params`` dict for ``pyscf_compute_adi``.

        Fills in the keys the Libra adapter reads -- ``atom_labels`` sized so
        that ``3 * natoms >= ndof``, ``nstates``, ``dt``, and this instance as
        the strategy -- so a model run needs no hand-assembled dictionary.
        Any keyword given here overrides the defaults.
        """
        params = {
            "atom_labels": tuple("X" for _ in range(self.natoms)),
            "nstates": self.nstates,
            "dt": float(dt),
            "es_strategy": self,
            "gradient_state": "all",
            "time_overlap": True,
            "nacv": False,
        }
        params.update(extra)
        return params
