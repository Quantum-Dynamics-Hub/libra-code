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
.. module:: pyscf.implementations.df_casscf
   :platform: Unix, Windows
   :synopsis: PySCF density-fitted CASSCF implementation for Libra ES strategy.
.. moduleauthor::
       Jieyang Gu <jieyanggu792@gmail.com>

"""

from __future__ import annotations
import copy as pycopy
from dataclasses import dataclass
import sys
from pathlib import Path
from typing import Any, List, Optional, Sequence, Tuple, Union

import numpy as np
from pyscf import fci, gto, mcscf, scf
from pyscf.grad.deriv_eri import DerivativeERICache
from libra_py.packages.pyscf.interfaces import ES_Strategy, ES_Request, MolecularGeometry

BOHR_TO_ANG = 0.529177210903


@dataclass
class DF_CASSCF_States:
    mol: Optional[Any] = None
    mf: Optional[Any] = None
    mc: Optional[Any] = None


@dataclass
class gradient_nacv_shared:
    #put shared variables reused by gradient and nacv within a geom here like eris and orbital response LHS
    # DF has no NAC path (see the class docstring), so only the gradient reads
    # these.  Keeping them out of DF_CASSCF_States is also what lets the state
    # be deep-copied at all: the DF _ERIS holds its ppaa/papa blocks in a
    # lib.H5TmpFile, and deepcopy raises "h5py objects cannot be pickled".
    eris: Optional[Any] = None           # mc.ao2mo(mc.mo_coeff), a pyscf.mcscf.df._ERIS
    cphf_lhs: Optional[Any] = None       # (Aop, Adiag): the CP-MCSCF left-hand side
    cphf_precond: Optional[Any] = None   # and its preconditioner
    AO_derivative_integrals: Optional[Any] = None  # provider for (nabla i,j|k,l)


class DF_CASSCF(ES_Strategy):
    """PySCF-based density-fitted CASSCF backend for the universal ES interface.

    Identical in sequencing to the conventional :class:`CASSCF` backend, with
    the two-electron integrals resolved by the identity throughout: the SCF
    guess is a DF-RHF/ROHF, and the MCSCF object carries the ``_DFCASSCF``
    mixin, so the JK build, the active-space integral transform and the
    nuclear gradients all go through ``mc.with_df`` rather than the exact ERI.
    The same ``auxbasis`` is handed to the SCF and to the MCSCF so that
    ``mcscf.df.density_fit`` reuses the SCF's ``with_df`` object instead of
    building a second one.

    NAC vectors are not implemented: PySCF has no density-fitted analogue of
    pyscf.nac.sacasscf, and dispatching to the exact-ERI one on top of DF
    orbitals would mix the two integral treatments inside a single Lagrangian.
    """

    # 1) when a geom is set the DF-HF is run; the MF is written to the member
    #    attribute `_mf`
    # 2) when the CASSCF energy is requested for a certain root, the DF-CASSCF
    #    is run (number of roots requested = self._nroots). The resulting MC
    #    object is stored in the member attribute `mc`, which also contains
    #    energies of other states.

    def __init__(
        self,
        mol: Optional[Any] = None,
        norbcas: int = 0,
        nelecas: int = 0,
        nroots: int = 1,
        basis: str = "sto-3g",
        unit: str = "Bohr",
        charge: int = 0,
        cas_list: Optional[List[int]] = None,
        spin_multiplicity: int = 1,
        auxbasis: Optional[Any] = None,
    ) -> None:
        #setting up the initial state
        self._mol: Optional[Any] = mol
        self._norbcas: int = norbcas
        self._nelecas: Union[int, Tuple[int, int]] = nelecas
        self._nroots: int = nroots    #can only be >1
        self._basis: str = basis
        self._unit: str = unit  #default to Bohr, but can be set to Angstrom
        self._charge: int = int(charge)
        self._cas_list: Optional[List[int]] = cas_list
        self._spin_multiplicity = int(spin_multiplicity)
        # None lets PySCF pick the fitting basis paired with self._basis
        # (weigend+etb when the orbital basis is a dict rather than a string).
        self._auxbasis: Optional[Any] = auxbasis
        if self._spin_multiplicity < 1:
            raise ValueError("spin_multiplicity must be a positive integer")
        self._spin = self._spin_multiplicity - 1
        self._geom: Optional[MolecularGeometry] = None
        self._state: DF_CASSCF_States | None = None
        self._previous_state: DF_CASSCF_States | None = None
        # rebuilt with mc: one per wavefunction, never snapshotted
        self._shared: gradient_nacv_shared = gradient_nacv_shared()
        self.eris_reuse: bool = False  
        self.cphf_reuse: bool = False  # does not work correctly
        self.deriv_eri_reuse: bool = False  

    @property
    def nroots(self) -> int:
        """Number of roots this strategy solves for.

        Fixed at construction, and the single source of truth for the state
        count: the adapter builds ``ES_Request.n_singlets`` from it rather than
        from a separate ``model_params["nstates"]`` entry.
        """
        return int(self._nroots)

    @property
    def spin_multiplicity(self) -> int:
        """Spin multiplicity ``2*S+1`` represented by this strategy."""
        return self._spin_multiplicity

    def _n_total(self) -> int:
        """Number of electronic states requested for the current calculation."""
        return int(self._nroots)

    def _make_scf(self) -> Any:
        """Build a fresh DF-RHF/ROHF object at ``self._mol``."""
        scf_cls = scf.RHF if self._spin == 0 else scf.ROHF
        return scf_cls(self._mol).density_fit(auxbasis=self._auxbasis)

    def _run_scf(self) -> None:
        """Run restricted DF-HF at ``self._mol``.

        When a previous geometry exists, the SCF is restarted from the previous
        HF 1-electron density matrix (same AO basis and ordering), which is the
        standard warm-start used in dynamics.  Falls back to an MO-based guess
        and then to the default SCF guess.
        """
        if self._state is None:
            self._state = DF_CASSCF_States(mol=self._mol, mf=None, mc=None)

        prev_state = self.get_previous_state()
        if prev_state is None or prev_state.mf is None:
            self._state.mf = self._make_scf().run(verbose=0)
            return

        try:
            prev_dm = prev_state.mf.make_rdm1()
        except Exception:
            prev_dm = None

        if prev_dm is not None:
            try:
                self._state.mf = self._make_scf().run(dm0=prev_dm, verbose=0)
                return
            except Exception:
                pass

        try:
            mf = self._make_scf()
            mf.init_guess_by_mo(prev_state.mf.mo_coeff)
            self._state.mf = mf.run(verbose=0)
        except Exception:
            self._state.mf = self._make_scf().run(verbose=0)

    def set_geom(self, geom: MolecularGeometry) -> None:
        self._geom = geom
        self._previous_state = pycopy.deepcopy(self._state)

        charge: int = self._charge
        coords = getattr(geom, "coords_bohr", None)
        if coords is None:
            coords = np.asarray(getattr(geom, "coords_angstrom"), dtype=float) / BOHR_TO_ANG
        else:
            coords = np.asarray(coords, dtype=float)

        self._mol = gto.M(
            atom=";".join(
                f"{label} {coord[0]} {coord[1]} {coord[2]}"
                for label, coord in zip(geom.atom_labels, coords)
            ),
            basis=self._basis,
            unit=self._unit,
            charge=charge,
            spin=self._spin,
        )

        self._run_scf()

        # Write the current SCF result back into the state object.
        if self._state is None:
            self._state = DF_CASSCF_States(mol=self._mol, mf=None, mc=None)
        self._state.mol = self._mol
        self._state.mc = None
        self._shared = gradient_nacv_shared()   # new geometry -> everything shared is stale

    def get_geom(self) -> MolecularGeometry:
        if self._geom is None:
            raise ValueError("Geometry has not been set.")
        return self._geom

    def get_state(self) -> Optional[DF_CASSCF_States]:
        return self._state

    def get_previous_state(self) -> Optional[DF_CASSCF_States]:
        return self._previous_state

    def snapshot_state(self) -> None:
        """Save the current state of the calculation to a previous-state snapshot."""
        self._previous_state = pycopy.deepcopy(self._state)

    def copy(self) -> "DF_CASSCF":
        """Return an independent snapshot of the full strategy state."""
        # The clone starts with empty shared intermediates: they are rebuilt on
        # demand, and the DF eris keeps its ppaa/papa in a lib.H5TmpFile, which
        # deepcopy refuses ("h5py objects cannot be pickled").
        shared, self._shared = self._shared, gradient_nacv_shared()
        try:
            return pycopy.deepcopy(self)
        finally:
            self._shared = shared

    def compute_H_el(self) -> np.ndarray:
        restarting = (
            self._previous_state is not None and self._previous_state.mc is not None
        )

        # .density_fit() attaches the _DFCASSCF mixin.  Because the SCF was
        # built with the same auxbasis, mcscf.df.density_fit adopts
        # mf.with_df rather than constructing a second fitting object, so the
        # SCF and the MCSCF share one set of 3-center integrals.
        mc = mcscf.CASSCF(self._state.mf, self._norbcas, self._nelecas)
        mc = mc.density_fit(auxbasis=self._auxbasis)
        if self._spin == 0:
            mc.fcisolver = fci.direct_spin0.FCI(self._state.mol)
        else:
            mc.fcisolver = fci.direct_spin1.FCI(self._state.mol)
            target_s = 0.5 * self._spin
            mc.fcisolver = fci.addons.fix_spin_(
                mc.fcisolver, ss=target_s * (target_s + 1.0)
            )
        mc.fcisolver.nroots = self._nroots
        if self._nroots > 1:
            mc = mc.state_average_([1.0 / self._nroots] * self._nroots)  # equal weights as default; required for gradients in pyscf

        if restarting:
            # Orbitals converged at the PREVIOUS geometry are orthonormal in
            # that geometry's AO metric, not in this one.  Handing them to
            # mc.kernel unprojected leaves mc.mo_coeff non-orthonormal here,
            # because the orbital optimizer only applies unitary rotations and
            # so preserves whatever non-orthonormality it was given.  When the
            # state-averaged surface is flat -- e.g. HeH+ CAS(2,2) with 3 roots,
            # where the roots span the complete active space and the SA energy
            # is a basis-independent trace -- the optimizer returns the guess
            # bit-for-bit and the error survives untouched into mc.mo_coeff.
            #
            # project_init_guess re-orthonormalizes against the current
            # molecule (SVD per orbital subspace) while keeping the active
            # space in the canonical ncore:ncore+ncas block, which is what
            # makes it safe to skip cas_list below.
            mocoeff = mcscf.project_init_guess(
                mc,
                self._previous_state.mc.mo_coeff,
                prev_mol=self._previous_state.mol,
            )
        else:
            mocoeff = self._state.mf.mo_coeff

            # cas_list indexes the *HF* orbitals, so it only applies to an HF
            # starting guess.  A restart already carries its active space in
            # the canonical block; re-sorting would pull a different set of
            # columns entirely -- e.g. cas_list = [4, 7, 11, 14, 17] selects
            # 0-indexed 3, 6, 10, 13, 16 while the active space already sits at
            # 5-9.  CASSCF often re-converges to the same energies from that
            # scrambled guess, which is what makes the bug easy to miss.
            if self._cas_list is not None:
                mocoeff = mcscf.sort_mo(mc, mocoeff, self._cas_list)

        mc.kernel(mocoeff)

        self._state.mc = mc
        self._shared = gradient_nacv_shared()   # new wavefunction -> everything shared is stale

        #write the energies to the H_el attribute of the state object and return the energies as a numpy array
        energies = getattr(mc, "e_states", None)
        if energies is None:
            energies = [mc.e_tot]
        self._state.H_el = np.asarray(energies, dtype=np.float64)
        return self._state.H_el

    #gradient

    #helpers for the parts every root shares

    def compute_eris(self) -> Any:
        """DF MO-basis ERI at the converged CASSCF orbitals.

        Depends on mo_coeff alone, so one transform serves every gradient root
        of this wavefunction.  mc.ao2mo returns a pyscf.mcscf.df._ERIS here,
        built from the 3-center (P|pq) integrals rather than the full 4-index
        tensor.
        """
        mc = self._state.mc
        return mc.ao2mo(mc.mo_coeff)

    def compute_deriv_eri(self) -> Any:
        """Provider for the DF derivative integrals, shared across roots.

        Density fitting differentiates the 3-center (i,j|P) rather than forming
        (nabla i,j|k,l), so the provider serves int3c2e_ip1/int3c2e_ip2 here.
        They depend on the geometry, basis and auxbasis alone -- not on the root
        or the CI vectors -- so one evaluation feeds every gradient root.
        """
        if self._shared.AO_derivative_integrals is None:
            self._shared.AO_derivative_integrals = DerivativeERICache()
        return self._shared.AO_derivative_integrals

    def compute_cphf_lhs(self, solver: Any) -> Any:
        """Build the CP-MCSCF left-hand side once, then feed it to every solve.

        Every root solves A z = -b with the same A = d2 E_SA/dp dq: PySCF builds
        it from make_fcasscf_sa at the converged (mo, ci), and the projection it
        wraps A in ignores the state.  Only b differs.  pyscf.df.grad.sacasscf
        inherits these hooks from the exact-ERI solver, so the DF path shares
        them the same way.
        """
        shared = self._shared
        build_lhs, build_precond = solver.get_Aop_Adiag, solver.get_lagrange_precond

        def get_Aop_Adiag(**kwargs):
            if shared.cphf_lhs is None:
                shared.cphf_lhs = build_lhs(**kwargs)
            return shared.cphf_lhs

        def get_lagrange_precond(Adiag, level_shift=None, **kwargs):
            if shared.cphf_precond is None:
                shared.cphf_precond = build_precond(Adiag, level_shift=level_shift, **kwargs)
            return shared.cphf_precond

        solver.get_Aop_Adiag = get_Aop_Adiag
        solver.get_lagrange_precond = get_lagrange_precond
        return solver

    def compute_gradient(self, roots: Sequence[int]) -> list[np.ndarray]:
        """Nuclear gradients for ``roots``, returned in the order requested.

        The state-averaged branch goes through pyscf.df.grad.sacasscf, which is
        the DF counterpart of the exact-ERI Lagrangian and includes the
        auxiliary-basis response term (``auxbasis_response=True`` by default),
        so the gradient is the exact derivative of the DF energy rather than an
        approximation to the conventional one.

        Note: PySCF's DF gradient JK build divides by the occupancy of each
        density matrix it is handed, so an active space with no closed-shell
        core (ncore == 0) raises ZeroDivisionError inside
        pyscf.df.grad.rhf.get_jk.  That is an upstream limitation, not one
        introduced here; use the conventional CASSCF backend for such cases.
        """
        mc = self._state.mc
        roots = list(roots)
        n_avail = self._nroots
        for root in roots:
            if not 0 <= root < n_avail:
                raise IndexError(
                    f"Requested root {root}, but only {n_avail} roots are available."
                )

        # 1. compute the intermediate shared by every root
        if self._nroots == 1:
            grad = mc.nuc_grad_method()
            if self.deriv_eri_reuse:
                grad.deriv_eri = self.compute_deriv_eri()
            return [
                np.asarray(grad.kernel(), dtype=np.float64)
                for _ in roots
            ]

        grad = mc.nuc_grad_method(state=0)
        if self.deriv_eri_reuse:
            grad.deriv_eri = self.compute_deriv_eri()
        if self.cphf_reuse:
            self.compute_cphf_lhs(grad)

        # 2. per root: build a shared part on first use, then read it back
        gradients = []
        for root in roots:
            eris = None
            if self.eris_reuse:
                if self._shared.eris is None:
                    self._shared.eris = self.compute_eris()
                eris = self._shared.eris
            gradients.append(
                np.asarray(grad.kernel(state=root, eris=eris), dtype=np.float64)
            )
        return gradients


    # time overlap

    def _compute_ao_overlap(self, prev_state: DF_CASSCF_States, curr_state: DF_CASSCF_States) -> np.ndarray:
        """Compute the AO overlap matrix between the previous and current geometries."""
        prev_mol = prev_state.mol
        curr_mol = curr_state.mol
        if prev_mol is None or curr_mol is None:
            raise ValueError("Both previous and current molecule objects are required.")
        return gto.intor_cross("int1e_ovlp", prev_mol, curr_mol)

    def _make_casci(self, state: DF_CASSCF_States, nroots: int) -> Any:
        """Build a DF-CASCI solver on top of ``state``'s DF-HF orbitals."""
        casci = mcscf.CASCI(state.mf, self._norbcas, self._nelecas)
        casci = casci.density_fit(auxbasis=self._auxbasis)
        solver_cls = fci.direct_spin0.FCISolver if self._spin == 0 else fci.direct_spin1.FCISolver
        casci.fcisolver = solver_cls(state.mol)
        if self._spin != 0:
            target_s = 0.5 * self._spin
            casci.fcisolver = fci.addons.fix_spin_(
                casci.fcisolver, ss=target_s * (target_s + 1.0)
            )
        casci.fcisolver.nroots = nroots
        return casci

    def compute_time_overlap(self, state1: object, state2: object) -> np.ndarray:
        """Compute the adiabatic time-overlap matrix between two state snapshots.

        state1: the current-state snapshot (``DF_CASSCF_States``).
        state2: the previous-state snapshot (``DF_CASSCF_States``).

        Returns
        -------
        np.ndarray
            Shape ``(n_total, n_total)``.
        """
        if not isinstance(state1, DF_CASSCF_States) or not isinstance(state2, DF_CASSCF_States):
            raise TypeError(
                f"state1/state2 must be DF_CASSCF_States, got "
                f"{type(state1).__name__} / {type(state2).__name__}."
            )
        prev_state = state2  # previous geometry
        curr_state = state1  # current geometry

        if prev_state.mf is None or curr_state.mf is None:
            raise ValueError("HF states are required for time-overlap computation.")
        if prev_state.mol is None or curr_state.mol is None:
            raise ValueError("Molecule states are required for time-overlap computation.")

        nroots = self._n_total()
        if nroots <= 0:
            raise ValueError(f"nroots must be positive, got {nroots}")

        ao_overlap = self._compute_ao_overlap(prev_state, curr_state)

        # Always build CASCI roots and matching active-space MO blocks for both
        # geometries.  State-averaged CASSCF `mc.ci` representations can differ
        # (and be incompatible with `fci.addons.overlap`'s expected transform),
        # so recomputing CASCI ensures consistent CI + one-particle overlap.
        # The DF CASCI returns get_h2cas in the compact (npair, npair) layout;
        # the FCI solvers restore it internally, so it is passed through as is.
        prev_casci = self._make_casci(prev_state, nroots)
        h1prev, _ = prev_casci.get_h1eff(prev_casci.mo_coeff)
        h2prev = prev_casci.get_h2cas(prev_casci.mo_coeff)
        _, prev_roots_raw = prev_casci.fcisolver.kernel(
            h1prev,
            h2prev,
            prev_casci.ncas,
            prev_casci.nelecas,
            nroots=nroots,
        )

        curr_casci = self._make_casci(curr_state, nroots)
        h1curr, _ = curr_casci.get_h1eff(curr_casci.mo_coeff)
        h2curr = curr_casci.get_h2cas(curr_casci.mo_coeff)
        _, curr_roots_raw = curr_casci.fcisolver.kernel(
            h1curr,
            h2curr,
            curr_casci.ncas,
            curr_casci.nelecas,
            nroots=nroots,
        )

        if nroots == 1:
            prev_roots = [np.asarray(prev_roots_raw)]
            curr_roots = [np.asarray(curr_roots_raw)]
        else:
            prev_roots = [np.asarray(vec) for vec in prev_roots_raw[:nroots]]
            curr_roots = [np.asarray(vec) for vec in curr_roots_raw[:nroots]]
        prev_act = np.asarray(prev_casci.mo_coeff)[
            :, prev_casci.ncore : prev_casci.ncore + prev_casci.ncas
        ]
        curr_act = np.asarray(curr_casci.mo_coeff)[
            :, curr_casci.ncore : curr_casci.ncore + curr_casci.ncas
        ]
        s12_mo = prev_act.T.conj() @ ao_overlap @ curr_act

        overlap = np.zeros((nroots, nroots), dtype=float)
        for i in range(nroots):
            for j in range(nroots):
                overlap[i, j] = fci.addons.overlap(
                    prev_roots[i],
                    curr_roots[j],
                    self._norbcas,
                    self._nelecas,
                    s=s12_mo,
                )

        overlap = np.asarray(np.real_if_close(overlap))
        for j in range(nroots):
            if overlap[j, j] < 0:
                overlap[:, j] = -overlap[:, j]
        return overlap

if __name__ == "__main__":

    def _smoke_test_df_casscf_lih_sequence() -> None:
        """Smoke test for LiH DF-CASSCF state sequencing and time-overlap.

        LiH rather than HeH+ because the DF gradient needs ncore > 0; see the
        note on compute_gradient.
        """
        geom1 = MolecularGeometry(
            atom_labels=['Li', 'H'],
            coords_bohr=np.array([
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 3.0],
            ], dtype=np.float64),
        )
        geom2 = MolecularGeometry(
            atom_labels=['Li', 'H'],
            coords_bohr=np.array([
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 3.2],
            ], dtype=np.float64),
        )

        df_casscf = DF_CASSCF(norbcas=2, nelecas=2, nroots=2, basis='sto-3g', unit='Bohr')

        request = ES_Request(
            n_singlets=2,
            gradient_state='all',
            hessian_state=None,
            nacv=False,
            time_overlap=True,
        )

        result1 = df_casscf.compute_result(geom1, request)
        print("\nResult 1 H_el:", result1.H_el)
        print("Result 1 gradients:\n", result1.gradients)
        print("Result 1 time_overlap:", result1.time_overlap)

        result2 = df_casscf.compute_result(geom2, request)
        print("\nResult 2 H_el:", result2.H_el)
        print("Result 2 gradients:\n", result2.gradients)
        print("Result 2 time_overlap:\n", result2.time_overlap)

    def _smoke_test_df_casscf_nacv_refused() -> None:
        """NACVs are deliberately unimplemented; the request must say so."""
        geom = MolecularGeometry(
            atom_labels=['Li', 'H'],
            coords_bohr=np.array([
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 3.0],
            ], dtype=np.float64),
        )
        df_casscf = DF_CASSCF(norbcas=2, nelecas=2, nroots=2, basis='sto-3g', unit='Bohr')
        request = ES_Request(n_singlets=2, gradient_state=None, nacv=True)
        try:
            df_casscf.compute_result(geom, request)
        except NotImplementedError as exc:
            print("\nNACV request refused as expected:\n ", exc)
        else:
            raise AssertionError("expected NotImplementedError for nacv=True")

    def _smoke_test_df_casscf_c2h4_eri_reuse_timing(basis: str = '6-31g*', repeats: int = 2) -> None:
        """Time the three SA-3 gradient calls of C2H4 CAS(4,4), eris_reuse on vs off.

        Only compute_gradient([0, 1, 2]) is inside the timer: DF-SCF and
        DF-CASSCF are converged once, beforehand, and both variants
        differentiate that same wavefunction.
          shared -- eris_reuse=True:  one DF mc.ao2mo, then three kernels fed from it
          naive  -- eris_reuse=False: three kernels, each running its own DF ao2mo
        CAS(4,4) rather than CAS(2,2): three singlet roots of (2e,2o) are the
        complete space, whose state-averaged surface is flat.  The variants
        alternate and the fastest of ``repeats`` is kept, so neither pays alone
        for first-call warm-up.
        """
        import time
        from pyscf import lib
        from pyscf.mcscf import df as mcscf_df

        geom = MolecularGeometry(
            atom_labels=['C', 'C', 'H', 'H', 'H', 'H'],
            coords_bohr=np.array([
                [0.0000,  0.0000,  0.6695],
                [0.0000,  0.0000, -0.6695],
                [0.0000,  0.9289,  1.2321],
                [0.0000, -0.9289,  1.2321],
                [0.0000,  0.9289, -1.2321],
                [0.0000, -0.9289, -1.2321],
            ], dtype=np.float64) / BOHR_TO_ANG,
        )
        roots = [0, 1, 2]
        df_casscf = DF_CASSCF(norbcas=4, nelecas=4, nroots=3, basis=basis, unit='Bohr')
        df_casscf.set_geom(geom)
        e_states = df_casscf.compute_H_el()

        # Count the DF AO->MO transforms that happen inside the timed calls.
        builds = {'n': 0, 't': 0.0}
        eris_init = mcscf_df._ERIS.__init__
        def counted_init(self, *args, **kwargs):
            t0 = time.perf_counter()
            eris_init(self, *args, **kwargs)
            builds['n'] += 1
            builds['t'] += time.perf_counter() - t0
        mcscf_df._ERIS.__init__ = counted_init

        best, grads = {}, {}
        try:
            for _ in range(repeats):
                for tag, reuse in (('naive', False), ('shared', True)):
                    df_casscf.eris_reuse = reuse
                    df_casscf._shared = gradient_nacv_shared()   # shared pays for its one transform
                    builds['n'], builds['t'] = 0, 0.0
                    t0 = time.perf_counter()
                    grads[tag] = df_casscf.compute_gradient(roots)
                    wall = time.perf_counter() - t0
                    if tag not in best or wall < best[tag][0]:
                        best[tag] = (wall, builds['n'], builds['t'])
        finally:
            mcscf_df._ERIS.__init__ = eris_init
            df_casscf.eris_reuse = True

        err = max(float(np.max(np.abs(a - b))) for a, b in zip(grads['naive'], grads['shared']))
        t_naive, t_shared = best['naive'][0], best['shared'][0]
        print(f"\nC2H4 DF-CAS(4,4) SA-3  basis={basis}  nao={df_casscf.get_state().mol.nao}"
              f"  naux={df_casscf.get_state().mf.with_df.get_naoaux()}"
              f"  threads={lib.num_threads()}  -- compute_gradient({roots}) only, best of {repeats}")
        print(f"  E_states = {e_states}")
        print(f"  {'variant':<8}{'wall/s':>9}{'ao2mo builds':>14}{'in ao2mo/s':>12}")
        for tag in ('naive', 'shared'):
            wall, n, t = best[tag]
            print(f"  {tag:<8}{wall:>9.2f}{n:>14d}{t:>12.2f}")
        print(f"  saved {t_naive - t_shared:.2f} s  ({100 * (t_naive - t_shared) / t_naive:.1f}%),"
              f"  speedup {t_naive / t_shared:.2f}x,  max |g_naive - g_shared| = {err:.1e}")
        assert err < 1e-8, f"shared and naive gradients differ by {err:.2e}"

    _smoke_test_df_casscf_lih_sequence()
    _smoke_test_df_casscf_c2h4_eri_reuse_timing()
    _smoke_test_df_casscf_nacv_refused()    
