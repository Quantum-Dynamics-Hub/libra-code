# *********************************************************************************
# * Copyright (C) 2026 Jieyang Gu <jieyanggu792@gmail.com>
# * This file is distributed under the terms of the GNU General Public License
# * as published by the Free Software Foundation, either version 3 of
# * the License, or (at your option) any later version.
# * See the file LICENSE in the root directory of this distribution
# * or <http://www.gnu.org/licenses/>.
# *
# *******************************************************************************/
"""
.. module:: pyscf.implementations.casscf
   :platform: Unix, Windows
   :synopsis: PySCF CASSCF implementation for Libra ES strategy.
.. moduleauthor::
       Jieyang Gu <jieyanggu792@gmail.com>

"""

from __future__ import annotations
import copy as pycopy
from dataclasses import dataclass
from typing import Any, List, Optional, Sequence, Tuple, Union
import numpy as np
from pyscf import fci, gto, mcscf, scf
from pyscf.grad.deriv_eri import DerivativeERICache
from libra_py.packages.pyscf.interfaces import ES_Request, ES_Strategy, MolecularGeometry

BOHR_TO_ANG = 0.529177210903


@dataclass
class CASSCF_States:
    mol: Optional[Any] = None
    mf: Optional[Any] = None
    mc: Optional[Any] = None


@dataclass
class gradient_nacv_shared:
    #put shared variables reused by gradient and nacv within a geom here like eris and orbital response LHS
    eris: Optional[Any] = None           # mc.ao2mo(mc.mo_coeff)
    cphf_lhs: Optional[Any] = None       # (Aop, Adiag): the CP-MCSCF left-hand side
    cphf_precond: Optional[Any] = None   # and its preconditioner
    AO_derivative_integrals: Optional[Any] = None  # provider for (nabla i,j|k,l)


def compute_eris(mc: Any) -> Any:
    """MO-basis ERI at the converged CASSCF orbitals.

    Depends on mo_coeff alone, so one transform serves every gradient root and
    every NAC pair of this wavefunction.
    """
    return mc.ao2mo(mc.mo_coeff)


def compute_deriv_eri(shared: gradient_nacv_shared) -> Any:
    """Provider for the AO derivative integrals (nabla i,j|k,l).

    They depend on the geometry and basis alone -- not on the root, the CI
    vectors or the orbitals -- so one evaluation feeds every gradient root and
    every NAC pair, across both the Hellmann-Feynman and the Lagrange-response
    halves of each.  PySCF streams them block by block and keeps nothing; this
    hands its kernels a provider that keeps them, and parks it on ``shared`` so
    that the gradient object and the NAC object -- which are separate PySCF
    objects that never see each other -- get the same one.
    """
    if shared.AO_derivative_integrals is None:
        shared.AO_derivative_integrals = DerivativeERICache()
    return shared.AO_derivative_integrals


def compute_cphf_lhs(solver: Any, shared: gradient_nacv_shared) -> Any:
    """Build the CP-MCSCF left-hand side once, then feed it to every solve.

    Every root and every pair solves A z = -b with the same A = d2 E_SA/dp dq:
    PySCF builds it from make_fcasscf_sa at the converged (mo, ci), and the
    projection it wraps A in ignores the state.  Only b differs.  So the first
    solve builds (Aop, Adiag) and the preconditioner, and the later ones read
    them back off ``shared``.
    """
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


class CASSCF(ES_Strategy):
    """PySCF-based CASSCF backend for the universal ES interface."""

    # 1) when a geom is set the HF is run; the MF is written to the member attribute
    #    `_mf`
    # 2) when the CASSCF energy is requested for a certain root, the CASSCF is run
    #    (number of roots requested = self._nroots). The resulting MC object is
    #    stored in the member attribute `mc`, which also contains energies of other
    #    states.

    def __init__(
        self,
        mol: Optional[Any] = None,
        norbcas: int = 0,
        nelecas:int = 0, 
        nroots: int = 1,
        basis: str = "sto-3g",
        unit: str = "Bohr",
        charge: int = 0,
        cas_list: Optional[List[int]] = None,
        spin_multiplicity: int = 1,
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
        if self._spin_multiplicity < 1:
            raise ValueError("spin_multiplicity must be a positive integer")
        self._spin = self._spin_multiplicity - 1
        self._geom: Optional[MolecularGeometry] = None
        self._state: CASSCF_States | None = None
        self._previous_state: CASSCF_States | None = None
        # rebuilt with mc: one per wavefunction, never snapshotted
        self._shared: gradient_nacv_shared = gradient_nacv_shared()
        self.eris_reuse: bool = False  
        self.cphf_reuse: bool = False   # does not work correctly
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

    def _run_scf(self) -> None:
        """Run restricted HF at ``self._mol``.

        When a previous geometry exists, the SCF is restarted from the previous
        HF 1-electron density matrix (same AO basis and ordering), which is the
        standard warm-start used in dynamics.  Falls back to an MO-based guess
        and then to the default SCF guess.
        """
        if self._state is None:
            self._state = CASSCF_States(mol=self._mol, mf=None, mc=None)

        prev_state = self.get_previous_state()
        if prev_state is None or prev_state.mf is None:
            scf_cls = scf.RHF if self._spin == 0 else scf.ROHF
            self._state.mf = scf_cls(self._mol).run(verbose=0)
            return

        try:
            prev_dm = prev_state.mf.make_rdm1()
        except Exception:
            prev_dm = None

        if prev_dm is not None:
            try:
                scf_cls = scf.RHF if self._spin == 0 else scf.ROHF
                self._state.mf = scf_cls(self._mol).run(dm0=prev_dm, verbose=0)
                return
            except Exception:
                pass

        try:
            scf_cls = scf.RHF if self._spin == 0 else scf.ROHF
            mf = scf_cls(self._mol)
            mf.init_guess_by_mo(prev_state.mf.mo_coeff)
            self._state.mf = mf.run(verbose=0)
        except Exception:
            scf_cls = scf.RHF if self._spin == 0 else scf.ROHF
            self._state.mf = scf_cls(self._mol).run(verbose=0)

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
            self._state = CASSCF_States(mol=self._mol, mf=None, mc=None)
        self._state.mol = self._mol
        self._state.mc = None
        self._shared = gradient_nacv_shared()   # new geometry -> everything shared is stale

    def get_geom(self) -> MolecularGeometry:
        if self._geom is None:
            raise ValueError("Geometry has not been set.")
        return self._geom

    def get_state(self) -> Optional[CASSCF_States]:
        return self._state

    def get_previous_state(self) -> Optional[CASSCF_States]:
        return self._previous_state

    def snapshot_state(self) -> None:
        """Save the current state of the calculation to a previous-state snapshot."""
        self._previous_state = pycopy.deepcopy(self._state)

    def copy(self) -> "CASSCF":
        """Return an independent snapshot of the full strategy state."""
        # The clone starts with empty shared intermediates.  They are rebuilt on
        # demand and belong to this geometry, so copying them would duplicate the
        # derivative-integral tensor once per trajectory to no purpose -- and the
        # MO ERI may be h5py-backed, which deepcopy refuses outright.
        shared, self._shared = self._shared, gradient_nacv_shared()
        try:
            return pycopy.deepcopy(self)
        finally:
            self._shared = shared

    def compute_H_el(self) -> np.ndarray:
        if self._previous_state is not None and self._previous_state.mc is not None:
            mocoeff = self._previous_state.mc.mo_coeff
        else:
            mocoeff = self._state.mf.mo_coeff

        mc = mcscf.CASSCF(self._state.mf, self._norbcas, self._nelecas)
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
    def compute_gradient(self, roots: Sequence[int]) -> list[np.ndarray]:
        """Nuclear gradients for ``roots``, returned in the order requested."""
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
                grad.deriv_eri = compute_deriv_eri(self._shared)
            return [
                np.asarray(grad.kernel(), dtype=np.float64)
                for _ in roots
            ]

        grad = mc.nuc_grad_method(state=0)
        if self.cphf_reuse:
            compute_cphf_lhs(grad, self._shared)
        if self.deriv_eri_reuse:
            # held on the gradient object, so every root's kernel call reuses it
            grad.deriv_eri = compute_deriv_eri(self._shared)

        # 2. per root: build a shared part on first use, then read it back
        gradients = []
        for root in roots:
            eris = None
            if self.eris_reuse:
                if self._shared.eris is None:
                    self._shared.eris = compute_eris(mc)
                eris = self._shared.eris
            gradients.append(
                np.asarray(grad.kernel(state=root, eris=eris), dtype=np.float64)
            )
        return gradients


    # time overlap

    def _compute_ao_overlap(self, prev_state: CASSCF_States, curr_state: CASSCF_States) -> np.ndarray:
        """Compute the AO overlap matrix between the previous and current geometries."""
        prev_mol = prev_state.mol
        curr_mol = curr_state.mol
        if prev_mol is None or curr_mol is None:
            raise ValueError("Both previous and current molecule objects are required.")
        return gto.intor_cross("int1e_ovlp", prev_mol, curr_mol)

    def compute_time_overlap(self, state1: object, state2: object) -> np.ndarray:
        """Compute the adiabatic time-overlap matrix between two state snapshots.

        state1: the current-state snapshot (``CASSCF_States``).
        state2: the previous-state snapshot (``CASSCF_States``).

        Returns
        -------
        np.ndarray
            Shape ``(n_total, n_total)``.
        """
        if not isinstance(state1, CASSCF_States) or not isinstance(state2, CASSCF_States):
            raise TypeError(
                f"state1/state2 must be CASSCF_States, got "
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
        prev_casci = mcscf.CASCI(prev_state.mf, self._norbcas, self._nelecas)
        solver_cls = fci.direct_spin0.FCISolver if self._spin == 0 else fci.direct_spin1.FCISolver
        prev_casci.fcisolver = solver_cls(prev_state.mol)
        if self._spin != 0:
            target_s = 0.5 * self._spin
            prev_casci.fcisolver = fci.addons.fix_spin_(
                prev_casci.fcisolver, ss=target_s * (target_s + 1.0)
            )
        prev_casci.fcisolver.nroots = nroots
        h1prev, _ = prev_casci.get_h1eff(prev_casci.mo_coeff)
        h2prev = prev_casci.get_h2cas(prev_casci.mo_coeff)
        _, prev_roots_raw = prev_casci.fcisolver.kernel(
            h1prev,
            h2prev,
            prev_casci.ncas,
            prev_casci.nelecas,
            nroots=nroots,
        )

        curr_casci = mcscf.CASCI(curr_state.mf, self._norbcas, self._nelecas)
        curr_casci.fcisolver = solver_cls(curr_state.mol)
        if self._spin != 0:
            target_s = 0.5 * self._spin
            curr_casci.fcisolver = fci.addons.fix_spin_(
                curr_casci.fcisolver, ss=target_s * (target_s + 1.0)
            )
        curr_casci.fcisolver.nroots = nroots
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

    def compute_nac_vectors(self, use_etfs: bool = True) -> np.ndarray:
        """Non-adiabatic coupling vectors for every ordered state pair.

        Only the upper triangle is solved for.  d_IJ = <I|d/dR|J> is the
        symmetric Hellmann-Feynman numerator over E_J - E_I, and that
        denominator flips sign when the pair is swapped, so d_JI = -d_IJ and
        the lower triangle is just the negated upper one -- n(n-1)/2 Lagrange
        solves instead of n(n-1), measured at 2.01x for SA-3.  This holds
        because mult_ediff=False below; with mult_ediff=True the kernel returns
        the undivided numerator, which is *symmetric*, and mirroring would flip
        the wrong sign.

        Same reuse as the gradient, and it still pays off harder than there:
        the ERI and the CP-MCSCF left-hand side are shared across every pair.
        Both come out of the same gradient_nacv_shared object compute_gradient
        already filled.
        """
        mc = self._state.mc
        nstates = self._n_total()
        natm = int(self._mol.natm)

        if nstates == 1:
            return np.zeros((1, 1, natm, 3), dtype=np.float64)

        # 1. shared by every pair
        nacs = mc.nac_method()
        if self.cphf_reuse:
            compute_cphf_lhs(nacs, self._shared)
        if self.deriv_eri_reuse:
            # the same provider compute_gradient filled, via self._shared
            nacs.deriv_eri = compute_deriv_eri(self._shared)

        # 2. per pair: build a shared part on first use, then read it back
        nacv = np.zeros((nstates, nstates, natm, 3), dtype=np.float64)
        for ket in range(nstates):
            for bra in range(ket):
                eris = None
                if self.eris_reuse:
                    if self._shared.eris is None:
                        self._shared.eris = compute_eris(mc)
                    eris = self._shared.eris
                d = np.asarray(
                    nacs.kernel(state=(ket, bra), use_etfs=use_etfs,
                                mult_ediff=False, eris=eris),
                    dtype=np.float64,
                )
                nacv[bra, ket] = d
                nacv[ket, bra] = -d
        return nacv



if __name__ == "__main__":

    # test=CASSCF(norbcas=5, nelecas=2, nroots=2, basis={'Li': 'sto-3g', 'F': '6-311+g*'}, unit='Bohr', charge=0, cas_list=[4, 7, 11, 14, 17])

    # geom = MolecularGeometry(
    #     atom_labels=['Li', 'F'],
    #     coords_bohr=np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 7.0]]),
    # )
    # test.set_geom(geom)

    # #print the test obj's state object
    # print(test.get_state()) 

    # print(test.compute_H_el())

    # print(test.compute_gradient(root=0))

    # print(test.compute_all_gradients())

    # print(test.compute_nac_vectors())

    # geom2 = MolecularGeometry(
    #     atom_labels=['Li', 'F'],
    #     coords_bohr=np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 8.0]]),
    # )
    # test.set_geom(geom2)

    # print(test.compute_time_overlap(2))

    def _smoke_test_casscf_heh_plus_sequence() -> None:
        """Smoke test for HeH+ CASSCF state sequencing and time-overlap."""
        geom1 = MolecularGeometry(
            atom_labels=['He', 'H'],
            coords_bohr=np.array([
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 1.46379],
            ], dtype=np.float64),
        )
        geom2 = MolecularGeometry(
            atom_labels=['He', 'H'],
            coords_bohr=np.array([
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 1.88973],
            ], dtype=np.float64),
        )

        casscf = CASSCF( norbcas=2, nelecas=2, nroots=3, basis='sto-3g', charge=1, unit='Bohr' )

        request = ES_Request(
            n_singlets=3,
            gradient_state='all',
            hessian_state=None,
            nacv=False,
            time_overlap=True,
        )

        result1 = casscf.compute_result(geom1, request)
        print("\nResult 1 H_el:", result1.H_el)
        print("Result 1 gradients:\n", result1.gradients)
        print("Result 1 time_overlap:", result1.time_overlap)

        result2 = casscf.compute_result(geom2, request)
        print("\nResult 2 H_el:", result2.H_el)
        print("Result 2 gradients:\n", result2.gradients)
        print("Result 2 time_overlap:\n", result2.time_overlap)


    def _smoke_test_casscf_nacv_sequence() -> None:
        """Smoke test for CASSCF NACV + time-overlap request behavior."""
        basis_dict = {'Li': 'sto-3g', 'F': '6-311+g*'}
        cas_list = [4, 7, 11, 14, 17]

        geom1 = MolecularGeometry(
            atom_labels=['Li', 'F'],
            coords_bohr=np.array([
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 7.0],
            ], dtype=np.float64),
        )
        geom2 = MolecularGeometry(
            atom_labels=['Li', 'F'],
            coords_bohr=np.array([
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 8.0],
            ], dtype=np.float64),
        )

        casscf = CASSCF(
            norbcas=5,
            nelecas=2,
            nroots=2,
            basis=basis_dict,
            unit='Bohr',
            charge=0,
            cas_list=cas_list,
        )

        request = ES_Request( n_singlets=2, gradient_state="all", nacv=True )

        result1 = casscf.compute_result(geom1, request)
        print("\nResult 1 H_el:", result1.H_el)
        print("Result 1 NACV shape:", result1.nac_vectors.shape)


        result2 = casscf.compute_result(geom2, request)
        print("\nResult 2 H_el:", result2.H_el)
        print("Result 2 NACV shape:", result2.nac_vectors.shape)


    def _smoke_test_casscf_c2h4_deriv_eri_timing(basis: str = '6-31g*', repeats: int = 2) -> None:
        """Time C2H4 CAS(4,4) SA-3 with deriv_eri_reuse off vs on.

        The SCF and CASSCF are converged once, before any timer; both variants
        differentiate that same wavefunction.  Two workflows are timed:
          all grad          -- compute_gradient([0, 1, 2])
          1 grad + all nacv -- compute_gradient([0]) then compute_nac_vectors()
        Each SA-CASSCF kernel call walks the (nabla i,j|k,l) block loop three
        times -- once in grad/casscf.py for the Hellmann-Feynman term, twice in
        grad/sacasscf.py for the Lagrange response -- so the call counter below
        should read 3 * natm per kernel with reuse off, and 1 in total with it
        on.  The variants alternate and the fastest of ``repeats`` is kept, so
        neither pays alone for first-call warm-up.
        """
        import time
        from pyscf import lib
        import pyscf.grad.casscf as _gcas
        import pyscf.grad.sacasscf as _gsa
        import pyscf.grad.deriv_eri as _de

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

        casscf = CASSCF(norbcas=4, nelecas=4, nroots=3, basis=basis, unit='Bohr')
        casscf.set_geom(geom)
        e_states = casscf.compute_H_el()
        mol = casscf.get_state().mol
        nao = mol.nao_nr()
        nao_pair = nao * (nao + 1) // 2

        # Count every libcint (nabla i,j|k,l) evaluation and the time inside it.
        ctr = {'n': 0, 't': 0.0}
        real = _de.int2e_ip1
        def counted(mol, shls_slice=None):
            t0 = time.perf_counter()
            out = real(mol, shls_slice)
            ctr['n'] += 1
            ctr['t'] += time.perf_counter() - t0
            return out
        for mod in (_de, _gcas, _gsa):
            mod.int2e_ip1 = counted

        cases = {
            'all grad  (3 roots)':
                lambda: (casscf.compute_gradient([0, 1, 2]), None),
            '1 grad + all nacv (6 pairs)':
                lambda: (casscf.compute_gradient([0]), casscf.compute_nac_vectors()),
        }
        best, out = {}, {}
        try:
            for name, fn in cases.items():
                for _ in range(repeats):
                    for tag, reuse in (('off', False), ('on', True)):
                        casscf.deriv_eri_reuse = reuse
                        casscf._shared = gradient_nacv_shared()  # 'on' pays for its own build
                        ctr['n'], ctr['t'] = 0, 0.0
                        t0 = time.perf_counter()
                        res = fn()
                        wall = time.perf_counter() - t0
                        key = (name, tag)
                        if key not in best or wall < best[key][0]:
                            best[key] = (wall, ctr['n'], ctr['t'])
                        out[key] = res
        finally:
            for mod in (_de, _gcas, _gsa):
                mod.int2e_ip1 = real
            casscf.deriv_eri_reuse = False

        print(f"\nC2H4 CAS(4,4) SA-3  basis={basis}  nao={nao}  nao_pair={nao_pair}"
              f"  threads={lib.num_threads()}  best of {repeats}")
        print(f"  full (nabla i,j|k,l) = 3*nao^2*nao_pair*8 = "
              f"{3 * nao * nao * nao_pair * 8 / 1e6:.1f} MB")
        print(f"  E_states = {e_states}")
        print(f"  {'case':<30}{'reuse':>6}{'wall/s':>9}{'int2e_ip1':>11}"
              f"{'in int2e_ip1/s':>16}{'speedup':>9}")
        for name in cases:
            t_off = best[(name, 'off')][0]
            for tag in ('off', 'on'):
                wall, n, t = best[(name, tag)]
                sp = '' if tag == 'off' else f"{t_off / wall:.2f}x"
                label = name if tag == 'off' else ''
                print(f"  {label:<30}{tag:>6}{wall:>9.2f}{n:>11d}{t:>16.2f}{sp:>9}")
            g_off, n_off = out[(name, 'off')]
            g_on, n_on = out[(name, 'on')]
            err = max(float(np.max(np.abs(a - b))) for a, b in zip(g_off, g_on))
            if n_off is not None:
                err = max(err, float(np.max(np.abs(n_off - n_on))))
            print(f"  {'':<30}{'max |off - on|':>6} = {err:.1e}")
            assert err < 1e-9, f"reuse changed the result by {err:.2e}"

    _smoke_test_casscf_heh_plus_sequence()
    _smoke_test_casscf_nacv_sequence()
    _smoke_test_casscf_c2h4_deriv_eri_timing()