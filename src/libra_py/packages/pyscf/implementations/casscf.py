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
import copy as pycopy
from dataclasses import dataclass
from typing import Any, List, Optional, Sequence, Tuple, Union
import numpy as np
from pyscf import fci, gto, mcscf, scf
from interface.interfaces import ES_Request, ES_Strategy, MolecularGeometry

BOHR_TO_ANG = 0.529177210903


@dataclass
class CASSCF_States:
    mol: Optional[Any] = None
    mf: Optional[Any] = None
    mc: Optional[Any] = None
    # MO-basis ERI at the converged CASSCF orbitals, mc.ao2mo(mc.mo_coeff).
    # Shared by compute_gradient and compute_nac_vectors; built lazily, so
    # energy-only steps never pay for it.
    eris: Optional[Any] = None 

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
        self._geom: Optional[MolecularGeometry] = None
        self._state: CASSCF_States | None = None
        self._previous_state: CASSCF_States | None = None

    @property
    def nroots(self) -> int:
        """Number of roots this strategy solves for.

        Fixed at construction, and the single source of truth for the state
        count: the adapter builds ``ES_Request.n_singlets`` from it rather than
        from a separate ``model_params["nstates"]`` entry.
        """
        return int(self._nroots)

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
            self._state.mf = scf.RHF(self._mol).run(verbose=0)
            return

        try:
            prev_dm = prev_state.mf.make_rdm1()
        except Exception:
            prev_dm = None

        if prev_dm is not None:
            try:
                self._state.mf = scf.RHF(self._mol).run(dm0=prev_dm, verbose=0)
                return
            except Exception:
                pass

        try:
            mf = scf.RHF(self._mol)
            mf.init_guess_by_mo(prev_state.mf.mo_coeff)
            self._state.mf = mf.run(verbose=0)
        except Exception:
            self._state.mf = scf.RHF(self._mol).run(verbose=0)

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
            spin=0,
        )

        self._run_scf()

        # Write the current SCF result back into the state object.
        if self._state is None:
            self._state = CASSCF_States(mol=self._mol, mf=None, mc=None)
        self._state.mol = self._mol
        self._state.mc = None
        self._state.eris = None   # new geometry -> the cached integrals are stale

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
        return pycopy.deepcopy(self)

    def compute_H_el(self) -> np.ndarray:
        if self._previous_state is not None and self._previous_state.mc is not None:
            mocoeff = self._previous_state.mc.mo_coeff
        else:
            mocoeff = self._state.mf.mo_coeff

        mc = mcscf.CASSCF(self._state.mf, self._norbcas, self._nelecas)
        mc.fcisolver = fci.direct_spin0.FCI(self._state.mol)
        mc.fcisolver.nroots = self._nroots
        if self._nroots > 1:
            mc = mc.state_average_([1.0 / self._nroots] * self._nroots)  # equal weights as default; required for gradients in pyscf

        if self._cas_list is not None:
            mocoeff = mcscf.sort_mo(mc, mocoeff, self._cas_list)

        mc.kernel(mocoeff)

        self._state.mc = mc
        self._state.eris = None   # new wavefunction -> the cached integrals are stale

        #write the energies to the H_el attribute of the state object and return the energies as a numpy array
        self._state.H_el = np.asarray(mc.e_states, dtype=np.float64)
        return self._state.H_el

    #gradient 

    def _share_eris(self, solver: Any) -> Any:
        """transform the integrals to the CASSCF MO basis and save in the state object"""
        st = self._state
        if st.eris is None:
            st.eris = st.mc.ao2mo(st.mc.mo_coeff)
        for name in ("make_fcasscf", "make_fcasscf_nacs"):
            builder = getattr(solver, name, None)
            if builder is None:
                continue
            def patched(*args, _build=builder, **kwargs):
                fcasscf = _build(*args, **kwargs)
                fcasscf.ao2mo = lambda mo_coeff=None: st.eris
                return fcasscf
            setattr(solver, name, patched)
        return solver

    def compute_gradient(self, roots: Sequence[int]) -> list[np.ndarray]:
        """Nuclear gradients for ``roots``, returned in the order requested."""
        mc = self._state.mc
        roots = list(roots)
        n_avail = len(mc.e_states)
        for root in roots:
            if not 0 <= root < n_avail:
                raise IndexError(
                    f"Requested root {root}, but only {n_avail} roots are available."
                )

        # 1. compute the intermediate shared by every root
        grad = self._share_eris(mc.nuc_grad_method(state=0))

        # 2. per root
        return [
            np.asarray(grad.kernel(state=root, eris=self._state.eris), dtype=np.float64)
            for root in roots
        ]


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
        prev_casci.fcisolver = fci.direct_spin0.FCISolver(prev_state.mol)
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
        curr_casci.fcisolver = fci.direct_spin0.FCISolver(curr_state.mol)
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

        prev_roots = [np.asarray(vec) for vec in prev_roots_raw[:nroots]]
        curr_roots = [np.asarray(vec) for vec in curr_roots_raw[:nroots]]
        prev_act = np.asarray(prev_casci.mo_coeff)[:, : prev_casci.ncas]
        curr_act = np.asarray(curr_casci.mo_coeff)[:, : curr_casci.ncas]
        s12_mo = prev_act.T.conj() @ ao_overlap @ curr_act
        overlap = np.zeros((nroots, nroots), dtype=float)
        for i in range(nroots):
            for j in range(nroots):
                overlap[i, j] = fci.addons.overlap(
                    prev_roots[i],
                    curr_roots[j],
                    self._norbcas,
                    self._nelecas if isinstance(self._nelecas, int) else sum(self._nelecas),
                    s=s12_mo,
                )

        overlap = np.asarray(np.real_if_close(overlap))
        for j in range(nroots):
            if overlap[j, j] < 0:
                overlap[:, j] = -overlap[:, j]
        return overlap

    def compute_nac_vectors(self, use_etfs: bool = True) -> np.ndarray:
        """Non-adiabatic coupling vectors for every ordered state pair.

        Same split as the gradient, and it pays off harder: this is one Lagrange
        solve per *ordered pair*, so an n-state run reuses the shared ERI and
        Jacobian across n(n-1) solves rather than n.  Both come out of the same
        CASSCF_States cache that compute_gradient already populated.
        """
        mc = self._state.mc
        nstates = self._n_total()
        natm = int(self._mol.natm)

        # 1. shared by every pair
        nacs = self._share_eris(mc.nac_method())

        # 2. per pair
        nacv = np.zeros((nstates, nstates, natm, 3), dtype=np.float64)
        for ket in range(nstates):
            for bra in range(nstates):
                if bra == ket:
                    continue
                nacv[bra, ket] = np.asarray(
                    nacs.kernel(state=(ket, bra), use_etfs=use_etfs,
                                mult_ediff=False, eris=self._state.eris),
                    dtype=np.float64,
                )
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

    _smoke_test_casscf_heh_plus_sequence()
    _smoke_test_casscf_nacv_sequence()
