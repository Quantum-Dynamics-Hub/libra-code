"""Small real-PySCF tests of selectable CASSCF time-overlap representations."""
import itertools
import unittest
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pytest

pytest.importorskip("pyscf")
from pyscf import fci, gto, lib, mcscf
from libra_py.packages.pyscf.implementations.casscf import CASSCF, CASSCF_States
from libra_py.packages.pyscf.interfaces import ES_Request, MolecularGeometry


def determinant_overlap(left, right, orbital_overlap, nelec):
    """Independent small-space determinant expansion (no fci.addons.overlap)."""
    occupations = [sorted(itertools.combinations(range(len(orbital_overlap)), n),
                          key=lambda occ: sum(1 << i for i in occ)) for n in nelec]
    value = 0.0
    for ia, a in enumerate(occupations[0]):
        for ib, b in enumerate(occupations[1]):
            for ic, c in enumerate(occupations[0]):
                for id_, d in enumerate(occupations[1]):
                    value += (left[ia, ib].conjugate() * right[ic, id_]
                              * np.linalg.det(orbital_overlap[np.ix_(a, c)])
                              * np.linalg.det(orbital_overlap[np.ix_(b, d)]))
    return value


def expected_overlap(previous, current, roots_previous, roots_current, mo_previous, mo_current):
    s = mo_previous.T @ gto.intor_cross("int1e_ovlp", previous.mol, current.mol) @ mo_current
    result = np.array([[determinant_overlap(a, b, s, previous.mc.nelecas)
                        for b in roots_current] for a in roots_previous])
    for j in range(len(result)):
        if result[j, j] < 0:
            result[:, j] *= -1
    return result


class CASSCFOverlapTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.old_threads = lib.num_threads()
        lib.num_threads(1)
        cls.strategy = CASSCF(norbcas=2, nelecas=2, nroots=2, overlap_orbitals="casscf")
        request = ES_Request(n_singlets=2, time_overlap=True)
        cls.first_result = cls.strategy.compute_result(cls.geometry(3.0), request)
        cls.previous = cls.strategy.copy().get_state()
        cls.second_result = cls.strategy.compute_result(cls.geometry(3.1), request)
        cls.current = cls.strategy.get_state()
        for state in (cls.previous, cls.current):
            assert state.mf.converged and state.mc.converged

    @classmethod
    def tearDownClass(cls):
        lib.num_threads(cls.old_threads)

    @staticmethod
    def geometry(distance):
        return MolecularGeometry(("Li", "H"), np.array([[0., 0., 0.], [0., 0., distance]]))

    def test_default_is_hf_and_invalid_selection_fails(self):
        self.assertEqual(CASSCF().overlap_orbitals, "hf")
        for invalid in ("dft", "CASSCF", None):
            with self.assertRaisesRegex(ValueError, "overlap_orbitals"):
                CASSCF(overlap_orbitals=invalid)

    def test_hf_default_matches_auxiliary_casci_determinants(self):
        roots, active = [], []
        for state in (self.previous, self.current):
            cas = mcscf.CASCI(state.mf, 2, 2)
            cas.fcisolver = fci.direct_spin0.FCI(state.mol)
            cas.fcisolver.nroots = 2
            cas.verbose = 0
            cas.kernel()
            roots.append(cas.ci)
            active.append(cas.mo_coeff[:, cas.ncore:cas.ncore+cas.ncas])
        expected = expected_overlap(self.previous, self.current, *roots, *active)
        default = CASSCF(norbcas=2, nelecas=2, nroots=2)
        explicit = CASSCF(norbcas=2, nelecas=2, nroots=2, overlap_orbitals="hf")
        actual = default.compute_time_overlap(self.current, self.previous)
        np.testing.assert_allclose(actual, explicit.compute_time_overlap(self.current, self.previous), atol=1e-12)
        # Independent CI diagonalizations can choose different root signs.
        np.testing.assert_allclose(abs(actual), abs(expected), atol=1e-9)

    def test_optimized_overlaps_match_determinants_without_resolving_ci(self):
        active = [s.mc.mo_coeff[:, s.mc.ncore:s.mc.ncore+s.mc.ncas]
                  for s in (self.previous, self.current)]
        expected = expected_overlap(self.previous, self.current, self.previous.mc.ci, self.current.mc.ci, *active)
        # Optimized mode must not create another CASCI calculation or reset roots.
        with patch("libra_py.packages.pyscf.implementations.casscf.mcscf.CASCI",
                   side_effect=AssertionError("Unexpected auxiliary CASCI")):
            actual = self.strategy.compute_time_overlap(self.current, self.previous)
        np.testing.assert_allclose(actual, expected, atol=1e-12)
        np.testing.assert_allclose(actual, self.second_result.time_overlap, atol=1e-12)
        self.assertIsNone(self.first_result.time_overlap)
        hf = CASSCF(norbcas=2, nelecas=2, nroots=2).compute_time_overlap(self.current, self.previous)
        self.assertGreater(np.linalg.norm(abs(actual)-abs(hf)), 1e-5)

    def test_initial_geometry_identity_and_independent_copy(self):
        # The first solve has orthonormal orbitals. The pre-existing sequential
        # warm start does not restore the MO metric at a displaced geometry;
        # displaced overlaps above are tested against their actual determinants.
        np.testing.assert_allclose(self.strategy.compute_time_overlap(self.previous, self.previous), np.eye(2), atol=1e-10)
        cloned = self.strategy.copy()
        self.assertEqual(cloned.overlap_orbitals, "casscf")
        cloned.get_state().mc.ci[0][0, 0] += 0.1
        self.assertFalse(np.array_equal(cloned.get_state().mc.ci[0], self.current.mc.ci[0]))

    def test_single_root_array_and_column_phase_alignment(self):
        # Use one normalized stored root, but represent it as a single-root array.
        def state(source, sign):
            mc = SimpleNamespace(ncas=source.mc.ncas, ncore=source.mc.ncore,
                                 nelecas=source.mc.nelecas, mo_coeff=source.mc.mo_coeff,
                                 ci=sign*source.mc.ci[0])
            return CASSCF_States(mol=source.mol, mc=mc)
        single = CASSCF(norbcas=2, nelecas=2, overlap_orbitals="casscf")
        actual = single.compute_time_overlap(state(self.current, -1), state(self.previous, 1))
        self.assertEqual(actual.shape, (1, 1))
        self.assertGreater(actual[0, 0], 0)
        np.testing.assert_allclose(actual[0, 0], self.second_result.time_overlap[0, 0], atol=1e-12)

    def test_missing_optimized_states_do_not_fall_back_to_hf(self):
        missing = CASSCF_States(mol=self.current.mol, mf=self.current.mf)
        with self.assertRaisesRegex(ValueError, "requires optimized CASSCF"):
            self.strategy.compute_time_overlap(missing, self.previous)
        bad = self.strategy.copy().get_state()
        bad.mc.ci = bad.mc.ci[:1]
        with self.assertRaisesRegex(ValueError, "requested determinant CI roots"):
            self.strategy.compute_time_overlap(bad, self.previous)

    def test_open_shell_uses_snapshot_electron_partition(self):
        # A synthetic highest-projection quartet: three alpha electrons, no beta.
        # This distinguishes mc.nelecas=(3,0) from an ambiguous integer count of 3.
        s = np.array([[.96, .02, .01], [.01, .94, .03], [0., .02, .97]])
        mc = SimpleNamespace(ncas=3, ncore=0, nelecas=(3, 0), mo_coeff=np.eye(3), ci=np.ones((1, 1)))
        state = CASSCF_States(mol=self.current.mol, mc=mc)
        strategy = CASSCF(norbcas=3, nelecas=3, spin_multiplicity=4, overlap_orbitals="casscf")
        with patch.object(strategy, "_compute_ao_overlap", return_value=s):
            actual = strategy.compute_time_overlap(state, state)
        np.testing.assert_allclose(actual, [[np.linalg.det(s)]], atol=1e-12)


if __name__ == "__main__":
    unittest.main()
