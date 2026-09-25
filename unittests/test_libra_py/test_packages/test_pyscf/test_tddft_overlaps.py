"""Tests for determinant overlaps of restricted TDDFT/TDA pseudo-states."""
import unittest
from types import SimpleNamespace

import numpy as np
import pytest

pytest.importorskip("pyscf")
from pyscf import gto, lib

from libra_py.packages.pyscf.implementations.tddft import TDDFT, TDDFT_States
from libra_py.packages.pyscf.interfaces import ES_Request, MolecularGeometry


def restricted_state(mo_coeff, amplitudes, nocc):
    """Construct the minimal cached state needed by the overlap routines."""
    nmo = np.asarray(mo_coeff).shape[1]
    mo_occ = np.zeros(nmo)
    mo_occ[:nocc] = 2.0
    return TDDFT_States(
        mf=SimpleNamespace(mo_occ=mo_occ),
        mo_coeff=np.asarray(mo_coeff),
        amplitudes=[np.asarray(x) for x in amplitudes],
    )


class RestrictedTDDFTOverlapTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.old_threads = lib.num_threads()
        lib.num_threads(1)

    @classmethod
    def tearDownClass(cls):
        lib.num_threads(cls.old_threads)

    def setUp(self):
        self.strategy = TDDFT(atom_labels=("H", "H"), nexc=1, phase_tol=2.0)

    def test_two_electron_rotation_matches_normalized_singlet_formula(self):
        theta = 0.1
        c, s = np.cos(theta), np.sin(theta)
        previous = restricted_state(np.eye(2), [np.ones((1, 1))], nocc=1)
        current = restricted_state(np.array([[c, s], [-s, c]]),
                                   [np.ones((1, 1))], nocc=1)

        actual = self.strategy._cis_overlap(previous, current, np.eye(2), 2)
        expected = np.array([
            [c * c, np.sqrt(2.0) * c * s],
            [-np.sqrt(2.0) * c * s, c * c - s * s],
        ])
        np.testing.assert_allclose(actual, expected, atol=1e-13)

        # These entries specifically distinguish the determinant result from
        # the former det(Soo), S_ov, and Soo*Svv contraction.
        self.assertNotAlmostEqual(actual[0, 0], c)
        self.assertNotAlmostEqual(actual[0, 1], s)

    def test_identical_geometry_gives_orthonormal_singlet_states(self):
        amplitudes = [np.array([[1.0, 0.0]]), np.array([[0.0, 1.0]])]
        state = restricted_state(np.eye(3), amplitudes, nocc=1)
        actual = self.strategy._cis_overlap(state, state, np.eye(3), 3)
        np.testing.assert_allclose(actual, np.eye(3), atol=1e-14)

    def test_multielectron_reference_is_squared_determinant(self):
        # A nonorthogonal occupied block makes det(Soo)^2 visibly different
        # from the single determinant used by the old restricted contraction.
        s_mo = np.array([
            [0.94, 0.03, 0.02],
            [-0.01, 0.91, 0.04],
            [0.05, -0.02, 0.89],
        ])
        state = restricted_state(np.eye(3), [np.array([[1.0], [0.0]])], nocc=2)
        actual = self.strategy._cis_overlap(state, state, s_mo, 2)
        det_occ = np.linalg.det(s_mo[:2, :2])
        np.testing.assert_allclose(actual[0, 0], det_occ ** 2, atol=1e-14)

    def test_bad_amplitude_shape_fails_explicitly(self):
        state = restricted_state(np.eye(3), [np.ones((1, 1))], nocc=1)
        with self.assertRaisesRegex(ValueError, "amplitude shape"):
            self.strategy._cis_overlap(state, state, np.eye(3), 2)

    def test_unrestricted_determinant_path_remains_compatible(self):
        theta = 0.1
        c, s = np.cos(theta), np.sin(theta)
        rotation = np.array([[c, s], [-s, c]])
        mf = SimpleNamespace(
            mo_occ=(np.array([1.0, 0.0]), np.array([1.0, 0.0]))
        )
        amplitude = (
            np.ones((1, 1)) / np.sqrt(2.0),
            np.ones((1, 1)) / np.sqrt(2.0),
        )
        previous = TDDFT_States(
            mf=mf, mo_coeff=np.array([np.eye(2), np.eye(2)]),
            amplitudes=[amplitude],
        )
        current = TDDFT_States(
            mf=mf, mo_coeff=np.array([rotation, rotation]),
            amplitudes=[amplitude],
        )
        actual = self.strategy._cis_overlap(previous, current, np.eye(2), 2)
        expected = np.array([
            [c * c, np.sqrt(2.0) * c * s],
            [-np.sqrt(2.0) * c * s, c * c - s * s],
        ])
        np.testing.assert_allclose(actual, expected, atol=1e-13)

    def test_real_tda_request_uses_determinant_overlap(self):
        strategy = TDDFT(
            atom_labels=("H", "H"), nexc=1, basis="sto-3g", xc="lda,vwn",
            grid_level=0, use_tda=True,
        )
        request = ES_Request(n_singlets=2, time_overlap=True)
        first = strategy.compute_result(self.geometry(1.4), request)
        previous = strategy.copy().get_state()
        second = strategy.compute_result(self.geometry(1.45), request)
        current = strategy.get_state()

        self.assertIsNone(first.time_overlap)
        self.assertTrue(strategy.get_state().mf.converged)
        self.assertTrue(np.all(np.atleast_1d(strategy.get_state().td.converged)))
        direct = strategy._cis_overlap(
            previous,
            current,
            gto.intor_cross("int1e_ovlp", previous.mol, current.mol),
            2,
        )
        expected = strategy._align_time_overlap_phases(direct, strategy._phase_tol)
        np.testing.assert_allclose(second.time_overlap, expected, atol=1e-12)
        self_overlap = strategy._cis_overlap(
            current, current, current.mf.get_ovlp(), 2
        )
        np.testing.assert_allclose(self_overlap, np.eye(2), atol=1e-8)

    @staticmethod
    def geometry(distance):
        return MolecularGeometry(
            ("H", "H"), np.array([[0.0, 0.0, 0.0], [0.0, 0.0, distance]])
        )


if __name__ == "__main__":
    unittest.main()
