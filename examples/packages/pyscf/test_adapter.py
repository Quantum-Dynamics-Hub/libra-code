# *********************************************************************************
# * Copyright (C) 2026 Jieyang Gu <jieyanggu792@gmail.com>
# * This file is distributed under the terms of the GNU General Public License
# * as published by the Free Software Foundation, either version 3 of
# * the License, or (at your option) any later version.
# * See the file LICENSE in the root directory of this distribution
# * or <http://www.gnu.org/licenses/>.
# *
# *********************************************************************************/
"""Unit test for the PySCF ES adapter using the CASSCF backend."""

from __future__ import annotations

import unittest

import numpy as np

from libra_py.packages.pyscf.adapter import LibraESAdapter
from libra_py.packages.pyscf.implementations.casscf import CASSCF
from libra_py.packages.pyscf.interfaces import ES_Request


class FakeNham:
    def __init__(self) -> None:
        self.values: dict[str, object] = {}

    def set_ham_adi_by_val(self, value):
        self.values["ham"] = value

    def set_ham_dia_by_val(self, value):
        self.values["soc"] = value

    def set_d1ham_adi_by_val(self, value):
        self.values["d1"] = value

    def set_dc1_adi_by_val(self, value):
        self.values["dc1"] = value

    def set_time_overlap_adi_by_val(self, value):
        self.values["time"] = value


class TestLibraESAdapter(unittest.TestCase):
    def test_adapter_produces_energies_gradients_nac_and_time_overlap(self) -> None:
        backend = CASSCF(
            norbcas=2,
            nelecas=2,
            nroots=3,
            basis="sto-3g",
            charge=1,
            unit="Bohr",
        )
        request = ES_Request(
            n_singlets=3,
            n_triplets=0,
            H_soc=False,
            gradient_state="all",
            hessian_state=None,
            time_overlap=True,
            nacv=True,
        )
        adapter = LibraESAdapter(
            backend=backend,
            atom_labels=["He", "H"],
            request_template=request,
        )

        q1 = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 1.46379], dtype=float)
        nham1 = FakeNham()
        result1 = adapter.compute_adiabatic_callback(q1, {"nham": nham1, "request": request})

        self.assertIsNotNone(result1.H_el)
        self.assertEqual(result1.H_el.shape, (3,))
        self.assertIsNotNone(result1.gradients)
        self.assertEqual(len(result1.gradients), 3)
        self.assertIsNotNone(result1.nac_vectors)
        self.assertEqual(result1.nac_vectors.shape, (3, 3, 2, 3))
        self.assertIsNotNone(result1.time_overlap)
        self.assertEqual(result1.time_overlap.shape, (3, 3))
        self.assertIn("ham", nham1.values)
        self.assertIn("d1", nham1.values)
        self.assertIn("dc1", nham1.values)
        self.assertIn("time", nham1.values)

        q2 = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 1.88973], dtype=float)
        previous_backend = backend.copy()
        nham2 = FakeNham()
        result2 = adapter.compute_adiabatic_callback(
            q2,
            {
                "nham": nham2,
                "request": request,
                "previous_strategy": previous_backend,
            },
        )

        self.assertIsNotNone(result2.H_el)
        self.assertIsNotNone(result2.time_overlap)
        self.assertTrue(np.all(np.isfinite(result2.time_overlap)))
        self.assertGreater(np.abs(result2.time_overlap[0, 0]), 0.9)
        self.assertGreater(np.abs(result2.time_overlap[1, 1]), 0.9)
        self.assertGreater(np.abs(result2.time_overlap[2, 2]), 0.9)


if __name__ == "__main__":
    unittest.main()
