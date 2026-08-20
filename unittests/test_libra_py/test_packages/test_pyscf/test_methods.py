import copy

import numpy as np

from liblibra_core import MATRIX
from libra_py.packages.pyscf.interfaces import ES_Result
from libra_py.packages.pyscf.methods import pyscf_compute_adi


class FakeStrategy:
    def __init__(self, multiplicity, energies, overlap):
        self.spin_multiplicity = multiplicity
        self.nroots = len(energies)
        self.energies = np.asarray(energies)
        self.overlap = np.asarray(overlap)

    def copy(self):
        return copy.deepcopy(self)

    def compute_result(self, geometry, request):
        gradients = [
            np.full((len(geometry.atom_labels), 3), root + 1.0)
            for root in range(self.nroots)
        ]
        nac = np.zeros((self.nroots, self.nroots, 1, 3))
        if self.nroots > 1:
            nac[0, 1, 0, 0] = 0.25
            nac[1, 0, 0, 0] = -0.25
        return ES_Result(
            H_el=self.energies,
            gradients=gradients,
            nac_vectors=nac,
            time_overlap=self.overlap,
        )


def test_mixed_spin_manifolds_expand_to_openmolcas_layout():
    q = MATRIX(3, 1)
    params = {
        "atom_labels": ["X"],
        "dt": 2.0,
        "energy_zero": -10.0,
        "spin_manifolds": [
            {
                "spin": 2,
                "nroots": 2,
                "es_strategy": FakeStrategy(2, [-9.0, -8.0], [[1.0, 0.2], [-0.1, 1.0]]),
            },
            {
                "spin": 4,
                "nroots": 1,
                "es_strategy": FakeStrategy(4, [-7.0], [[1.0]]),
            },
        ],
    }

    obj = pyscf_compute_adi(q, params, [0])

    assert obj.spin_labels == [
        (2, 0, 1), (2, 0, -1),
        (2, 1, 1), (2, 1, -1),
        (4, 0, 3), (4, 0, 1), (4, 0, -1), (4, 0, -3),
    ]
    assert [obj.ham_adi.get(i, i).real for i in range(8)] == [
        1.0, 1.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0,
    ]

    # Root couplings are copied only between equal-Ms components.
    assert obj.time_overlap_adi.get(0, 2).real == 0.2
    assert obj.time_overlap_adi.get(1, 3).real == 0.2
    assert obj.time_overlap_adi.get(0, 3) == 0.0j
    assert obj.time_overlap_adi.get(0, 4) == 0.0j
    assert obj.dc1_adi[0].get(0, 2).real == 0.25
    assert obj.dc1_adi[0].get(0, 3) == 0.0j

    # Gradients are repeated for every Ms component of a spin-free root.
    assert [obj.d1ham_adi[0].get(i, i).real for i in range(8)] == [
        1.0, 1.0, 2.0, 2.0, 1.0, 1.0, 1.0, 1.0,
    ]
