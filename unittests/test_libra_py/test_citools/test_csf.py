import numpy as np

from libra_py.citools import csf


def _assert_csf_equal(actual, expected):
    actual_dict = {det: coeff for det, coeff in actual}
    expected_dict = {det: coeff for det, coeff in expected}

    assert actual_dict.keys() == expected_dict.keys()
    for det, coeff in expected_dict.items():
        assert np.isclose(actual_dict[det], coeff, atol=1e-12)


def test_two_unpaired_electrons_generate_expected_singlet_and_triplet_csfs():
    dets_with_parity = [
        ((-1, 2), 1),
        ((1, -2), 1),
    ]

    csfs = csf.generate_CSFs_grouped(dets_with_parity)

    assert len(csfs[(0.0, 0.0)]) == 1
    assert len(csfs[(1.0, 0.0)]) == 1

    inv_sqrt2 = 1.0 / np.sqrt(2.0)
    _assert_csf_equal(
        csfs[(0.0, 0.0)][0],
        [
            ((-1, 2), inv_sqrt2),
            ((1, -2), -inv_sqrt2),
        ],
    )
    _assert_csf_equal(
        csfs[(1.0, 0.0)][0],
        [
            ((-1, 2), inv_sqrt2),
            ((1, -2), inv_sqrt2),
        ],
    )


def test_closed_shell_determinant_is_a_pure_singlet_csf():
    csfs = csf.generate_CSFs_grouped([((1, -1, 2, -2), 1)])

    assert list(csfs.keys()) == [(0.0, 0.0)]
    assert csfs[(0.0, 0.0)] == [[((1, -1, 2, -2), 1.0)]]


def test_csf_groups_are_separated_by_spatial_occupation_pattern():
    dets_with_parity = [
        ((1, -1, -2, 3), 1),
        ((1, -1, 2, -3), 1),
        ((-1, 2, -2, 3), 1),
        ((1, 2, -2, -3), 1),
    ]

    csfs = csf.generate_CSFs_grouped(dets_with_parity)

    assert len(csfs[(0.0, 0.0)]) == 2
    for singlet_csf in csfs[(0.0, 0.0)]:
        assert np.isclose(
            sum(coeff * coeff for _, coeff in singlet_csf),
            1.0,
            atol=1e-12,
        )

