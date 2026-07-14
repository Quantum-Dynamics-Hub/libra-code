import numpy as np
import pytest

from libra_py.citools import slatdet as sd


def _spin_orbital_indices(det, nspatial):
    return [abs(o) - 1 if o > 0 else nspatial + abs(o) - 1 for o in det]


def _full_spin_orbital_overlap(dets, spin_orbital_overlap):
    res = np.zeros((len(dets), len(dets)))
    nspatial = spin_orbital_overlap.shape[0] // 2

    for i, bra in enumerate(dets):
        bra_indx = _spin_orbital_indices(bra, nspatial)
        for j, ket in enumerate(dets):
            ket_indx = _spin_orbital_indices(ket, nspatial)
            res[i, j] = np.linalg.det(
                spin_orbital_overlap[np.ix_(bra_indx, ket_indx)]
            )

    return res


@pytest.mark.parametrize(
    "det,expected",
    [
        ((1, -1), 1),
        ((1, -1, 2, -2), -1),
        ((1, -1, 2, -2, 3, -3), -1),
        ((1, -1, -2, 3), 1),
        ((1, -1, 2, -3), -1),
        ((-1, 2, -2, 3), -1),
        ((1, 2, -2, -3), 1),
    ],
)
def test_alpha_beta_ordering_phase_known_cases(det, expected):
    assert sd.alpha_beta_ordering_phase(det) == expected


def test_slater_overlap_uses_beta_block_of_doubled_spin_matrix():
    alpha_overlap = np.array([
        [1.0, 0.2],
        [0.3, 1.1],
    ])
    beta_overlap = np.array([
        [0.7, 0.4],
        [0.5, 1.3],
    ])
    spin_overlap = np.block([
        [alpha_overlap, np.zeros_like(alpha_overlap)],
        [np.zeros_like(beta_overlap), beta_overlap],
    ])

    dets = [(1, -1), (-1, 2), (1, -2)]
    phases = [sd.alpha_beta_ordering_phase(det) for det in dets]

    st_sd = sd.slater_overlap_matrix(
        dets,
        dets,
        spin_overlap,
        phases_A=phases,
        phases_B=phases,
        spin_orbital_matrix=True,
    )
    ref_sd = _full_spin_orbital_overlap(dets, spin_overlap)

    assert np.allclose(st_sd, ref_sd, atol=1e-12)
    assert not np.allclose(alpha_overlap, beta_overlap)


def test_slater_overlap_spatial_matrix_matches_doubled_equal_blocks():
    spatial_overlap = np.array([
        [1.0, 0.2],
        [0.3, 1.1],
    ])
    spin_overlap = np.kron(np.eye(2), spatial_overlap)
    dets = [(1, -1), (-1, 2), (1, -2)]
    phases = [sd.alpha_beta_ordering_phase(det) for det in dets]

    spatial_result = sd.slater_overlap_matrix(
        dets,
        dets,
        spatial_overlap,
        phases_A=phases,
        phases_B=phases,
    )
    doubled_result = sd.slater_overlap_matrix(
        dets,
        dets,
        spin_overlap,
        phases_A=phases,
        phases_B=phases,
        spin_orbital_matrix=True,
    )

    assert np.allclose(spatial_result, doubled_result, atol=1e-12)


def test_slater_overlap_rejects_bad_doubled_matrix_shape():
    with pytest.raises(ValueError, match="Doubled spin-orbital overlap"):
        sd.slater_overlap_matrix(
            [(1, -1)],
            [(1, -1)],
            np.eye(3),
            spin_orbital_matrix=True,
        )


def test_make_ref_det_uses_homo_aligned_closed_shell_window():
    assert sd.make_ref_det(nelec=4, homo_indx=5) == [4, -4, 5, -5]


def test_make_excitation_preserves_reference_order_before_canonicalization():
    ref = [1, -1, 2, -2]
    assert sd.make_excitation(ref, 1, 3) == [3, -1, 2, -2]

