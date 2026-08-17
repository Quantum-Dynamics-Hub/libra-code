import numpy as np

from libra_py.citools.interfaces import sd_and_csf_overlaps_singlet
from libra_py.citools import interfaces
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


def _run_equivalent_active_space_examples():
    results = []

    s = np.array([
        [1.0, 0.2],
        [0.2, 1.2],
    ])
    results.append(sd_and_csf_overlaps_singlet(
        np.kron(np.eye(2), s),
        lowest_orbital=1,
        highest_orbital=2,
        nelec=2,
        homo_indx=1,
        common_sd_basis=[[1, 2]],
        _active_space=None,
    ))

    s = np.array([
        [-1.0, 0.0, 0.0],
        [0.0, 1.0, 0.2],
        [0.0, 0.2, 1.2],
    ])
    results.append(sd_and_csf_overlaps_singlet(
        np.kron(np.eye(2), s),
        lowest_orbital=1,
        highest_orbital=3,
        nelec=4,
        homo_indx=2,
        common_sd_basis=[[2, 3]],
        _active_space=None,
    ))
    results.append(sd_and_csf_overlaps_singlet(
        np.kron(np.eye(2), s),
        lowest_orbital=1,
        highest_orbital=3,
        nelec=2,
        homo_indx=2,
        common_sd_basis=[[2, 3]],
        _active_space=[2, 3],
    ))

    return results


def test_equivalent_active_space_constructions_produce_same_overlaps():
    results = _run_equivalent_active_space_examples()
    ref_csf, ref_sd = results[0]

    for st_csf, st_sd in results[1:]:
        assert np.allclose(st_sd, ref_sd, atol=1e-12)
        assert np.allclose(st_csf, ref_csf, atol=1e-12)


def test_same_time_overlap_matrices_are_symmetric_and_square():
    results = _run_equivalent_active_space_examples()

    for st_csf, st_sd in results:
        assert st_sd.ndim == 2
        assert st_csf.ndim == 2
        assert st_sd.shape[0] == st_sd.shape[1]
        assert st_csf.shape[0] == st_csf.shape[1]
        assert np.allclose(st_sd, st_sd.T, atol=1e-12)
        assert np.allclose(st_csf, st_csf.T, atol=1e-12)


def test_configs_and_T_matrix_singlet_is_orthonormal_for_multiple_excitations():
    nelec = 4
    homo_indx = 2
    active_space = [1, 2, 3]
    common_sd_basis = [[2, 3], [1, 3]]

    gs = sd.make_ref_det(nelec, homo_indx)
    configs0_raw = [tuple(gs)]
    for occ, vir in common_sd_basis:
        configs0_raw.append(tuple(sd.make_excitation(gs, occ, vir)))

    mapped_basis, T = interfaces.configs_and_T_matrix_singlet(
        configs0_raw,
        active_space=active_space,
        orbital_space=active_space,
        nelec=nelec,
        S=0,
        Ms=0,
    )

    T_dense = T.toarray()

    assert len(mapped_basis) == 2 * len(common_sd_basis) + 1
    assert len(set(mapped_basis)) == len(mapped_basis)
    assert np.allclose(
        T_dense.conj().T @ T_dense,
        np.eye(T_dense.shape[1]),
        atol=1e-12,
    )


def test_factorized_sd_overlap_matches_full_spin_orbital_determinant():
    alpha_overlap = np.array([
        [1.00, 0.11, 0.02],
        [0.13, 0.91, 0.03],
        [0.04, 0.05, 1.08],
    ])
    beta_overlap = np.array([
        [0.97, 0.07, 0.06],
        [0.09, 1.03, 0.04],
        [0.08, 0.02, 0.94],
    ])
    spin_orbital_overlap = np.block([
        [alpha_overlap, np.zeros_like(alpha_overlap)],
        [np.zeros_like(beta_overlap), beta_overlap],
    ])

    nelec = 4
    homo_indx = 2
    common_sd_basis = [[2, 3], [1, 3]]

    st_csf, st_sd = sd_and_csf_overlaps_singlet(
        spin_orbital_overlap,
        lowest_orbital=1,
        highest_orbital=3,
        nelec=nelec,
        homo_indx=homo_indx,
        common_sd_basis=common_sd_basis,
        _active_space=None,
    )

    gs = sd.make_ref_det(nelec, homo_indx)
    configs0_raw = [tuple(gs)]
    for occ, vir in common_sd_basis:
        configs0_raw.append(tuple(sd.make_excitation(gs, occ, vir)))

    mapped_basis, T = interfaces.configs_and_T_matrix_singlet(
        configs0_raw,
        active_space=[1, 2, 3],
        orbital_space=[1, 2, 3],
        nelec=nelec,
        S=0,
        Ms=0,
    )

    ref_sd = _full_spin_orbital_overlap(mapped_basis, spin_orbital_overlap)
    ref_csf = T.T @ ref_sd @ T

    assert np.allclose(st_sd, ref_sd, atol=1e-12)
    assert np.allclose(st_csf, ref_csf, atol=1e-12)


def test_triplet_overlap_keeps_singlet_ground_state_and_triplet_excitation():
    spin_orbital_overlap = np.eye(4)

    st_csf, st_sd = interfaces.sd_and_csf_overlaps(
        spin_orbital_overlap,
        lowest_orbital=1,
        highest_orbital=2,
        nelec=2,
        homo_indx=1,
        common_sd_basis=[[1, 2]],
        S=1,
        Ms=0,
    )

    assert st_sd.shape == (3, 3)
    assert st_csf.shape == (2, 2)
    assert np.allclose(st_csf, np.eye(2), atol=1e-12)


def test_explicit_doublet_reference_supports_odd_electron_count():
    st_csf, st_sd = interfaces.sd_and_csf_overlaps(
        np.eye(4),
        lowest_orbital=1,
        highest_orbital=2,
        nelec=1,
        homo_indx=1,
        common_sd_basis=[[1, 2]],
        S=0.5,
        Ms=0.5,
        reference_det=[1],
    )

    assert st_sd.shape == (2, 2)
    assert np.allclose(st_csf, np.eye(2), atol=1e-12)


def test_explicit_triplet_reference_does_not_add_closed_shell_state():
    raw_configs = [(1, 2)]
    mapped_basis, _ = interfaces.configs_and_T_matrix(
        raw_configs, [1, 2], [1, 2], nelec=2, S=1, Ms=0,
    )
    st_csf, st_sd = interfaces.sd_and_csf_overlaps(
        np.eye(4),
        lowest_orbital=1,
        highest_orbital=2,
        nelec=2,
        homo_indx=2,
        common_sd_basis=[],
        S=1,
        Ms=0,
        reference_det=[1, 2],
    )

    assert all(not (orb in det and -orb in det) for det in mapped_basis for orb in det)
    assert st_sd.shape == (2, 2)
    assert st_csf.shape == (1, 1)
    assert np.allclose(st_csf, np.eye(1), atol=1e-12)


def test_shared_spin_helpers_cover_open_shell_package_interfaces():
    assert interfaces.spin_quantum_numbers(2, -0.5) == (0.5, -0.5)
    assert interfaces.spin_quantum_numbers(3, None) == (1.0, 1.0)
    assert interfaces.make_open_shell_reference(1, 4) == [1, -1, 2, 3, 4]
    assert interfaces.reference_from_electron_count(5, 4, 4) == [1, -1, 2, 3, 4]


def test_reference_builder_rejects_incompatible_electron_parity():
    with np.testing.assert_raises_regex(ValueError, "incompatible"):
        interfaces.reference_from_electron_count(4, 2)
