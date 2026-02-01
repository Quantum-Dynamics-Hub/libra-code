import numpy as np
import pytest

from libra_py.citools.interfaces import sd_and_csf_overlaps_singlet


# ============================================================
# Example runner (can be used outside pytest)
# ============================================================

def run_examples():
    results = []

    # ================= Case 1 ======================
    print("Case 1")
    s = np.array([
        [1.0, 0.2],
        [0.2, 1.2],
    ])
    S = np.kron(np.eye(2), s)

    params = dict(
        S=S,
        lowest_orbital=1,
        highest_orbital=2,
        nelec=2,
        homo_indx=1,
        common_sd_basis=[[1, 2]],
        active=None,
    )

    st_csf, st_sd = sd_and_csf_overlaps_singlet(
        params["S"],
        params["lowest_orbital"],
        params["highest_orbital"],
        params["nelec"],
        params["homo_indx"],
        params["common_sd_basis"],
        params["active"],
    )

    print("CSF overlap:\n", st_csf)
    print("SD overlap:\n", st_sd)
    results.append((st_csf, st_sd))

    # ================= Case 2 ======================
    print("\nCase 2")
    s = np.array([
        [-1.0, 0.0, 0.0],
        [0.0,  1.0, 0.2],
        [0.0,  0.2, 1.2],
    ])
    S = np.kron(np.eye(2), s)

    params = dict(
        S=S,
        lowest_orbital=1,
        highest_orbital=3,
        nelec=4,
        homo_indx=2,
        common_sd_basis=[[2, 3]],
        active=None,
    )

    st_csf, st_sd = sd_and_csf_overlaps_singlet(
        params["S"],
        params["lowest_orbital"],
        params["highest_orbital"],
        params["nelec"],
        params["homo_indx"],
        params["common_sd_basis"],
        params["active"],
    )

    print("CSF overlap:\n", st_csf)
    print("SD overlap:\n", st_sd)
    results.append((st_csf, st_sd))

    # ================= Case 3 ======================
    print("\nCase 3")
    s = np.array([
        [-1.0, 0.0, 0.0],
        [0.0,  1.0, 0.2],
        [0.0,  0.2, 1.2],
    ])
    S = np.kron(np.eye(2), s)

    params = dict(
        S=S,
        lowest_orbital=1,
        highest_orbital=3,
        nelec=2,
        homo_indx=2,
        common_sd_basis=[[2, 3]],
        active=[2, 3],
    )

    st_csf, st_sd = sd_and_csf_overlaps_singlet(
        params["S"],
        params["lowest_orbital"],
        params["highest_orbital"],
        params["nelec"],
        params["homo_indx"],
        params["common_sd_basis"],
        params["active"],
    )

    print("CSF overlap:\n", st_csf)
    print("SD overlap:\n", st_sd)
    results.append((st_csf, st_sd))

    return results


# ============================================================
# Pytest tests
# ============================================================

def test_all_cases_produce_identical_overlaps():
    """
    All three examples are physically equivalent constructions
    and must produce identical SD and CSF overlap matrices.
    """
    results = run_examples()

    ref_csf, ref_sd = results[0]

    for i, (st_csf, st_sd) in enumerate(results[1:], start=2):
        assert np.allclose(
            st_sd, ref_sd, atol=1e-12
        ), f"SD overlap differs in case {i}"

        assert np.allclose(
            st_csf, ref_csf, atol=1e-12
        ), f"CSF overlap differs in case {i}"


def test_overlap_matrices_are_symmetric_and_square():
    results = run_examples()

    for st_csf, st_sd in results:
        assert st_sd.ndim == 2
        assert st_csf.ndim == 2

        assert st_sd.shape[0] == st_sd.shape[1]
        assert st_csf.shape[0] == st_csf.shape[1]

        assert np.allclose(st_sd, st_sd.T, atol=1e-12)
        assert np.allclose(st_csf, st_csf.T, atol=1e-12)


def test_example_runner_produces_output(capsys):
    run_examples()
    out = capsys.readouterr().out
    assert "Case 1" in out
    assert "Case 2" in out
    assert "Case 3" in out

