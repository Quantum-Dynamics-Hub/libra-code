import numpy as np
import pytest

from libra_py.citools.ci import overlap


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

    # lowest = homo - nocc
    # highest = homo + nvirt
    params = dict(
        S=S,
        nelec=2,
        nocc=0,
        nvirt=1,
        homo_indx=1,
        nstates = 2,
        data = [ [],
                 [ [[1, 2]] ],
                 [ [1.0 ] ]
               ],
        active_space=None,
    )

    st_ci = overlap(params["S"], params["data"], params["data"], params)

    print("CI overlap:\n", st_ci)
    results.append(st_ci)

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
        nelec=4,
        nocc=1,
        nvirt=1,
        homo_indx=2,
        nstates = 2,
        data = [ [],
                 [ [[2, 3]] ],
                 [ [1.0 ] ]
               ],
        active_space=None,
    )

    st_ci = overlap(params["S"], params["data"], params["data"], params)

    print("CI overlap:\n", st_ci)
    results.append(st_ci)

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
        nelec=2,
        nocc=1,
        nvirt=1,
        homo_indx=2,
        nstates = 2,
        data = [ [],
                 [ [[2, 3]] ],
                 [ [1.0 ] ]
               ],
        active_space=[2, 3],
    )

    st_ci = overlap(params["S"], params["data"], params["data"], params)

    print("CI overlap:\n", st_ci)
    results.append(st_ci)

    return results


# Uncomment to run the example
#run_examples()


# ============================================================
# Pytest tests
# ============================================================

def test_all_cases_produce_identical_overlaps():
    """
    All three examples are physically equivalent constructions
    and must produce identical CI overlap matrices.
    """
    results = run_examples()

    ref_ci = results[0]

    for i, st_ci in enumerate(results[1:], start=2):
        assert np.allclose(
            st_ci, ref_ci, atol=1e-12
        ), f"CI overlap differs in case {i}"



