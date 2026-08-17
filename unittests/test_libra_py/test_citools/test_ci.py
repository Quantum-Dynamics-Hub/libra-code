import numpy as np

from libra_py.citools.ci import overlap


def _run_equivalent_ci_examples():
    results = []

    s = np.array([
        [1.0, 0.2],
        [0.2, 1.2],
    ])
    params = dict(
        S=np.kron(np.eye(2), s),
        nelec=2,
        nocc=0,
        nvirt=1,
        homo_indx=1,
        nstates=2,
        data=[[], [[[1, 2]]], [[1.0]]],
        active_space=None,
    )
    results.append(overlap(params["S"], params["data"], params["data"], params))

    s = np.array([
        [-1.0, 0.0, 0.0],
        [0.0, 1.0, 0.2],
        [0.0, 0.2, 1.2],
    ])
    params = dict(
        S=np.kron(np.eye(2), s),
        nelec=4,
        nocc=1,
        nvirt=1,
        homo_indx=2,
        nstates=2,
        data=[[], [[[2, 3]]], [[1.0]]],
        active_space=None,
    )
    results.append(overlap(params["S"], params["data"], params["data"], params))

    params = dict(
        S=np.kron(np.eye(2), s),
        nelec=2,
        nocc=1,
        nvirt=1,
        homo_indx=2,
        nstates=2,
        data=[[], [[[2, 3]]], [[1.0]]],
        active_space=[2, 3],
    )
    results.append(overlap(params["S"], params["data"], params["data"], params))

    return results


def test_equivalent_active_space_constructions_produce_same_ci_overlap():
    results = _run_equivalent_ci_examples()
    ref_ci = results[0]

    for st_ci in results[1:]:
        assert np.allclose(st_ci, ref_ci, atol=1e-12)


def test_triplet_ci_overlap_uses_spin_adapted_excited_state():
    params = dict(
        nelec=2,
        nocc=0,
        nvirt=1,
        homo_indx=1,
        nstates=2,
        active_space=None,
        spin=1,
        spin_projection=0,
    )
    data = [[], [[[1, 2]]], [[1.0]]]

    st_ci = overlap(np.eye(4), data, data, params)

    assert st_ci.shape == (2, 2)
    assert np.allclose(st_ci, np.eye(2), atol=1e-12)
