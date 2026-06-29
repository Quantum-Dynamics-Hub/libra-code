# *********************************************************************************
# * Copyright (C) 2026 Alexey V. Akimov
# *
# * This file is distributed under the terms of the GNU General Public License
# * as published by the Free Software Foundation, either version 3 of
# * the License, or (at your option) any later version.
# * See the file LICENSE in the root directory of this distribution
# * or <http://www.gnu.org/licenses/>.
# *
# *********************************************************************************/

import pytest
import numpy as np
from libra_py.initial_conditions import remove_translation, remove_rotation, cleanup_momenta


# ======================================================
# Helpers
# ======================================================

def total_linear_momentum(p):
    natoms = p.shape[0] // 3
    P = p.reshape(natoms, 3, -1)

    return np.sum(P, axis=0)


def total_angular_momentum(q, p, masses):

    ndof, ntraj = q.shape
    natoms = ndof // 3

    m = masses.reshape(natoms, 3)[:, 0]

    R = q.reshape(natoms, 3, ntraj)
    P = p.reshape(natoms, 3, ntraj)

    L_all = []

    for itraj in range(ntraj):

        r = R[:, :, itraj]

        r_com = (
            np.sum(m[:, None] * r, axis=0)
            / np.sum(m)
        )

        r = r - r_com

        mom = P[:, :, itraj]

        L = np.sum(np.cross(r, mom), axis=0)

        L_all.append(L)

    return np.array(L_all).T


def toy_system(ntraj=1):

    natoms = 3
    ndof = 9

    masses = np.repeat([1.0, 2.0, 3.0], 3)

    q = np.array([
        [0],[0],[0],
        [1],[0],[0],
        [0],[1],[0]
    ], dtype=float)

    q = np.tile(q, (1, ntraj))

    return q, masses


# ======================================================
# Translation tests
# ======================================================

def test_remove_translation_zero_total_momentum():

    q, masses = toy_system()

    p = np.array([
        [1],[0],[0],
        [2],[0],[0],
        [3],[0],[0]
    ], dtype=float)

    p2 = remove_translation(p)

    Ptot = total_linear_momentum(p2)

    assert np.allclose(Ptot, 0.0, atol=1e-12)


def test_remove_translation_preserves_shape():

    p = np.random.random((12, 7))

    out = remove_translation(p)

    assert out.shape == p.shape


def test_remove_translation_no_change_if_clean():

    p = np.array([
        [-1],[0],[0],
        [1],[0],[0]
    ] * 3)

    out = remove_translation(p)

    assert np.allclose(out, p)


# ======================================================
# Rotation tests
# ======================================================

def test_remove_rotation_zero_total_angular_momentum():

    q, masses = toy_system()

    p = np.array([
        [0],[1],[0],
        [0],[1],[0],
        [0],[1],[0]
    ], dtype=float)

    p2 = remove_rotation(q, p, masses)

    L = total_angular_momentum(q, p2, masses)

    assert np.allclose(L, 0.0, atol=1e-10)


def test_remove_rotation_preserves_shape():

    q = np.random.random((15, 4))
    p = np.random.random((15, 4))

    masses = np.repeat(np.ones(5), 3)

    out = remove_rotation(q, p, masses)

    assert out.shape == p.shape


# ======================================================
# Combined cleanup
# ======================================================

def test_cleanup_removes_translation_and_rotation():

    q, masses = toy_system()

    p = np.array([
        [1],[1],[0],
        [2],[1],[0],
        [3],[1],[0]
    ], dtype=float)

    out = cleanup_momenta(
        q,
        p,
        masses,
        stages=(
            "translation",
            "rotation",
            "translation"
        )
    )

    P = total_linear_momentum(out)
    L = total_angular_momentum(q, out, masses)

    assert np.allclose(P, 0.0, atol=1e-10)
    assert np.allclose(L, 0.0, atol=1e-10)


def test_cleanup_multiple_trajectories():

    ntraj = 5

    q, masses = toy_system(ntraj)

    rng = np.random.default_rng(123)

    p = rng.random((9, ntraj))

    out = cleanup_momenta(q, p, masses)

    P = total_linear_momentum(out)
    L = total_angular_momentum(q, out, masses)

    assert np.allclose(P, 0.0, atol=1e-10)
    assert np.allclose(L, 0.0, atol=1e-10)


# ======================================================
# Error handling
# ======================================================

def test_cleanup_unknown_stage():

    q = np.zeros((9, 1))
    p = np.zeros((9, 1))

    masses = np.ones(9)

    with pytest.raises(ValueError):
        cleanup_momenta(
            q,
            p,
            masses,
            stages=("banana",)
        )


