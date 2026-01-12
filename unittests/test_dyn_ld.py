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
from libra_py.dyn.local_diabatization import orthogonalized_T


# 1. Helper utilities (used in multiple tests)
def is_unitary(U, tol=1e-10):
    I = np.eye(U.shape[0], dtype=U.dtype)
    return np.linalg.norm(U.conj().T @ U - I) < tol


# 2. Identity input (sanity check)
# If T = I, the result must be exactly I.
def test_identity_matrix():
    n = 5
    T = np.eye(n, dtype=np.complex128)

    U = orthogonalized_T(T)

    assert np.allclose(U, np.eye(n))
    assert is_unitary(U)


# 3. Already-unitary matrix (should be unchanged)
# If T is unitary, orthogonalization should not modify it (up to numerical noise).
def test_unitary_input():
    rng = np.random.default_rng(0)
    n = 4

    # Random unitary via QR
    A = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    Q, _ = np.linalg.qr(A)

    U = orthogonalized_T(Q)

    assert is_unitary(U)
    assert np.allclose(U, Q, atol=1e-10)


# 4. Random non-unitary matrix
# This is the main expected use case.
def test_random_nonunitary_matrix():
    rng = np.random.default_rng(1)
    n = 6

    T = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    U = orthogonalized_T(T)

    # Must be unitary
    assert is_unitary(U)

    # U should span the same column space as T
    # i.e., T = U H for some Hermitian H
    H = U.conj().T @ T
    assert np.allclose(H, H.conj().T, atol=1e-10)


# 5. Phase consistency (diagonal matrix)
# For diagonal matrices, the result should be pure phases.
def test_diagonal_matrix():
    phases = np.exp(1j * np.array([0.1, 1.3, -2.0]))
    scales = np.array([2.0, 0.5, 3.0])

    T = np.diag(scales * phases)

    U = orthogonalized_T(T)

    # Result must be diagonal unitary with same phases
    assert is_unitary(U)
    assert np.allclose(np.diag(U), phases)



# 6. Gauge optimality (closest unitary)
# Checks Frobenius optimality property
def test_closest_unitary_property():
    rng = np.random.default_rng(3)
    n = 4

    T = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    U = orthogonalized_T(T)

    # Compare against a random unitary
    A = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    Q, _ = np.linalg.qr(A)

    dist_U = np.linalg.norm(T - U)
    dist_Q = np.linalg.norm(T - Q)

    assert dist_U <= dist_Q + 1e-10



# 7. Full-rank → unitary
def test_full_rank_unitary():
    rng = np.random.default_rng(0)
    n = 5
    T = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))

    U = orthogonalized_T(T)
    assert np.allclose(U.conj().T @ U, np.eye(n), atol=1e-10)


# 8. Near-singular → raises by default
def test_near_singular_raises():
    rng = np.random.default_rng(1)
    n = 5

    U0, _ = np.linalg.qr(rng.normal(size=(n, n)))
    V0, _ = np.linalg.qr(rng.normal(size=(n, n)))
    s = np.ones(n)
    s[2] = 1e-8

    T = U0 @ np.diag(s) @ V0.conj().T

    with pytest.raises(ValueError):
        orthogonalized_T(T, tol=1e-6)

# 9. Near-singular → reduced subspace
def test_near_singular_drop():
    rng = np.random.default_rng(2)
    n = 5

    U0, _ = np.linalg.qr(rng.normal(size=(n, n)))
    V0, _ = np.linalg.qr(rng.normal(size=(n, n)))
    s = np.ones(n)
    s[-1] = 1e-8

    T = U0 @ np.diag(s) @ V0.conj().T

    U, info = orthogonalized_T(T, tol=1e-6, drop=True)

    P = info["projector"]
    assert info["rank"] == 4
    assert np.allclose(P @ P, P, atol=1e-10)

