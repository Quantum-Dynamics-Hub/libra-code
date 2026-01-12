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

"""
.. module:: local_diabatization
   :platform: Unix, Windows
   :synopsis: Implements various functions needed for local diabatization
.. moduleauthor:: Alexey V. Akimov

"""
import numpy as np

def orthogonalized_T(T, tol=1e-10, drop=False):
    """
    Orthonormalize a matrix T via (right) polar decomposition.

    Parameters
    ----------
    T : (N, N) complex ndarray
        Typically T = P^{-1}(t, t+dt), where P is an adiabatic overlap matrix.
    tol : float
        Eigenvalue cutoff for detecting ill-conditioning.
    drop : bool, optional
        If False (default), require T to be full rank and return a unitary matrix.
        If True, drop ill-conditioned states and return a partial isometry.

    Returns
    -------
    U : (N, N) complex ndarray
        Orthonormalized transformation matrix.
        - If drop=False: U is unitary (U† U = I)
        - If drop=True:  U† U is a projector

    info : dict (only if drop=True)
        Dictionary with diagnostic information:
        - 'rank'       : number of retained states
        - 'projector'  : U† U
        - 'eigvals'    : eigenvalues of T† T

    Raises
    ------
    ValueError
        If T is ill-conditioned and drop=False.


    Notes
    -----
    Is meant to replace: `CMATRIX orthogonalized_T(CMATRIX& T)` from
    dyn/dyn_variables_electronic.cpp


    Related theory
    --------------

    We are given a generally non-unitary matrix

        X = P^{-1}(t, t+dt)

    where:
        P(t, t+dt) = ⟨ψ_adi(t) | ψ_adi(t+dt)⟩

    is the time-overlap matrix between adiabatic states at consecutive time steps.

    We want to define a transformed adiabatic basis:

        |ψ̃_adi(t+dt)⟩ = |ψ_adi(t+dt)⟩ T

    such that the overlap with the previous basis is exactly the identity:

        ⟨ψ̃_adi(t) | ψ̃_adi(t+dt)⟩ = I

    This formally suggests T ≈ P^{-1}, but P^{-1} is not guaranteed to be unitary.
    Therefore, we seek a *unitary matrix unitarily equivalent to P^{-1}*.

    ### Polar decomposition
    Let:

        S = X† X

    which is Hermitian and positive-definite.

    Define:

        U = X S^{-1/2}

    Then:

        U† U = S^{-1/2} X† X S^{-1/2} = I

    Thus, U is unitary.

    This construction corresponds to the right polar decomposition:

        X = U H
        H = (X† X)^{1/2}

    Among all unitary matrices, U is the closest to X
    in the Frobenius norm sense.

    ### Final result

    The function returns:

        T_orth = P^{-1} (P^{-1†} P^{-1})^{-1/2}

    which is unitary and preserves the adiabatic overlap structure
    while eliminating numerical drift.

    ### Handling ill-conditioned or rank-deficient inputs

    In practical nonadiabatic dynamics, the overlap matrix P(t, t+dt) may
    become nearly singular due to:
      - large nuclear time steps,
      - near-degeneracies or avoided crossings,
      - numerical noise in the electronic structure.

    As a result, T = P^{-1} may be ill-conditioned or rank-deficient, and
    the inverse square root (T† T)^{-1/2} may not exist.

    This function supports two distinct behaviors controlled by `drop`:

    1. drop = False (default)
       ----------------------
       The function *requires* T to be full rank.

       If one or more eigenvalues of T† T are smaller than `tol`, a
       ValueError is raised. This corresponds to a physically fatal loss
       of linear independence in the adiabatic basis and mirrors the
       strict behavior of the original C++ implementation.

       In this mode, the returned matrix U is strictly unitary:

           U† U = I.

    2. drop = True
       -------------
       Ill-conditioned states with eigenvalues λ < tol of T† T are removed
       from the orthonormalization.

       The transformation is constructed only in the well-conditioned
       subspace, leading to a partial isometry:

           U† U = Π,

       where Π is the projector onto the retained (well-conditioned)
       subspace.

       This mode explicitly reduces the effective electronic subspace and
       should only be used when state pruning or adaptive basis reduction
       is intended.

    """

    # Hermitian overlap
    S = T.conj().T @ T

    eigvals, eigvecs = np.linalg.eigh(S)

    # Identify well-conditioned subspace
    mask = eigvals > tol
    rank = int(np.sum(mask))

    if rank < T.shape[0] and not drop:
        raise ValueError(
            "orthogonalized_T: T is rank-deficient or ill-conditioned; "
            "set drop=True to allow state reduction"
        )

    # Construct S^{-1/2} (full or reduced)
    if drop:
        V = eigvecs[:, mask]
        D_inv_sqrt = np.diag(eigvals[mask] ** (-0.5))
        S_inv_sqrt = V @ D_inv_sqrt @ V.conj().T
    else:
        S_inv_sqrt = eigvecs  @ np.diag(eigvals ** (-0.5)) @ eigvecs.conj().T

    U = T @ S_inv_sqrt

    # Diagnostics
    if drop:
        info = {
            "rank": rank,
            "projector": U.conj().T @ U,
            "eigvals": eigvals,
        }
        return U, info

    return U



