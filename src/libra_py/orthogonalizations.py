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
.. module:: orthogonalizations
   :platform: Unix, Windows
   :synopsis: this module implements various general-purpose orthogonalization procedures

.. moduleauthor:: Alexey V. Akimov, ChatGPT

"""

import numpy as np

def lowdin_inverse_sqrt(S, thresh=1e-10):
    """
    Compute the inverse square root of an overlap matrix using Löwdin symmetric orthogonalization.

    This function evaluates S^{-1/2} via eigenvalue decomposition:
        S = U diag(λ) U†  →  S^{-1/2} = U diag(λ^{-1/2}) U†

    Small eigenvalues below a specified threshold are treated as zero to ensure
    numerical stability and avoid divergences. The input matrix is explicitly
    symmetrized to enforce Hermiticity before diagonalization.

    Parameters
    ----------
    S : array_like (N, N)
        Real or complex overlap matrix. It is assumed to be Hermitian
        (or close to Hermitian within numerical precision).
    thresh : float, optional
        Eigenvalue cutoff below which eigenvalues are treated as zero.
        This prevents instabilities due to near-linear dependencies.
        Default is 1e-10.

    Returns
    -------
    S_inv_sqrt : ndarray (N, N)
        The inverse square root of the overlap matrix, S^{-1/2}, computed
        in the Löwdin (symmetric) orthogonalization scheme.

    Notes
    -----
    - The function enforces Hermiticity via (S + S†)/2 before diagonalization.
    - Input arrays are cast to float64 or complex128 to ensure compatibility
      with NumPy linear algebra routines.
    - Eigenvalues below `thresh` are set to zero, effectively projecting out
      linearly dependent components of the basis.
    - The resulting matrix can be used to orthonormalize a non-orthogonal basis:
          C_orth = S^{-1/2} C

    Raises
    ------
    LinAlgError
        If the eigenvalue decomposition fails.

    Examples
    --------
    >>> S = np.array([[1.0, 0.2], [0.2, 1.0]])
    >>> S_inv_sqrt = lowdin_inverse_sqrt(S)
    >>> np.allclose(S_inv_sqrt @ S @ S_inv_sqrt, np.eye(2))
    True
    """

    # Ensure numpy array
    S = np.array(S)

    # Force supported dtype
    if np.iscomplexobj(S):
        S = S.astype(np.complex128)
    else:
        S = S.astype(np.float64)

    # Enforce Hermiticity 
    S = 0.5 * (S + S.conj().T)

    eigvals, eigvecs = np.linalg.eigh(S)

    eigvals_inv_sqrt = np.zeros_like(eigvals)

    for i, val in enumerate(eigvals):
        if val > thresh:
            eigvals_inv_sqrt[i] = 1.0 / np.sqrt(val)
        else:
            eigvals_inv_sqrt[i] = 0.0

    S_inv_sqrt = eigvecs @ np.diag(eigvals_inv_sqrt) @ eigvecs.conj().T

    return S_inv_sqrt
