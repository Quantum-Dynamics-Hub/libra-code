
import numpy as np
from scipy.linalg import cholesky, solve_triangular
from scipy.linalg import eigh, eig

def generalized_eigensolve_scipy(A, B, hermitian=True, sort=True):
    """
    Solve the generalized eigenvalue problem A v = λ B v using SciPy.

    Args:
        A (np.ndarray): Matrix A (N x N)
        B (np.ndarray): Matrix B (N x N)
        hermitian (bool): Set True if A and B are Hermitian/symmetric.
        sort (bool): If True, sort eigenpairs by ascending eigenvalue.

    Returns:
        eigvals (np.ndarray): Eigenvalues (N,)
        eigvecs (np.ndarray): Corresponding eigenvectors (N, N)

    Example:

    import numpy as np

    # Symmetric A and B
    A = np.array([[6., 2.], [2., 3.]])
    B = np.array([[1., 0.], [0., 2.]])
    eigvals, eigvecs = generalized_eigensolve_scipy(A, B)
    print("Eigenvalues:", eigvals)
    print("Eigenvectors (columns):\n", eigvecs)

    """

    if hermitian:
        eigvals, eigvecs = eigh(A, B)
    else:
        eigvals, eigvecs = eig(A, B)

    if sort:
        idx = np.argsort(eigvals.real if not hermitian else eigvals)
        eigvals = eigvals[idx]
        eigvecs = eigvecs[:, idx]

    # Normalize eigenvectors with respect to B: v† B v = 1
    for i in range(eigvecs.shape[1]):
        vec = eigvecs[:, i]
        #norm = np.sqrt(np.conj(vec).T @ B @ vec)
        norm = np.sqrt(np.vdot(vec, B @ vec))  # vdot does conjugation
        eigvecs[:, i] = vec / norm


    return eigvals, eigvecs



def generalized_eigensolve_scipy_cholesky(A, B, hermitian=True, sort=True):
    """
    Robust solver for A v = λ B v, with B-orthonormal eigenvectors:
        v_i^† B v_j = δ_ij (up to numerical error)
    
    For Hermitian A, B > 0.
    """
    if not hermitian:
        raise NotImplementedError("Only Hermitian positive-definite B is supported in stable form")

    # Cholesky decomposition: B = L L^H
    L = cholesky(B, lower=True)

    # Transform to standard eigenvalue problem: C = L^{-1} A L^{-H}
    Linv_A = solve_triangular(L, A, lower=True, trans='N')
    C = solve_triangular(L, Linv_A.T, lower=True, trans='C').T

    # Solve C y = λ y
    eigvals, y = eigh(C)

    if sort:
        idx = np.argsort(eigvals)
        eigvals = eigvals[idx]
        y = y[:, idx]

    # Recover v = L^{-H} y
    eigvecs = solve_triangular(L, y, lower=True, trans='C')

    # At this point, eigvecs should satisfy v.T.conj() @ B @ v ≈ I
    return eigvals, eigvecs


"""

Tests and conventions for this module:

from liblibra_core import *
from libra_py import eigensolvers, data_conv
import numpy as np

# ===== Do the solution ==========
A = np.array([[0, 1.0], 
              [1.0, 1.0]])
B = np.array([[1.0, 0.0], 
              [0.0, 1.0]])

eigvals, eigvecs = eigensolvers.generalized_eigensolve_scipy(A, B, hermitian=False)


# ======= Numpy multiplications ==========
E = np.diag(eigvals)
C = eigvecs
print( E )
print( C )
print("LHS = ", A @ C)
print("RHS = ", B @ C @ E)

Selection deleted
E = np.diag(eigvals)
C = eigvecs
print( E )
print( C )
print("LHS = ", A @ C)

print("RHS = ", B @ C @ E)

>>[[-0.61803399+0.j  0.        +0.j]
>> [ 0.        +0.j  1.61803399+0.j]]
>>[[ 0.85065081  0.52573111]
>> [-0.52573111  0.85065081]]
>>LHS =  [[-0.52573111  0.85065081]
>> [ 0.3249197   1.37638192]]
>>RHS =  [[-0.52573111+0.j  0.85065081+0.j]
>> [ 0.3249197 +0.j  1.37638192+0.j]]


# ======= Conversion to MATRIX ===========

H = data_conv.nparray2MATRIX(A)
U = data_conv.nparray2MATRIX(C)
S = data_conv.nparray2MATRIX(B)
e = data_conv.nparray2MATRIX(E)

# ======== Libra multiplications ========

(H * U).show_matrix()

>>-0.52573111  0.85065081  
>>0.32491970  1.3763819   

(S * U * e).show_matrix()

>>-0.52573111  0.85065081
>>0.32491970  1.3763819

# ========= Eigenvectors ==============
print("First eigenvector")
print(U.get(0, 0))
print(U.get(1, 0))

>>First eigenvector
>>0.8506508083520399
>>-0.5257311121191337

print("First eigenvector")
print(C[0, 0])
print(C[1, 0])

>>First eigenvector
>>0.8506508083520399
>>-0.5257311121191337


"""

