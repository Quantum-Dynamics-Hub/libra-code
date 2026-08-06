from __future__ import annotations

import numpy as np


class NumpyBackend:
    """
    Small NumPy facade used by the new dyn prototype.

    The class intentionally mirrors the array operations used by the
    backend-agnostic layers without hiding the underlying ndarray type.
    """

    array = staticmethod(np.array)
    asarray = staticmethod(np.asarray)
    zeros = staticmethod(np.zeros)
    ones = staticmethod(np.ones)
    eye = staticmethod(np.eye)
    diag = staticmethod(np.diag)
    einsum = staticmethod(np.einsum)
    matmul = staticmethod(np.matmul)
    inverse = staticmethod(np.linalg.inv)
    solve = staticmethod(np.linalg.solve)

    @staticmethod
    def conjugate_transpose(x):
        return np.swapaxes(np.conjugate(x), -1, -2)

    @staticmethod
    def expm(x):
        """
        Matrix exponential for one matrix or a batch of square matrices.

        This uses an eigendecomposition-based implementation to avoid adding
        a hard dependency on SciPy in the lightweight backend scaffold.
        """

        values, vectors = np.linalg.eig(x)
        exp_values = np.exp(values)
        scaled_vectors = vectors * exp_values[..., None, :]
        return np.matmul(scaled_vectors, np.linalg.inv(vectors))
