"""Optional PyTorch facade for TensorStorage and HamiltonianEngine.

Importing dyn does not import Torch. Construct this facade explicitly to select
float32/float64 and a CPU/CUDA device; complex allocations follow that precision.
"""
from __future__ import annotations


class TorchBackend:
    def __init__(self, device='cpu', dtype='float64'):
        import torch
        self.xp = torch
        self.device = torch.device(device)
        self.dtype = getattr(torch, dtype) if isinstance(dtype, str) else dtype
        if self.dtype not in (torch.float32, torch.float64):
            raise ValueError('TorchBackend requires float32 or float64')

    def _dtype(self, dtype):
        torch = self.xp
        if dtype is float:
            return self.dtype
        if dtype is complex:
            return torch.complex128 if self.dtype == torch.float64 else torch.complex64
        if dtype is int:
            return torch.int64
        if dtype is bool:
            return torch.bool
        return dtype

    def asarray(self, values, dtype=None):
        if dtype is None and not hasattr(values, 'dtype'):
            inferred = self.xp.as_tensor(values)
            if inferred.is_floating_point():
                dtype = float
            elif inferred.is_complex():
                dtype = complex
        return self.xp.as_tensor(values, dtype=self._dtype(dtype) if dtype is not None else None,
                                 device=self.device)

    def array(self, values, dtype=None):
        return self.asarray(values, dtype).clone()

    def zeros(self, shape, dtype=float):
        return self.xp.zeros(shape, dtype=self._dtype(dtype), device=self.device)

    def ones(self, shape, dtype=float):
        return self.xp.ones(shape, dtype=self._dtype(dtype), device=self.device)

    def eye(self, n, dtype=float):
        return self.xp.eye(n, dtype=self._dtype(dtype), device=self.device)

    def diag(self, values):
        return self.xp.diag(values)

    def einsum(self, equation, *operands):
        # NumPy promotes mixed real/complex operands; Torch contractions
        # require a common dtype (notably velocity times complex DC1).
        dtype = operands[0].dtype
        for operand in operands[1:]:
            dtype = self.xp.promote_types(dtype, operand.dtype)
        return self.xp.einsum(equation, *(x.to(dtype=dtype) for x in operands))

    def matmul(self, a, b):
        return a @ b

    def inverse(self, a):
        return self.xp.linalg.inv(a)

    def solve(self, a, b):
        return self.xp.linalg.solve(a, b)

    def conjugate_transpose(self, a):
        return a.conj().swapaxes(-1, -2)

    def expm(self, a):
        return self.xp.linalg.matrix_exp(a)
