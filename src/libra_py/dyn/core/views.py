"""
 no data duplication
 tensor rank transparent

"""

from .storage import TensorStorage

class TensorView:
    """
    Lightweight views into TensorStorage for a single TBF.
    """

    def __init__(self, storage: TensorStorage, idx: int):
        self._s = storage
        self._i = idx

    @property
    def R(self): return self._s.R[self._i]
    @property
    def P(self): return self._s.P[self._i]
    @property
    def C(self): return self._s.C[self._i]

    @property
    def H_adi(self): return self._s.H_adi[self._i]
    @property
    def H_dia(self): return self._s.H_dia[self._i]

    @property
    def forces(self): return self._s.forces[self._i]

