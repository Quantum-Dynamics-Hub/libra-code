from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ._helpers import merged_params, zeros_from
from .base import AnalyticalHamiltonianModel


@dataclass
class EschLevineLinearModel(AnalyticalHamiltonianModel):
    """
    General M-state linear diabatic crossing model.

    ``Q[0]=x``. Each matrix element is linear in the coordinate:
    ``H_ij(x) = V_ij + w_ij x`` and ``dH_ij/dx = w_ij``. The number of states
    is taken from ``params["nstates"]`` and the matrices ``V`` and ``w``.

    This is the flexible linear-crossing form in the legacy Esch-Levine module.
    Legacy source: ``libra_py.models.Esch_Levine.general``.
    """

    ndof: int = 1

    def diabatic_with_derivatives(self, Q, params):
        p = merged_params(params, {"nstates": 2, "V": [[0.0, 0.0], [0.0, 0.0]], "w": [[0.0, 0.0], [0.0, 0.0]]})
        x = Q[0]
        V = np.asarray(p["V"])
        w = np.asarray(p["w"])
        n = int(p["nstates"])
        self.nstates = n
        H = zeros_from(x, (*tuple(x.shape), n, n))
        dH = zeros_from(x, (*tuple(x.shape), 1, n, n))
        for i in range(n):
            for j in range(n):
                H[..., i, j] = V[i, j] + w[i, j] * x
                dH[..., 0, i, j] = w[i, j]
        return H, dH


@dataclass
class EschLevineJCP2020Model(AnalyticalHamiltonianModel):
    """
    Esch-Levine JCP 2020 multi-state crossing model.

    ``Q[0]=x``. State 0 has slope ``-w0``. States ``i>=1`` have slope ``w1``
    and offsets ``-i*delta`` with an optional extra shift ``eps`` for
    ``i >= i_crit``. State 0 is coupled uniformly to every other state by
    ``V``; other off-diagonal elements are zero.

    Parameters
    ----------
    ``w0``, ``w1``, ``V``, ``eps``, ``i_crit``, ``nstates``, and ``delta``.
    ``esch_levine_jcp2020_params`` returns the legacy parameter sets.

    Reference: M. P. Esch and B. G. Levine, J. Chem. Phys. 2020,
    153, 114104, https://doi.org/10.1063/5.0022529.
    Legacy source: ``libra_py.models.Esch_Levine``.
    """

    ndof: int = 1

    def diabatic_with_derivatives(self, Q, params):
        p = merged_params(params, {"w0": 0.25, "w1": 0.025, "V": 0.005, "eps": 0.0, "i_crit": 2, "nstates": 2, "delta": 0.01})
        x = Q[0]
        n = int(p["nstates"])
        self.nstates = n
        H = zeros_from(x, (*tuple(x.shape), n, n))
        dH = zeros_from(x, (*tuple(x.shape), 1, n, n))
        H[..., 0, 0] = -p["w0"] * x
        dH[..., 0, 0, 0] = -p["w0"]
        for i in range(1, n):
            shift = p["eps"] if i >= p["i_crit"] else 0.0
            H[..., i, i] = p["w1"] * x - i * p["delta"] - shift
            H[..., 0, i] = p["V"]
            H[..., i, 0] = p["V"]
            dH[..., 0, i, i] = p["w1"]
        return H, dH


def esch_levine_jcp2020_params(set_index):
    """
    Return parameter sets from the legacy Esch-Levine JCP 2020 examples.

    ``set_index`` in ``1..6`` selects the number of states, state spacing, and
    shifted manifold described in the original helper.
    """
    if set_index == 1:
        return {"w0": 0.25, "w1": 0.025, "V": 0.005, "eps": 0.0, "i_crit": 0, "nstates": 9, "delta": 0.01}
    if set_index == 2:
        return {"w0": 0.25, "w1": 0.025, "V": 0.005, "eps": 0.08, "i_crit": 5, "nstates": 9, "delta": 0.01}
    if set_index == 3:
        return {"w0": 0.25, "w1": 0.025, "V": 0.005, "eps": 0.0, "i_crit": 0, "nstates": 17, "delta": 0.005}
    if set_index == 4:
        return {"w0": 0.25, "w1": 0.025, "V": 0.005, "eps": 0.0, "i_crit": 0, "nstates": 17, "delta": 0.0025}
    if set_index == 5:
        return {"w0": 0.25, "w1": 0.025, "V": 0.005, "eps": 0.0, "i_crit": 0, "nstates": 17, "delta": 0.00125}
    if set_index == 6:
        return {"w0": 0.25, "w1": 0.025, "V": 0.005, "eps": 0.0, "i_crit": 0, "nstates": 17, "delta": 0.000625}
    raise ValueError("set_index must be in 1..6")
