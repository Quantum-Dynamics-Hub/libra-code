from __future__ import annotations

from dataclasses import dataclass

from ._helpers import merged_params, two_state_with_derivatives, zeros_from
from .base import AnalyticalHamiltonianModel


@dataclass
class GranucciPersicoModel1(AnalyticalHamiltonianModel):
    """
    One-dimensional Granucci-Persico avoided-crossing model.

    Coordinates
    -----------
    ``Q[0] = x`` is a one-dimensional nuclear coordinate in Bohr.

    Math
    ----
    ``H00 = a1 exp(-alp1 x) + dE``
    ``H11 = a2 exp(-alp2 x)``
    ``H01 = H10 = b exp[-beta (x-x_c)^2] + gamma sin^2(x)``

    The returned derivative tensor stores the analytical derivatives of these
    three scalar functions with respect to ``x``. The diabatic overlap is the
    identity and diabatic derivative couplings are zero, as in the original
    model interface.

    Parameters
    ----------
    ``a1``, ``a2``, ``alp1``, ``alp2``, ``dE``, ``b``, ``beta``, ``gamma``,
    and ``x_c`` follow the legacy defaults in atomic units.

    References
    ----------
    G. Granucci and M. Persico, J. Chem. Phys. 2007, 126, 134114.
    """

    ndof: int = 1
    nstates: int = 2

    def diabatic_with_derivatives(self, Q, params):
        xp = self.xp
        p = merged_params(
            params,
            {
                "a1": 0.006,
                "a2": 0.5,
                "alp1": 0.2,
                "alp2": 0.5,
                "dE": 0.03,
                "b": 0.013,
                "beta": 0.2,
                "gamma": 0.0002,
                "x_c": 5.0,
            },
        )
        x = Q[0]
        e1 = p["a1"] * xp.exp(-p["alp1"] * x)
        e2 = p["a2"] * xp.exp(-p["alp2"] * x)
        e3 = p["b"] * xp.exp(-p["beta"] * (x - p["x_c"]) ** 2)
        sin_x = xp.sin(x)
        cos_x = xp.cos(x)

        h00 = e1 + p["dE"]
        h11 = e2
        h01 = e3 + p["gamma"] * sin_x * sin_x
        dh00 = -p["alp1"] * e1
        dh11 = -p["alp2"] * e2
        dh01 = -2.0 * p["beta"] * (x - p["x_c"]) * e3 + 2.0 * p["gamma"] * sin_x * cos_x

        return two_state_with_derivatives(self, x, h00, h11, h01, dh00, dh11, dh01)


@dataclass
class GranucciPersicoModel2(AnalyticalHamiltonianModel):
    """
    Two-dimensional Granucci-Persico-Zoccante conical-intersection model.

    Coordinates
    -----------
    ``Q[0] = x`` and ``Q[1] = y``.

    Math
    ----
    The diagonal states are Morse-like functions plus a shared harmonic
    ``y`` term:

    ``H00 = D1 [exp(-2 alp1 (x-x1)) - 2 exp(-alp1 (x-x1))] + delta1 + 1/2 K y^2``
    ``H11 = D2 [exp(-2 alp2 (x-x2)) - 2 exp(-alp2 (x-x2))] + delta2 + 1/2 K y^2``

    The coupling is odd in ``y``:

    ``H01 = H10 = gamma y exp[-beta1 (x-x3)^2 - beta2 y^2]``.

    Both ``x`` and ``y`` derivatives are returned. The legacy Python source had
    obvious transcription bugs in this routine (undefined ``beta`` and only one
    derivative matrix); this translation follows the documented equations.

    References
    ----------
    G. Granucci, M. Persico, and A. Zoccante, J. Chem. Phys. 2010,
    133, 134111.
    """

    ndof: int = 2
    nstates: int = 2

    def diabatic_with_derivatives(self, Q, params):
        xp = self.xp
        p = merged_params(
            params,
            {
                "D1": 0.015,
                "D2": 0.11,
                "delta1": 0.05,
                "delta2": 0.11,
                "alp1": 1.0,
                "alp2": 0.674,
                "beta1": 0.5,
                "beta2": 1.5,
                "gamma": 2.2295e-2,
                "K": 0.09,
                "x1": 3.9,
                "x2": 3.0,
                "x3": 5.0,
            },
        )
        x, y = Q[0], Q[1]
        shape = tuple(x.shape)
        H = zeros_from(x, (*shape, 2, 2))
        dH = zeros_from(x, (*shape, 2, 2, 2))

        e1 = xp.exp(-p["alp1"] * (x - p["x1"]))
        e2 = xp.exp(-p["alp2"] * (x - p["x2"]))
        e3 = xp.exp(-p["beta1"] * (x - p["x3"]) ** 2 - p["beta2"] * y * y)
        h00 = p["D1"] * (e1 * e1 - 2.0 * e1) + p["delta1"] + 0.5 * p["K"] * y * y
        h11 = p["D2"] * (e2 * e2 - 2.0 * e2) + p["delta2"] + 0.5 * p["K"] * y * y
        h01 = p["gamma"] * y * e3

        H[..., 0, 0] = h00
        H[..., 1, 1] = h11
        H[..., 0, 1] = h01
        H[..., 1, 0] = h01

        de1_dx = -p["alp1"] * e1
        de2_dx = -p["alp2"] * e2
        dh00_dx = 2.0 * p["D1"] * (e1 - 1.0) * de1_dx
        dh11_dx = 2.0 * p["D2"] * (e2 - 1.0) * de2_dx
        dh01_dx = -2.0 * p["beta1"] * (x - p["x3"]) * h01

        dh00_dy = p["K"] * y
        dh11_dy = dh00_dy
        dh01_dy = p["gamma"] * e3 - 2.0 * p["beta2"] * y * h01

        dH[..., 0, 0, 0] = dh00_dx
        dH[..., 0, 1, 1] = dh11_dx
        dH[..., 0, 0, 1] = dh01_dx
        dH[..., 0, 1, 0] = dh01_dx
        dH[..., 1, 0, 0] = dh00_dy
        dH[..., 1, 1, 1] = dh11_dy
        dH[..., 1, 0, 1] = dh01_dy
        dH[..., 1, 1, 0] = dh01_dy

        return H, dH
