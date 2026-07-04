from __future__ import annotations

from dataclasses import dataclass

from libra_py.units import Angst, ev2Ha

from ._helpers import zeros_from
from .base import AnalyticalHamiltonianModel


@dataclass
class PhenolModel(AnalyticalHamiltonianModel):
    """
    Three-state two-dimensional phenol model of Pollien/Arribas/Agostini.

    Coordinates
    -----------
    ``Q[0] = r``
        O-H distance in Bohr.
    ``Q[1] = theta``
        Angular coordinate in radians.

    Math
    ----
    The diabatic Hamiltonian is a real symmetric 3x3 matrix:

    ``H = [[V11, V12, 0], [V12, V22, V23], [0, V23, V33]]``.

    The diagonal terms combine Morse-like O-H stretch potentials with angular
    modulation factors ``1 - cos(2 theta)``. Several avoided-crossing helper
    surfaces are smoothed by ``sqrt((Va-Vb)^2 + chi)`` terms. The off-diagonal
    couplings are radial switching functions multiplied by ``sin(theta)``:

    ``V12 = lambda12(r) sin(theta)``
    ``V23 = lambda23(r) sin(theta)``.

    The derivative tensor contains analytical derivatives with respect to
    ``r`` and ``theta``. The square-root derivative is protected at exact
    degeneracy by returning zero when the denominator is zero.

    Parameters
    ----------
    The current implementation uses the published constants from the legacy
    ``libra_py.models.Phenol.Pollien_Arribas_Agostini`` function. User
    parameter overrides are not yet exposed because the original function also
    hard-coded this parameter set.

    References
    ----------
    A. Pollien, E. Villaseco Arribas, D. Lauvergnat, and F. Agostini,
    "Exact-Factorisation Study of the Photochemistry of Phenol",
    Molecular Physics, e2378960,
    https://doi.org/10.1080/00268976.2024.2378960.
    """

    ndof: int = 2
    nstates: int = 3

    def diabatic_with_derivatives(self, Q, params):
        xp = self.xp
        r, theta = Q[0], Q[1]
        shape = tuple(r.shape)
        H = zeros_from(r, (*shape, 3, 3))
        dH = zeros_from(r, (*shape, 2, 3, 3))

        D1e = 4.26302 * ev2Ha
        D3e = 4.47382 * ev2Ha
        a1 = 2.66021 / Angst
        a3 = 2.38671 / Angst
        r1 = 0.96944 * Angst
        r3 = 0.96304 * Angst
        a30 = 4.85842 * ev2Ha

        A1 = 0.27037 * ev2Ha
        A2 = 1.96606 * Angst
        A3 = 0.685264 * Angst
        C1 = 0.110336 * ev2Ha
        C2 = 1.21724 * Angst
        C3 = 0.06778 * Angst

        B201 = 0.192205 * ev2Ha
        B202 = 5.67356 / Angst
        B203 = 1.03171 * Angst
        B204 = 5.50696 * ev2Ha
        B205 = 4.70601 * ev2Ha
        B206 = 2.49826 / Angst
        B207 = 0.988188 * Angst
        B208 = 3.3257 * ev2Ha

        lambda_12_max = 1.47613 * ev2Ha
        d12 = 1.96984 * Angst
        beta12 = 0.494373 * Angst
        lambda_23_max = 0.327204 * ev2Ha
        d23 = 1.22594 * Angst
        beta23 = 0.0700604 * Angst

        chi20 = 0.326432 * ev2Ha**2
        chi21 = 0.021105 * ev2Ha**2
        chi22 = 0.0

        B211 = -0.2902 * ev2Ha
        B212 = 2.05715 * Angst
        B213 = 1.01574 * Angst
        B214 = -73.329 * ev2Ha
        B215 = 1.48285 * Angst
        B216 = -0.1111 * Angst
        B217 = -0.00055 * ev2Ha

        B221 = 27.3756 * ev2Ha
        B222 = 1.66881 * Angst
        B223 = 0.20557 * Angst
        B224 = 0.35567 * ev2Ha
        B225 = 1.43492 * Angst
        B226 = 0.56968 * Angst

        x = xp.exp(-a1 * (r - r1))
        dx_dr = -a1 * x
        v10 = D1e * (1.0 - x) ** 2
        dv10_dr = -2.0 * D1e * (1.0 - x) * dx_dr

        x = xp.tanh((r - A2) / A3)
        dx_dr = (1.0 - x * x) / A3
        v11 = 0.5 * A1 * (1.0 - x)
        dv11_dr = -0.5 * A1 * dx_dr

        x = xp.exp(-B202 * (r - B203))
        dx_dr = -B202 * x
        v201 = B201 * (1.0 - x) ** 2 + B204
        dv201_dr = -2.0 * B201 * (1.0 - x) * dx_dr

        x = xp.exp(-B206 * (r - B207))
        dx_dr = -B206 * x
        v202 = B205 * x + B208
        dv202_dr = B205 * dx_dr

        x = xp.sqrt((v201 - v202) ** 2 + chi20)
        dx_dr = (v201 - v202) * (dv201_dr - dv202_dr) / x
        v20 = 0.5 * (v201 + v202) - 0.5 * x
        dv20_dr = 0.5 * (dv201_dr + dv202_dr) - 0.5 * dx_dr

        x = xp.tanh((r - B212) / B213)
        dx_dr = (1.0 - x * x) / B213
        v211 = 0.5 * B211 * (1.0 - x)
        dv211_dr = -0.5 * B211 * dx_dr

        x = xp.tanh((r - B215) / B216)
        dx_dr = (1.0 - x * x) / B216
        v212 = 0.5 * B214 * (1.0 - x) + B217
        dv212_dr = -0.5 * B214 * dx_dr

        x = xp.sqrt((v211 - v212) ** 2 + chi21)
        dx_dr = (v211 - v212) * (dv211_dr - dv212_dr) / x
        v21 = 0.5 * (v211 + v212) + 0.5 * x
        dv21_dr = 0.5 * (dv211_dr + dv212_dr) + 0.5 * dx_dr

        x = xp.tanh((r - B222) / B223)
        dx_dr = (1.0 - x * x) / B223
        v221 = 0.5 * B221 * (1.0 + x)
        dv221_dr = 0.5 * B221 * dx_dr

        x = xp.tanh((r - B225) / B226)
        dx_dr = (1.0 - x * x) / B226
        v222 = 0.5 * B224 * (1.0 - x)
        dv222_dr = -0.5 * B224 * dx_dr

        x = xp.sqrt((v221 - v222) ** 2 + chi22)
        dx_dr = _safe_divide(xp, (v221 - v222) * (dv221_dr - dv222_dr), x)
        v22 = 0.5 * (v221 + v222) - 0.5 * x
        dv22_dr = 0.5 * (dv221_dr + dv222_dr) - 0.5 * dx_dr

        x = xp.exp(-a3 * (r - r3))
        dx_dr = -a3 * x
        v30 = D3e * (1.0 - x) ** 2 + a30
        dv30_dr = -2.0 * D3e * (1.0 - x) * dx_dr

        x = xp.tanh((r - C2) / C3)
        dx_dr = (1.0 - x * x) / C3
        v31 = 0.5 * C1 * (1.0 - x)
        dv31_dr = -0.5 * C1 * dx_dr

        cs = xp.cos(theta)
        si = xp.sin(theta)
        cs2 = xp.cos(2.0 * theta)
        si2 = xp.sin(2.0 * theta)

        x = xp.tanh((r - d12) / beta12)
        dx_dr = (1.0 - x * x) / beta12
        lambda_12_r = 0.5 * lambda_12_max * (1.0 - x)
        dlambda_12_r_dr = -0.5 * lambda_12_max * dx_dr
        V12 = lambda_12_r * si
        dV12_dr = dlambda_12_r_dr * si
        dV12_dtheta = lambda_12_r * cs

        V13 = 0.0 * r
        dV13_dr = 0.0 * r
        dV13_dtheta = 0.0 * r

        x = xp.tanh((r - d23) / beta23)
        dx_dr = (1.0 - x * x) / beta23
        lambda_23_r = 0.5 * lambda_23_max * (1.0 - x)
        dlambda_23_r_dr = -0.5 * lambda_23_max * dx_dr
        V23 = lambda_23_r * si
        dV23_dr = dlambda_23_r_dr * si
        dV23_dtheta = lambda_23_r * cs

        bend = 1.0 - cs2
        V11 = v10 + v11 * bend
        V22 = v20 + v21 * bend + v22 * bend * bend
        V33 = v30 + v31 * bend
        dV11_dr = dv10_dr + dv11_dr * bend
        dV11_dtheta = 2.0 * v11 * si2
        dV22_dr = dv20_dr + dv21_dr * bend + dv22_dr * bend * bend
        dV22_dtheta = 2.0 * v21 * si2 + 4.0 * v22 * bend * si2
        dV33_dr = dv30_dr + dv31_dr * bend
        dV33_dtheta = 2.0 * v31 * si2

        _fill_three_state(H, V11, V22, V33, V12, V13, V23)
        _fill_three_state(dH[..., 0, :, :], dV11_dr, dV22_dr, dV33_dr, dV12_dr, dV13_dr, dV23_dr)
        _fill_three_state(
            dH[..., 1, :, :],
            dV11_dtheta,
            dV22_dtheta,
            dV33_dtheta,
            dV12_dtheta,
            dV13_dtheta,
            dV23_dtheta,
        )

        return H, dH


def _fill_three_state(matrix, v11, v22, v33, v12, v13, v23):
    matrix[..., 0, 0] = v11
    matrix[..., 1, 1] = v22
    matrix[..., 2, 2] = v33
    matrix[..., 0, 1] = v12
    matrix[..., 1, 0] = v12
    matrix[..., 0, 2] = v13
    matrix[..., 2, 0] = v13
    matrix[..., 1, 2] = v23
    matrix[..., 2, 1] = v23


def _safe_divide(xp, numerator, denominator):
    if xp.__name__ == "torch":
        return xp.where(denominator != 0.0, numerator / denominator, xp.zeros_like(numerator))
    return xp.divide(numerator, denominator, out=xp.zeros_like(numerator), where=denominator != 0.0)
