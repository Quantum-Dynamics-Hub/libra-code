from __future__ import annotations

from dataclasses import dataclass

from libra_py.units import Angst, ev2Ha

from ._helpers import merged_params, one_state_with_derivatives
from .base import AnalyticalHamiltonianModel


@dataclass
class BeswickJortnerModel(AnalyticalHamiltonianModel):
    """
    Beswick-Jortner one-state ICN dissociation model.

    Coordinates
    -----------
    ``Q[0] = r``
        C-N bond length.
    ``Q[1] = R``
        Distance between iodine and the CN center of mass.

    Math
    ----
    The single diabatic state is

    ``H = T0 + 1/2 K (r-r0)^2 + A exp[a (m_C/(m_C+m_N) r - R)]``.

    The first term is an electronic offset, the second is the CN stretch, and
    the exponential term is the short-range I-CN repulsion. The returned
    derivative tensor contains ``dH/dr`` and ``dH/dR``.

    Parameters
    ----------
    ``K``, ``r0``, ``T0``, ``A``, ``a``, ``m_c``, ``m_n`` may be overridden.
    Defaults are converted to atomic units from the legacy
    ``libra_py.models.Beswick_Jortner`` constants. The old file used undefined
    mass symbols, so this translation supplies explicit ``m_c=12`` and
    ``m_n=14`` defaults.

    References
    ----------
    A. Beswick and J. Jortner, Chem. Phys. 1977, 24, 1.
    R. C. Brown and E. J. Heller, J. Chem. Phys. 1981, 75, 186-188.
    """

    ndof: int = 2
    nstates: int = 1

    def diabatic_with_derivatives(self, Q, params):
        xp = self.xp
        p = merged_params(
            params,
            {
                "K": 74.4434 * ev2Ha / (Angst * Angst),
                "r0": 1.2327 * Angst,
                "T0": 4.36994 * ev2Ha,
                "A": 200000.0 * ev2Ha,
                "a": 6.68 / Angst,
                "m_c": 12.0,
                "m_n": 14.0,
            },
        )
        r, R = Q[0], Q[1]
        cn_fraction = p["m_c"] / (p["m_c"] + p["m_n"])
        expo = xp.exp(p["a"] * (cn_fraction * r - R))
        repulsive = p["A"] * expo
        value = p["T0"] + 0.5 * p["K"] * (r - p["r0"]) ** 2 + repulsive
        d_dr = p["K"] * (r - p["r0"]) + p["a"] * cn_fraction * repulsive
        d_dR = -p["a"] * repulsive

        return one_state_with_derivatives(self, value, (d_dr, d_dR))
