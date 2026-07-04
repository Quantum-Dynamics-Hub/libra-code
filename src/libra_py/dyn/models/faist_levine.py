from __future__ import annotations

from dataclasses import dataclass

from libra_py.units import Angst, ev2au

from ._helpers import two_state_with_derivatives
from .base import AnalyticalHamiltonianModel


@dataclass
class FaistLevineModel(AnalyticalHamiltonianModel):
    """Faist-Levine two-state alkali-halogen collision model."""

    def diabatic_with_derivatives(self, Q, params):
        xp = self.xp
        R = Q[0]
        p = params
        inv = 1.0 / R
        p2 = inv * inv
        p4 = p2 * p2
        p6 = p4 * p2
        p8 = p6 * p2
        p12 = p6 * p6

        e = xp.exp(-R / p["rho_cov"])
        b12 = p["B_cov"] ** 12
        pb12 = b12 * p12
        h00 = (p["A_cov"] + pb12) * e - p["C_cov"] * p6
        dh00 = (-12.0 * pb12 / R) * e - (p["A_cov"] + pb12) * e / p["rho_cov"] + 6.0 * p["C_cov"] * p6 / R

        e = xp.exp(-R / p["rho_ion"])
        b8 = p["B_ion"] ** 8
        pb8 = b8 * p8
        alpha_sum = p["alp_M+"] + p["alp_X-"]
        alpha_prod = p["alp_M+"] * p["alp_X-"]
        h11 = (p["A_ion"] + pb8) * e - p["C_ion"] * p6 - inv - 0.5 * alpha_sum * p4 - 2.0 * alpha_prod * p6 / R + p["E_th"]
        dh11 = (-8.0 * pb8 / R) * e - (p["A_ion"] + pb8) * e / p["rho_ion"] + 6.0 * p["C_ion"] * p6 / R + p2 + 2.0 * alpha_sum * p4 / R + 14.0 * alpha_prod * p8

        h01 = p["A"] * xp.exp(-R / p["rho"])
        dh01 = -h01 / p["rho"]
        return two_state_with_derivatives(self, R, h00, h11, h01, dh00, dh11, dh01)


def faist_levine_nai_params():
    A = Angst
    eV = ev2au
    return {
        "A_cov": 3150.0 * eV,
        "A_ion": 2760.0 * eV,
        "B_cov": 2.647 * (eV ** (1.0 / 12.0)) * A,
        "B_ion": 2.398 * (eV ** (1.0 / 8.0)) * A,
        "C_cov": 1000.0 * eV * (A**6),
        "C_ion": 11.3 * eV * (A**6),
        "rho_cov": 0.435 * A,
        "rho_ion": 0.3489 * A,
        "alp_M+": 0.408 * (A**3),
        "alp_X-": 6.431 * (A**3),
        "E_th": 2.075 * eV,
        "A": 17.08 * eV,
        "rho": 1.239 * A,
    }


def faist_levine_lii_params():
    params = faist_levine_nai_params()
    A = Angst
    eV = ev2au
    params.update(
        {
            "A_ion": 1052.0 * eV,
            "B_cov": 2.996 * (eV ** (1.0 / 12.0)) * A,
            "B_ion": 1.839 * (eV ** (1.0 / 8.0)) * A,
            "C_cov": 1191.2 * eV * (A**6),
            "C_ion": 0.823 * eV * (A**6),
            "rho_cov": 0.44 * A,
            "rho_ion": 0.3786 * A,
            "alp_M+": 0.029 * (A**3),
            "E_th": 2.326 * eV,
        }
    )
    return params
