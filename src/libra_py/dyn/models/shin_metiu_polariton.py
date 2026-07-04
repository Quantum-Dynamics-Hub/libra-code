from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .shin_metiu import (
    ShinMetiuDVRData,
    bundled_shin_metiu_dvr_path,
    load_shin_metiu_dvr_data,
)


@dataclass
class ShinMetiuPolaritonModel:
    """
    Shin-Metiu molecule in a quantized single-mode optical cavity.

    The model starts from cavity-free Shin-Metiu DVR quantities in the
    electronic adiabatic basis and constructs a light-matter Hamiltonian in an
    electronic-state/Fock-state product basis. The coordinate is the nuclear
    Shin-Metiu coordinate ``R = Q[0]``.

    The cavity Hamiltonian follows the Pauli-Fierz long-wavelength form used
    in the legacy source:

    ``H = H_el + omega_c (a^dagger a + 1/2) + g_c epsilon mu (a^dagger + a)
          + (g_c^2 epsilon^2 / omega_c) mu^2``.

    ``g_c`` is the light-matter coupling strength, ``omega_c`` is the photon
    frequency, ``epsilon`` is the scalar polarization projection, ``mu`` is the
    Shin-Metiu dipole matrix, and the last term is the dipole self-energy. All
    quantities are in atomic units.

    ``model='2-state'`` uses the resonance subspace ``|e,0>`` and ``|g,1>``.
    ``model='4-state'`` uses ``|g,0>``, ``|e,0>``, ``|g,1>``, ``|e,1>``. The
    4-state derivative coupling matrix duplicates the electronic NAC inside
    each photon block. ``dH_adi`` follows the covariant derivative convention
    used by the Hamiltonian engine, and ``DC1_adi`` stores derivative
    couplings, not velocity-projected NACs.

    Legacy source: ``libra_py.models.Shin_Metiu.polariton``.
    Reference: J. Chem. Phys. 157, 104118 (2022), Eqs. 34-40.
    """

    params: dict | None = None
    dvr_data: ShinMetiuDVRData | None = None
    dvr_file: str | None = None
    model: str = "4-state"
    ndof: int = 1

    def __post_init__(self):
        self.params = dict(self.params or {})
        if self.dvr_data is None:
            filename = self.dvr_file
            if filename is None:
                source_model = self.params.get("shin_metiu_model", 1)
                filename = bundled_shin_metiu_dvr_path(source_model)
            self.dvr_data = load_shin_metiu_dvr_data(filename)
        if self.model not in ("2-state", "4-state"):
            raise ValueError("model must be '2-state' or '4-state'")
        self.nstates = 2 if self.model == "2-state" else 4

    def evaluate(self, Q, params=None):
        effective = dict(self.params)
        if params:
            effective.update(params)
        R = np.asarray(Q, dtype=float)[0]
        flat = R.reshape(-1)
        values = [
            polariton_info(
                point,
                self.dvr_data,
                model=self.model,
                g_c=effective.get("g_c", 0.005),
                omega_c=effective.get("omega_c", 0.1),
                epsilon=effective.get("epsilon", 1.0),
            )
            for point in flat
        ]
        shape = R.shape
        nstates = self.nstates
        H = np.stack([item["H_adi"] for item in values]).reshape(*shape, nstates, nstates)
        dH = np.stack([item["dH_adi"] for item in values]).reshape(*shape, nstates, nstates)
        dc = np.stack([item["DC1_adi"] for item in values]).reshape(*shape, nstates, nstates)
        eye = np.broadcast_to(np.eye(nstates), (*shape, nstates, nstates)).copy()
        return {
            "H_adi": H.astype(complex),
            "dH_adi": dH[..., None, :, :].astype(complex),
            "DC1_adi": dc[..., None, :, :].astype(complex),
            "basis_transform": eye.astype(complex),
            "time_overlap_adi": eye.astype(complex),
        }

    def __call__(self, R, P=None, storage=None, traj=None, params=None):
        arr = np.asarray(R, dtype=float)
        if arr.ndim == 1:
            arr = arr.reshape(1, self.ndof)
        return self.evaluate(np.swapaxes(arr, 0, 1), params=params)


def polariton_info(R, dvr_data, model="4-state", g_c=0.005, omega_c=0.1, epsilon=1.0):
    """
    Build Shin-Metiu polaritonic Hamiltonian data at one nuclear coordinate.

    The returned dictionary contains ``H_adi``, ``dH_adi``, and ``DC1_adi`` in
    the selected electronic/Fock product basis. The helper mirrors the legacy
    ``polariton_info`` routine but returns named NumPy arrays and keeps the
    derivative-coupling convention explicit.
    """

    values = dvr_data.interpolate(float(R))
    eigvals = values["eigvals"]
    d_V = values["d_V"]
    nac = values["nac"]
    mu = values["mu"]
    mu_deri = values["mu_deri"]
    H_el = np.diag(eigvals)
    H_deri = d_V + H_el @ nac - nac @ H_el
    dse = dipole_self_energy(mu, mu_deri, g_c=g_c, omega_c=omega_c, epsilon=epsilon)

    if model == "2-state":
        H, derivative, dc = build_two_state_polariton_hamiltonian(
            eigvals,
            H_deri,
            mu,
            mu_deri,
            dse["D_square"],
            dse["D_square_deri"],
            g_c=g_c,
            omega_c=omega_c,
            epsilon=epsilon,
        )
    elif model == "4-state":
        H, derivative, dc = build_four_state_polariton_hamiltonian(
            eigvals,
            H_deri,
            nac,
            mu,
            mu_deri,
            dse["D_square"],
            dse["D_square_deri"],
            g_c=g_c,
            omega_c=omega_c,
            epsilon=epsilon,
        )
    else:
        raise ValueError("model must be '2-state' or '4-state'")
    return {"H_adi": H, "dH_adi": derivative, "DC1_adi": dc}


def dipole_self_energy(mu, mu_deri, g_c=0.005, omega_c=0.1, epsilon=1.0):
    """
    Return dipole self-energy matrices for the Shin-Metiu polariton model.

    ``D_square = epsilon^2 g_c^2 mu^2 / omega_c`` and
    ``D_square_deri = epsilon^2 g_c^2 (mu mu_deri + mu_deri mu) / omega_c``.
    The derivative uses ``mu_deri``, i.e. the derivative of the dipole matrix
    elements in the electronic adiabatic basis.
    """

    prefactor = epsilon**2 * g_c**2 / omega_c
    return {
        "D_square": prefactor * (mu @ mu),
        "D_square_deri": prefactor * (mu @ mu_deri + mu_deri @ mu),
    }


def build_two_state_polariton_hamiltonian(
    eigvals,
    H_deri,
    mu,
    mu_deri,
    D_square,
    D_square_deri,
    g_c=0.005,
    omega_c=0.1,
    epsilon=1.0,
):
    """
    Build the ``|e,0>``, ``|g,1>`` Shin-Metiu polariton block.

    The matrix elements are
    ``H00 = E_e + D_ee``, ``H11 = E_g + D_gg + omega_c``, and
    ``H01 = g_c epsilon mu_eg``. The legacy implementation leaves the
    two-state derivative-coupling matrix at zero because this reduced block is
    a resonance subspace, not a full photon-block product basis.
    """

    H = np.zeros((2, 2), dtype=float)
    H[0, 0] = eigvals[1] + D_square[1, 1]
    H[0, 1] = g_c * epsilon * mu[1, 0]
    H[1, 0] = g_c * epsilon * mu[0, 1]
    H[1, 1] = eigvals[0] + D_square[0, 0] + omega_c

    derivative = np.zeros((2, 2), dtype=float)
    derivative[0, 0] = H_deri[1, 1] + D_square_deri[1, 1]
    derivative[0, 1] = g_c * epsilon * mu_deri[1, 0]
    derivative[1, 0] = g_c * epsilon * mu_deri[0, 1]
    derivative[1, 1] = H_deri[0, 0] + D_square_deri[0, 0]
    dc = np.zeros((2, 2), dtype=float)
    return H, derivative, dc


def build_four_state_polariton_hamiltonian(
    eigvals,
    H_deri,
    nac,
    mu,
    mu_deri,
    D_square,
    D_square_deri,
    g_c=0.005,
    omega_c=0.1,
    epsilon=1.0,
):
    """
    Build the ``|g,0>``, ``|e,0>``, ``|g,1>``, ``|e,1>`` polariton block.

    The diagonal contains electronic energies, dipole self-energy, and photon
    energies ``0.5 omega_c`` or ``1.5 omega_c``. Off-diagonal one-photon
    couplings are ``g_c epsilon mu_ij`` and intra-photon-block couplings are
    the off-diagonal dipole self-energy elements. The derivative returned is
    the covariant derivative ``dH = dV - H d + d H`` used in the original
    polariton source.
    """

    H = np.zeros((4, 4), dtype=float)
    H[0, 0] = eigvals[0] + D_square[0, 0] + 0.5 * omega_c
    H[1, 1] = eigvals[1] + D_square[1, 1] + 0.5 * omega_c
    H[2, 2] = eigvals[0] + D_square[0, 0] + 1.5 * omega_c
    H[3, 3] = eigvals[1] + D_square[1, 1] + 1.5 * omega_c
    H[0, 1] = D_square[0, 1]
    H[2, 3] = D_square[0, 1]
    H[0, 2] = g_c * epsilon * mu[0, 0]
    H[0, 3] = g_c * epsilon * mu[0, 1]
    H[1, 2] = g_c * epsilon * mu[1, 0]
    H[1, 3] = g_c * epsilon * mu[1, 1]
    _symmetrize_from_upper(H)

    dc = np.zeros((4, 4), dtype=float)
    dc[0, 1] = nac[0, 1]
    dc[2, 3] = nac[0, 1]
    dc[1, 0] = nac[1, 0]
    dc[3, 2] = nac[1, 0]

    derivative_plain = np.zeros((4, 4), dtype=float)
    derivative_plain[0, 0] = H_deri[0, 0] + D_square_deri[0, 0]
    derivative_plain[1, 1] = H_deri[1, 1] + D_square_deri[1, 1]
    derivative_plain[2, 2] = H_deri[0, 0] + D_square_deri[0, 0]
    derivative_plain[3, 3] = H_deri[1, 1] + D_square_deri[1, 1]
    derivative_plain[0, 1] = D_square_deri[0, 1]
    derivative_plain[2, 3] = D_square_deri[0, 1]
    derivative_plain[0, 2] = g_c * epsilon * mu_deri[0, 0]
    derivative_plain[0, 3] = g_c * epsilon * mu_deri[0, 1]
    derivative_plain[1, 2] = g_c * epsilon * mu_deri[1, 0]
    derivative_plain[1, 3] = g_c * epsilon * mu_deri[1, 1]
    _symmetrize_from_upper(derivative_plain)
    derivative = derivative_plain - H @ dc + dc @ H
    return H, derivative, dc


def _symmetrize_from_upper(matrix):
    i_upper, j_upper = np.triu_indices(matrix.shape[0], k=1)
    matrix[j_upper, i_upper] = matrix[i_upper, j_upper]
