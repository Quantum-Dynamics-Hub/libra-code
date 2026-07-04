from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.special import erf


@dataclass(frozen=True)
class ShinMetiuDVRData:
    """
    Interpolatable Shin-Metiu electronic DVR data.

    ``R_grid`` stores the nuclear coordinate. The remaining arrays are sampled
    on that grid and use electronic adiabatic-state indices as their last two
    axes. ``eigvals`` has shape ``(n_R, nstates)``; ``d_V``, ``nac``, ``mu``,
    ``d_mu``, and ``mu_deri`` have shape ``(n_R, nstates, nstates)``.

    The notation follows the legacy ``libra_py.models.Shin_Metiu`` sources:
    ``d_V`` is the adiabatic matrix of ``<phi_i|dV/dR|phi_j>``, ``nac`` is
    ``<phi_i|d/dR|phi_j>``, ``mu`` is the dipole matrix, ``d_mu`` is
    ``<phi_i|dmu/dR|phi_j>``, and ``mu_deri`` is the derivative of the dipole
    matrix elements including the basis-derivative contribution.

    Reference data source: ``libra_py.models.Shin_Metiu`` DVR files used for
    the Shin-Metiu polariton models of J. Chem. Phys. 157, 104118 (2022).
    """

    R_grid: np.ndarray
    eigvals: np.ndarray
    d_V: np.ndarray
    nac: np.ndarray
    mu: np.ndarray
    d_mu: np.ndarray
    mu_deri: np.ndarray

    @property
    def nstates(self) -> int:
        return int(self.eigvals.shape[-1])

    def interpolate(self, R):
        """
        Linearly interpolate all DVR quantities at one or more nuclear points.

        ``R`` may be scalar or array-like. The returned dictionary uses a
        leading batch shape matching ``np.asarray(R).shape``. ``dH_adi`` is
        reconstructed from the covariant derivative relation
        ``dH = d_V + H d - d H`` with ``H = diag(eigvals)`` and
        ``d = nac``.
        """

        points = np.asarray(R, dtype=float)
        flat = points.reshape(-1)
        nstates = self.nstates
        eigvals = np.empty((flat.size, nstates), dtype=float)
        matrices = {}
        for name in ("d_V", "nac", "mu", "d_mu", "mu_deri"):
            matrices[name] = np.empty((flat.size, nstates, nstates), dtype=float)

        for i in range(nstates):
            eigvals[:, i] = np.interp(flat, self.R_grid, self.eigvals[:, i])
            for j in range(nstates):
                for name in matrices:
                    source = getattr(self, name)
                    matrices[name][:, i, j] = np.interp(flat, self.R_grid, source[:, i, j])

        h = _diag_batch(eigvals)
        d_h = matrices["d_V"] + np.einsum("...ij,...jk->...ik", h, matrices["nac"])
        d_h -= np.einsum("...ij,...jk->...ik", matrices["nac"], h)

        out_shape = points.shape
        result = {
            "eigvals": eigvals.reshape(*out_shape, nstates),
            "H_adi": h.reshape(*out_shape, nstates, nstates),
            "dH_adi": d_h.reshape(*out_shape, nstates, nstates),
        }
        for name, value in matrices.items():
            result[name] = value.reshape(*out_shape, nstates, nstates)
        return result


@dataclass
class ShinMetiuDVRModel:
    """
    Cavity-free Shin-Metiu electronic DVR model.

    The Shin-Metiu model is a one-electron, one-nuclear-coordinate problem in
    which the electronic coordinate ``r`` is represented on a sine DVR grid and
    the nuclear coordinate is ``R = Q[0]``. The electronic Hamiltonian is
    ``H_el(R) = T_r + V_en(r; R)``. ``T_r`` is the Colbert-Miller/sine-DVR
    kinetic-energy matrix and ``V_en`` is a softened one-dimensional Coulomb
    potential. The model has two variants:

    ``model=1``
        Symmetric soft-Coulomb Shin-Metiu potential with parameters ``L_sm1``
        and ``Rc_sm1``.

    ``model=2``
        Variant with side-dependent smoothing lengths ``a_sm2[+1]`` and
        ``a_sm2[-1]`` plus central smoothing ``af_sm2``.

    This dyn translation returns the native adiabatic electronic quantities:
    ``H_adi = diag(E_i(R))``, ``dH_adi[..., 0, :, :]``, ``DC1_adi[..., 0, :, :]``,
    identity ``basis_transform``, and identity ``time_overlap_adi``. The
    derivative coupling convention is ``d_ij = <phi_i|d/dR|phi_j>`` with
    ``d_ij = <phi_i|dV/dR|phi_j> / (E_j - E_i)`` for ``i != j``.

    A precomputed ``ShinMetiuDVRData`` object or HDF5 file can be supplied. If
    neither is supplied, the DVR Hamiltonian is built and diagonalized on the
    fly from ``N``, ``r_min``, ``r_max``, and the model parameters.

    Legacy source: ``libra_py.models.Shin_Metiu.shin_metiu_dvr``.
    References: Shin and Metiu, J. Chem. Phys. 102, 9285 (1995); polaritonic
    DVR data usage in J. Chem. Phys. 157, 104118 (2022).
    """

    params: dict | None = None
    dvr_data: ShinMetiuDVRData | None = None
    dvr_file: str | Path | None = None
    nstates: int = 2
    ndof: int = 1

    def __post_init__(self):
        self.params = dict(self.params or {})
        if self.dvr_data is None and self.dvr_file is not None:
            self.dvr_data = load_shin_metiu_dvr_data(self.dvr_file)
        if self.dvr_data is not None:
            self.nstates = min(self.nstates, self.dvr_data.nstates)

    def evaluate(self, Q, params=None):
        effective = dict(self.params)
        if params:
            effective.update(params)
        R = np.asarray(Q, dtype=float)[0]
        if self.dvr_data is not None:
            values = self.dvr_data.interpolate(R)
            return _adiabatic_result_from_interpolated(values, self.nstates)
        return _compute_on_the_fly_result(R, effective, self.nstates)

    def __call__(self, R, P=None, storage=None, traj=None, params=None):
        arr = np.asarray(R, dtype=float)
        if arr.ndim == 1:
            arr = arr.reshape(1, self.ndof)
        return self.evaluate(np.swapaxes(arr, 0, 1), params=params)


def bundled_shin_metiu_dvr_path(model: int) -> Path:
    """Return the legacy bundled HDF5 path for Shin-Metiu model ``1`` or ``2``."""

    if model not in (1, 2):
        raise ValueError("model must be 1 or 2")
    root = Path(__file__).resolve().parents[2]
    return root / "models" / "Shin_Metiu" / f"dvr_sm{model}.h5"


def load_shin_metiu_dvr_data(filename: str | Path) -> ShinMetiuDVRData:
    """Load a Shin-Metiu DVR HDF5 file into an immutable data object."""

    import h5py

    with h5py.File(filename, "r") as handle:
        return ShinMetiuDVRData(
            R_grid=handle["R_grid"][:],
            eigvals=handle["eigvals"][:],
            d_V=handle["d_V"][:],
            nac=handle["nac"][:],
            mu=handle["mu"][:],
            d_mu=handle["d_mu"][:],
            mu_deri=handle["mu_deri"][:],
        )


def kinetic_energy_matrix(npoints: int, spacing: float) -> np.ndarray:
    """
    Construct the one-dimensional sine-DVR kinetic energy matrix.

    The legacy Shin-Metiu implementation uses atomic units with electronic mass
    set to one, so the prefactor is ``1/2``. ``spacing`` is the uniform DVR grid
    spacing in the electronic coordinate ``r``.
    """

    i = np.arange(npoints)[:, None]
    j = np.arange(npoints)[None, :]
    delta = j - i
    matrix = np.empty((npoints, npoints), dtype=float)
    diagonal = delta == 0
    matrix[diagonal] = 0.5 * np.pi**2 / (3.0 * spacing**2) * (1.0 + 2.0 / npoints**2)
    off = ~diagonal
    matrix[off] = (
        0.5
        * 2.0
        * (-1.0) ** delta[off]
        * np.pi**2
        / (spacing * npoints * np.sin(np.pi * delta[off] / npoints)) ** 2
    )
    return matrix


def shin_metiu_potential_sm1(r_grid, R, L_sm1=19.0, Rc_sm1=5.0):
    """Return the model-1 softened electron-nuclear potential on ``r_grid``."""

    r_grid = np.asarray(r_grid)
    value = np.zeros_like(r_grid, dtype=float)
    for sigma in (1.0, -1.0):
        R_shift = R + sigma * L_sm1 / 2.0
        r_shift = r_grid + sigma * L_sm1 / 2.0
        value += 1.0 / np.abs(R_shift) - erf(np.abs(r_shift) / Rc_sm1) / np.abs(r_shift)
    value -= erf(np.abs(R - r_grid) / Rc_sm1) / np.abs(R - r_grid)
    return value


def shin_metiu_potential_derivative_sm1(r_grid, R, L_sm1=19.0, Rc_sm1=5.0):
    """Return ``dV_en/dR`` for Shin-Metiu model 1 on ``r_grid``."""

    r_grid = np.asarray(r_grid)
    value = np.zeros_like(r_grid, dtype=float)
    for sigma in (1.0, -1.0):
        R_shift = R + sigma * L_sm1 / 2.0
        value += -np.sign(R_shift) / np.abs(R_shift) ** 2
    diff = R - r_grid
    value -= 2.0 / np.sqrt(np.pi) / Rc_sm1 * np.exp(-(diff / Rc_sm1) ** 2) * np.sign(diff) / np.abs(diff)
    value += erf(np.abs(diff) / Rc_sm1) * np.sign(diff) / np.abs(diff) ** 2
    return value


def shin_metiu_potential_sm2(r_grid, R, L_sm2=10.0, a_sm2=None, af_sm2=5.0):
    """Return the model-2 softened electron-nuclear potential on ``r_grid``."""

    a_sm2 = {1: 4.0, -1: 3.1} if a_sm2 is None else a_sm2
    r_grid = np.asarray(r_grid)
    value = np.zeros_like(r_grid, dtype=float)
    for sigma in (1, -1):
        R_shift = R + sigma * L_sm2 / 2.0
        r_shift = r_grid + sigma * L_sm2 / 2.0
        value += 1.0 / np.abs(R_shift) - erf(np.abs(r_shift) / a_sm2[sigma]) / np.abs(r_shift)
    value -= erf(np.abs(R - r_grid) / af_sm2) / np.abs(R - r_grid)
    return value


def shin_metiu_potential_derivative_sm2(r_grid, R, L_sm2=10.0, a_sm2=None, af_sm2=5.0):
    """Return ``dV_en/dR`` for Shin-Metiu model 2 on ``r_grid``."""

    del a_sm2
    r_grid = np.asarray(r_grid)
    value = np.zeros_like(r_grid, dtype=float)
    for sigma in (1.0, -1.0):
        R_shift = R + sigma * L_sm2 / 2.0
        value += -np.sign(R_shift) / np.abs(R_shift) ** 2
    diff = R - r_grid
    value -= 2.0 / np.sqrt(np.pi) / af_sm2 * np.exp(-(diff / af_sm2) ** 2) * np.sign(diff) / np.abs(diff)
    value += erf(np.abs(diff) / af_sm2) * np.sign(diff) / np.abs(diff) ** 2
    return value


def dipole_matrix_elements(r_grid, R):
    """Return the position-basis dipole diagonal ``mu(r; R) = R - r``."""

    return R - np.asarray(r_grid)


def _compute_on_the_fly_result(R, params, nstates):
    flat = np.asarray(R, dtype=float).reshape(-1)
    shape = np.asarray(R).shape
    values = [_single_point_dvr(point, params, nstates) for point in flat]
    H = np.stack([item["H_adi"] for item in values]).reshape(*shape, nstates, nstates)
    dH = np.stack([item["dH_adi"] for item in values]).reshape(*shape, nstates, nstates)
    dc = np.stack([item["DC1_adi"] for item in values]).reshape(*shape, nstates, nstates)
    eye = np.broadcast_to(np.eye(nstates), (*shape, nstates, nstates)).copy()
    return {
        "H_adi": H.astype(complex),
        "dH_adi": dH[:, None, :, :].reshape(*shape, 1, nstates, nstates).astype(complex),
        "DC1_adi": dc[:, None, :, :].reshape(*shape, 1, nstates, nstates).astype(complex),
        "basis_transform": eye.astype(complex),
        "time_overlap_adi": eye.astype(complex),
    }


def _single_point_dvr(R, params, nstates):
    npoints = params.get("N", 201)
    r_min = params.get("r_min", -20.0)
    r_max = params.get("r_max", 20.0)
    model = params.get("model", 1)
    r_grid = np.linspace(r_min, r_max, npoints)
    spacing = (r_max - r_min) / (npoints - 1)

    if model == 1:
        potential = shin_metiu_potential_sm1(
            r_grid,
            R,
            params.get("L_sm1", 19.0),
            params.get("Rc_sm1", 5.0),
        )
        derivative = shin_metiu_potential_derivative_sm1(
            r_grid,
            R,
            params.get("L_sm1", 19.0),
            params.get("Rc_sm1", 5.0),
        )
    elif model == 2:
        potential = shin_metiu_potential_sm2(
            r_grid,
            R,
            params.get("L_sm2", 10.0),
            params.get("a_sm2", {1: 4.0, -1: 3.1}),
            params.get("af_sm2", 5.0),
        )
        derivative = shin_metiu_potential_derivative_sm2(
            r_grid,
            R,
            params.get("L_sm2", 10.0),
            params.get("a_sm2", {1: 4.0, -1: 3.1}),
            params.get("af_sm2", 5.0),
        )
    else:
        raise ValueError("model must be 1 or 2")

    h_el = kinetic_energy_matrix(npoints, spacing) + np.diag(potential)
    eigvals, eigvecs = np.linalg.eigh(h_el)
    eigvals = eigvals[:nstates]
    eigvecs = eigvecs[:, :nstates]
    d_v = eigvecs.T @ np.diag(derivative) @ eigvecs
    dc = np.zeros_like(d_v)
    for i in range(nstates):
        for j in range(nstates):
            if i != j:
                dc[i, j] = d_v[i, j] / (eigvals[j] - eigvals[i])
    H = np.diag(eigvals)
    dH = d_v + H @ dc - dc @ H
    return {"H_adi": H, "dH_adi": dH, "DC1_adi": dc}


def _adiabatic_result_from_interpolated(values, nstates):
    H = values["H_adi"][..., :nstates, :nstates]
    dH = values["dH_adi"][..., :nstates, :nstates]
    dc = values["nac"][..., :nstates, :nstates]
    shape = H.shape[:-2]
    eye = np.broadcast_to(np.eye(nstates), (*shape, nstates, nstates)).copy()
    return {
        "H_adi": H.astype(complex),
        "dH_adi": dH[..., None, :, :].astype(complex),
        "DC1_adi": dc[..., None, :, :].astype(complex),
        "basis_transform": eye.astype(complex),
        "time_overlap_adi": eye.astype(complex),
    }


def _diag_batch(values):
    out = np.zeros((*values.shape[:-1], values.shape[-1], values.shape[-1]), dtype=values.dtype)
    idx = np.arange(values.shape[-1])
    out[..., idx, idx] = values
    return out
