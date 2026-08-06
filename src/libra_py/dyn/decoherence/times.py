"""Decoherence-rate and coherence-time models translated from C++.

The source is ``src/dyn/dyn_decoherence_time.cpp``.  Arrays use the new
Python convention: trajectory is the first axis, so amplitudes are
``(ntraj, nstates)`` and rates are ``(ntraj, nstates, nstates)``.
"""

from __future__ import annotations

import numpy as np

_KB_HARTREE_PER_K = 3.166811563e-6


def edc_rates(hvib, kinetic_energy, C_param=1.0, eps_param=0.1):
    r"""Return Granucci--Persico energy-based decoherence rates.

    .. math:: \tau_{ij}^{-1}=|E_i-E_j|/(C+\epsilon/E_{kin}).

    Parameters
    ----------
    hvib : array_like
        Vibronic Hamiltonian with shape ``(nstates,nstates)`` or
        ``(ntraj,nstates,nstates)``. Only real diagonal energies are used.
    kinetic_energy : float or array_like
        Nuclear kinetic energy in hartree, scalar or one value per trajectory.
    C_param, eps_param : float
        EDC empirical parameters in hartree; common values are 1.0 and 0.1.

    Returns
    -------
    numpy.ndarray
        Inverse decoherence times in atomic units of inverse time. The output
        has the same single/batched trajectory convention as ``hvib``.

    Notes
    -----
    ``hvib`` may be one matrix or a trajectory batch. Unlike the C++ code,
    zero kinetic energy is handled explicitly and gives zero rates instead of
    relying on floating-point division by zero.
    """
    h = np.asarray(hvib)
    single = h.ndim == 2
    h = h[None] if single else h
    ek = np.broadcast_to(np.asarray(kinetic_energy, dtype=float), (len(h),))
    energies = np.diagonal(h.real, axis1=-2, axis2=-1)
    gaps = np.abs(energies[:, :, None] - energies[:, None, :])
    denominator = np.full_like(ek, np.inf)
    np.divide(eps_param, ek, out=denominator, where=ek != 0.0)
    denominator += C_param
    result = gaps / denominator[:, None, None]
    return result[0] if single else result


def dephasing_informed_correction(rates, hvib, average_gaps):
    """Apply the Sifain--Wang--Teritiak--Prezhdo gap correction.

    Parameters are an uncorrected rate matrix/batch, the corresponding
    vibronic Hamiltonian matrix/batch, and ``average_gaps[i,j] =
    <|E_i-E_j|>``. Shapes follow :func:`edc_rates`.

    Returns a corrected copy; inputs are not modified. The correction is
    ``rate_ij *= |E_i-E_j| / average_gap_ij``. Where an average gap is
    non-positive, the
    C++ sentinel rate ``1e25`` is retained.
    """
    r = np.array(rates, dtype=float, copy=True)
    h = np.asarray(hvib)
    single = h.ndim == 2
    h = h[None] if single else h
    if r.ndim == 2:
        r = r[None]
    e = np.diagonal(h.real, axis1=-2, axis2=-1)
    gaps = np.abs(e[:, :, None] - e[:, None, :])
    avg = np.asarray(average_gaps, dtype=float)
    r = np.where(avg > 0.0, r * gaps / np.where(avg > 0.0, avg, 1.0), 1e25)
    return r[0] if single else r


def coherence_intervals(amplitudes, rates):
    r"""Compute population-dependent DISH coherence intervals.

    .. math:: \tau_i^{-1}=\sum_{j\ne i}|c_j|^2\,\tau_{ij}^{-1}.

    This equation follows the C++ implementation.  Its old doc comment says
    ``rho_ii`` inside the sum, but the executable C++ code correctly uses the
    population of state ``j``; this Python documentation makes that historical
    C++ documentation bug explicit.

    Parameters
    ----------
    amplitudes : array_like
        Complex coefficients, ``(nstates,)`` or ``(ntraj,nstates)``.
    rates : array_like
        Pair rates, ``(nstates,nstates)`` or trajectory batch.

    Returns
    -------
    numpy.ndarray
        State-resolved intervals, ``(nstates,)`` or ``(ntraj,nstates)``.
        An inverse rate of zero is represented by the C++ sentinel ``1e25``.
    """
    c = np.asarray(amplitudes)
    single = c.ndim == 1
    c = c[None] if single else c
    r = np.asarray(rates, dtype=float)
    r = np.broadcast_to(r if r.ndim == 3 else r[None], (len(c),) + r.shape[-2:])
    pop = np.abs(c) ** 2
    inverse = np.einsum("tij,tj->ti", r, pop) - np.diagonal(r, axis1=1, axis2=2) * pop
    result = np.full_like(inverse, 1e25)
    np.divide(1.0, inverse, out=result, where=inverse > 0.0)
    return result[0] if single else result


def schwartz_1(amplitudes, state_forces, inverse_alpha, mean_field_forces=None):
    """Return diagonal Schwartz-1 rates from mean-field/state force differences.

    ``amplitudes`` is ``(ntraj,nstates)``, ``state_forces`` is
    ``(ntraj,nstates,ndof)``, and ``inverse_alpha`` supplies one Gaussian
    inverse-width parameter per nuclear DOF. If ``mean_field_forces`` with
    shape ``(ntraj,ndof)`` is omitted, it is built from electronic
    populations, as in the C++ Hamiltonian call.

    Returns a trajectory batch of ``(nstates,nstates)`` matrices. Schwartz-1
    rates occupy only the diagonal because they are state-resolved rather
    than pair-resolved.
    """
    c = np.asarray(amplitudes)
    f = np.asarray(state_forces, dtype=float)
    inv = np.asarray(inverse_alpha, dtype=float).reshape(-1)
    mf = (np.einsum("ts,tsd->td", np.abs(c) ** 2, f) if mean_field_forces is None
          else np.asarray(mean_field_forces, dtype=float))
    diag = np.sqrt(0.25 * np.einsum("tsd,d,tsd->ts", mf[:, None] - f, inv, mf[:, None] - f))
    result = np.zeros(f.shape[:2] + (f.shape[1],), dtype=float)
    i = np.arange(f.shape[1]); result[:, i, i] = diag
    return result


def schwartz_1_interaction_width(amplitudes, state_forces, momentum, frequencies):
    """Compute the interaction-width Schwartz-1 rates (C++ option 4).

    ``momentum`` has shape ``(ntraj,ndof)`` and ``frequencies`` has one value
    per DOF. The effective inverse width is
    ``(4*pi/(frequency**2 * momentum))**2``. Zero momentum therefore produces
    an infinite rate, consistently with the algebra in the C++ source.

    Returns diagonal rate matrices shaped ``(ntraj,nstates,nstates)``.
    """
    p = np.asarray(momentum, dtype=float)
    w = np.asarray(frequencies, dtype=float).reshape(-1)
    with np.errstate(divide="ignore", invalid="ignore"):
        inv_alpha = (4.0 * np.pi / (w * w * p)) ** 2
    c = np.asarray(amplitudes); f = np.asarray(state_forces, dtype=float)
    mf = np.einsum("ts,tsd->td", np.abs(c) ** 2, f)
    diag = np.sqrt(0.25 * np.sum(inv_alpha[:, None, :] * (mf[:, None] - f) ** 2, axis=2))
    out = np.zeros((len(c), f.shape[1], f.shape[1])); i=np.arange(f.shape[1]); out[:,i,i]=diag
    return out


def schwartz_2(state_forces, inverse_alpha):
    r"""Return pairwise Schwartz-2 force-difference rates.

    ``state_forces`` is ``(ntraj,nstates,ndof)`` and ``inverse_alpha`` contains
    one Gaussian inverse width per DOF. The symmetric result satisfies
    :math:`\tau_{ij}^{-1}=[\frac14\sum_k\alpha_k^{-1}
    (F_{ik}-F_{jk})^2]^{1/2}` and has shape
    ``(ntraj,nstates,nstates)``.
    """
    f = np.asarray(state_forces, dtype=float)
    inv = np.asarray(inverse_alpha, dtype=float).reshape(-1)
    df = f[:, :, None, :] - f[:, None, :, :]
    return np.sqrt(0.25 * np.einsum("tijd,d,tijd->tij", df, inv, df))


def gu_franco(amplitudes, reorganization_energy, temperature):
    r"""Return Gu--Franco rates including the published 2020 correction.

    The rate is :math:`|\rho_{ij}|\sqrt{4\lambda k_B T}`. ``amplitudes`` may
    be ``(nstates,)`` or ``(ntraj,nstates)``; reorganization energy is in
    hartree and temperature in kelvin. The returned matrix/batch is expressed
    in atomic units of inverse time.
    """
    c = np.asarray(amplitudes)
    single = c.ndim == 1
    c = c[None] if single else c
    dm_abs = np.abs(c[:, :, None] * c[:, None, :].conj())
    result = dm_abs * np.sqrt(4.0 * reorganization_energy * _KB_HARTREE_PER_K * temperature)
    return result[0] if single else result
