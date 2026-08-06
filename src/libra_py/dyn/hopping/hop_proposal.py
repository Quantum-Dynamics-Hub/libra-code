"""Surface-hop proposal algorithms.

This module carries the algorithms and mathematical notes from
``src/dyn/dyn_hop_proposal.cpp`` into the Python dynamics implementation.  A
probability matrix uses the convention ``g[i, j] = P(i → j)``; a probability
vector contains all outcomes from its supplied active state.  Staying
probabilities are stored in ``g[i]`` or ``g[i, i]`` so that a proposal row sums
to one, subject to the same clipping rules as the corresponding C++ routine.

All energies and times are in atomic units.  Temperature is in kelvin, with
``k_B = 3.166811429e-6 Hartree/K``.  When ``use_boltz_factor`` is enabled, an
uphill proposal ``i → j`` is multiplied by

    exp[-(H_jj - H_ii)/(k_B T)],   H_jj > H_ii.

Downhill proposals are unchanged.  This scaling is intended primarily for
NBRA calculations where hop proposal and acceptance are combined.

The functions accept NumPy-compatible arrays.  As in the C++ overloads, FSSH,
GFSH, and MSSH return a complete transition matrix when passed an amplitude
vector without an active state, or one probability vector when passed a
density matrix and an active-state index.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from .hop_proposal_fssh3 import hopping_probabilities_fssh3


KB_HARTREE_PER_K = 3.166811429e-6

# Keep these numeric identifiers synchronized with dyn_hop_proposal.cpp and
# DynControlParams.tsh_method.  Do not renumber when adding implementations.
TSH_METHODS = {
    -1: "adiabatic/no hops",
    0: "FSSH",
    1: "GFSH",
    2: "MSSH",
    3: "Landau-Zener",
    4: "Zhu-Nakamura",
    5: "DISH (handled by the decoherence event scheduler)",
    6: "MASH",
    7: "FSSH2",
    8: "FSSH3",
    9: "original GFSH",
}


def hop_proposal_probabilities(
    params,
    density,
    hvib,
    active_states,
    density_old=None,
    *,
    ham=None,
    ham_prev=None,
    momentum=None,
    inverse_mass=None,
    fssh3_errors=None,
):
    """Compute proposal vectors for one trajectory or a batch.

    ``density`` and ``hvib`` may have shape ``(nstates, nstates)`` or
    ``(ntraj, nstates, nstates)``.  ``active_states`` is a scalar for the
    single-trajectory form and a vector for the batch form.  Methods based on
    finite population differences additionally require ``density_old``.

    LZ (method 3) and ZN (method 4) use the keyword-only ``ham``, ``ham_prev``,
    ``momentum``, and ``inverse_mass`` inputs.  Their representation is read
    from ``params.rep_lz``. For a batch, ``ham`` and ``ham_prev`` are sequences
    of Hamiltonian records and ``momentum`` has shape ``(ntraj, ndof)`` or the
    C++-style ``(ndof, ntraj)``. ``inverse_mass`` is shared by trajectories.

    The ``tsh_method`` mapping follows :class:`DynControlParams`: 0 FSSH,
    1 GFSH, 2 MSSH, 3 LZ, 4 ZN, 6 MASH, 7 FSSH2, 8 FSSH3, and 9 original
    GFSH. Method -1 returns stay-on-state probabilities. Method 5 is an
    explicit placeholder because, as noted in C++, DISH proposals are produced
    by the decoherence event scheduler rather than this probability interface.
    """

    method = int(_param(params, "tsh_method", -1))
    if method not in TSH_METHODS:
        choices = ", ".join(str(value) for value in TSH_METHODS)
        raise ValueError(f"tsh_method must be one of {choices}; got {method}")
    if method in (3, 4):
        missing = [
            name
            for name, value in (
                ("ham", ham),
                ("ham_prev", ham_prev),
                ("momentum", momentum),
                ("inverse_mass", inverse_mass),
            )
            if value is None
        ]
        if missing:
            raise ValueError(f"tsh_method {method} requires {', '.join(missing)}")
        return _history_probabilities(
            hopping_probabilities_lz if method == 3 else hopping_probabilities_zn,
            ham,
            ham_prev,
            active_states,
            int(_param(params, "rep_lz", 0)),
            momentum,
            inverse_mass,
        )

    if method == 5:
        raise NotImplementedError(
            "tsh_method 5 (DISH) is handled by the decoherence event scheduler; "
            "its event-time proposal interface is not implemented yet"
        )

    densities = np.asarray(density, dtype=complex)
    hamiltonians = np.asarray(hvib, dtype=complex)
    states = np.asarray(active_states, dtype=int)
    single = densities.ndim == 2
    if single:
        densities = densities[None, ...]
        hamiltonians = hamiltonians[None, ...]
        states = states.reshape(1)
    if densities.ndim != 3 or hamiltonians.ndim != 3:
        raise ValueError("density and hvib must be matrices or batches of matrices")
    if len(densities) != len(hamiltonians) or len(densities) != len(states):
        raise ValueError("density, hvib, and active_states batch sizes must match")

    old_batch = None
    if density_old is not None:
        old_batch = np.asarray(density_old, dtype=complex)
        if old_batch.ndim == 2:
            old_batch = old_batch[None, ...]
        if len(old_batch) != len(densities):
            raise ValueError("density_old batch size must match density")

    result = []
    for index, (denmat, ham, state) in enumerate(zip(densities, hamiltonians, states)):
        if method == -1:
            probabilities = np.zeros(denmat.shape[0], dtype=float)
            probabilities[_state_index(state, denmat.shape[0])] = 1.0
        elif method == 0:
            probabilities = hopping_probabilities_fssh(params, denmat, ham, state)
        elif method == 1:
            probabilities = hopping_probabilities_gfsh(params, denmat, ham, state)
        elif method == 2:
            probabilities = hopping_probabilities_mssh(params, denmat, ham, state)
        elif method == 6:
            probabilities = hopping_probabilities_mash(params, denmat)
        elif method in (7, 8, 9):
            if old_batch is None:
                raise ValueError(f"tsh_method {method} requires density_old")
            if method == 7:
                probabilities = hopping_probabilities_fssh2(
                    params, denmat, old_batch[index], state
                )
            elif method == 8:
                errors = None
                if fssh3_errors is not None:
                    error_batch = np.asarray(fssh3_errors)
                    errors = error_batch if error_batch.ndim == 1 else error_batch[index]
                probabilities = hopping_probabilities_fssh3(
                    params, denmat, old_batch[index], state, errors
                )
            else:
                probabilities = hopping_probabilities_gfsh_orig(
                    params, denmat, old_batch[index], state
                )
        result.append(probabilities)
    result = np.asarray(result)
    return result[0] if single else result


def hopping_probabilities_fssh(params, state, hvib, active_state=None):
    """Compute fewest-switches surface-hopping (FSSH) probabilities.

    Starting with ``ρ = c c†`` and ``dc/dt = -i H_vib c`` (ħ = 1),

        dρ/dt = i(ρ H_vib - H_vib ρ).

    The population flux assigned to a transition from active state ``i`` to
    state ``j`` gives Tully's proposal

        g_ij = (dt / ρ_ii) Im(ρ_ij H_ji - H_ij ρ_ji),   i != j.

    Negative fluxes are set to zero.  The density-matrix/active-state form
    caps each off-diagonal probability at one and caps their sum at one, as in
    the C++ vector overload.  A negligible active population
    (``ρ_ii < 1e-8``) produces no outgoing proposals.

    Parameters
    ----------
    params
        :class:`~libra_py.dyn.control_params.DynControlParams` or mapping.
        Uses ``dt``, ``Temperature``, and ``use_boltz_factor``.
    state
        Electronic amplitudes ``c`` when ``active_state`` is omitted, or the
        density matrix ``ρ`` otherwise.
    hvib
        Dynamically consistent vibronic Hamiltonian in the same representation
        as ``state``.
    active_state
        Index of the initial state, or ``None`` for the all-pairs form.

    References
    ----------
    J. C. Tully, J. Chem. Phys. 93, 1061–1071 (1990); E. Fabiano,
    T. W. Keal, and W. Thiel, Chem. Phys. 349, 334–347 (2008);
    A. V. Akimov, J. Comput. Chem. 37, 1626–1649 (2016).
    """

    hvib = _square_matrix(hvib, "hvib")
    if active_state is None:
        coeff = _amplitudes(state, hvib.shape[0])
        density = np.outer(coeff, coeff.conj())
        return np.vstack(
            [_fssh_vector(params, density, hvib, i, cap=False) for i in range(len(coeff))]
        )
    density = _density_matrix(state, hvib.shape[0])
    return _fssh_vector(params, density, hvib, active_state, cap=True)


def _fssh_vector(params, density, hvib, active_state, *, cap):
    i = _state_index(active_state, density.shape[0])
    result = np.zeros(density.shape[0], dtype=float)
    population = density[i, i].real
    if population >= 1.0e-8:
        for j in range(density.shape[0]):
            if j == i:
                continue
            flux = (density[i, j] * hvib[j, i] - hvib[i, j] * density[j, i]).imag
            probability = _param(params, "dt", 41.0) * flux / population
            probability *= _boltzmann_factor(params, hvib[j, j].real - hvib[i, i].real)
            result[j] = max(0.0, min(1.0, probability) if cap else probability)
    result[i] = 1.0 - min(1.0, result.sum()) if cap else 1.0 - result.sum()
    return result


def hopping_probabilities_gfsh(params, state, hvib, active_state=None):
    """Compute global-flux surface-hopping (GFSH) probabilities.

    The instantaneous density derivative is

        ρ_dot = i(ρ H_vib† - H_vib ρ).

    Let ``a_k = Re(ρ_kk)``, ``a_dot_k = Re(ρ_dot_kk)``, and let
    ``N = sum(a_dot_k for a_dot_k < 0)`` be the total decaying flux.  GFSH
    assigns flux only from a decaying active state to growing states:

        g_ij = dt (a_dot_j / a_i) (a_dot_i / N),

    provided ``a_i > 1e-12``, ``a_dot_i < 0``, ``a_dot_j > 0``, and
    ``|N| > 1e-12``.  All other outgoing probabilities are zero.

    ``state`` is an amplitude vector when ``active_state`` is ``None`` and a
    density matrix otherwise.  ``params`` supplies ``dt``, ``Temperature``,
    and ``use_boltz_factor``.

    References
    ----------
    L. Wang, D. Trivedi, and O. V. Prezhdo, J. Chem. Theory Comput. 10,
    3598–3605 (2014); A. V. Akimov, J. Comput. Chem. 37, 1626–1649 (2016).
    """

    hvib = _square_matrix(hvib, "hvib")
    if active_state is None:
        coeff = _amplitudes(state, hvib.shape[0])
        density = np.outer(coeff, coeff.conj())
        return np.vstack(
            [_gfsh_vector(params, density, hvib, i) for i in range(len(coeff))]
        )
    return _gfsh_vector(params, _density_matrix(state, hvib.shape[0]), hvib, active_state)


def _gfsh_vector(params, density, hvib, active_state):
    i = _state_index(active_state, density.shape[0])
    populations = np.diag(density).real
    density_dot = 1j * (density @ hvib.conj().T - hvib @ density)
    rates = np.diag(density_dot).real
    total_decay = rates[rates < 0.0].sum()
    result = np.zeros(density.shape[0], dtype=float)

    if populations[i] >= 1.0e-12 and abs(total_decay) > 1.0e-12 and rates[i] < 0.0:
        for j in range(density.shape[0]):
            if j != i and rates[j] > 0.0:
                probability = (
                    _param(params, "dt", 41.0)
                    * (rates[j] / populations[i])
                    * rates[i]
                    / total_decay
                )
                probability *= _boltzmann_factor(
                    params, hvib[j, j].real - hvib[i, i].real
                )
                result[j] = max(0.0, probability)
    result[i] = max(0.0, 1.0 - result.sum())
    return result


def hopping_probabilities_gfsh_orig(params, density, density_old, active_state):
    """Compute the original finite-difference form of GFSH.

    This form replaces the instantaneous rates in GFSH by population changes
    ``Δa_k = Re[ρ_kk(t + dt) - ρ_kk(t)]``.  With
    ``N = sum(Δa_k for Δa_k < 0)``, transitions from a decaying active state
    ``i`` to a growing state ``j`` are

        g_ij = (Δa_j / a_i(t + dt)) (Δa_i / N).

    No explicit ``dt`` appears because it cancels between the finite-difference
    rates.  The same population and flux thresholds as GFSH are used.

    Reference: L. Wang, D. Trivedi, and O. V. Prezhdo,
    J. Chem. Theory Comput. 10, 3598–3605 (2014).
    """

    density = _square_matrix(density, "density")
    density_old = _density_matrix(density_old, density.shape[0])
    i = _state_index(active_state, density.shape[0])
    populations = np.diag(density).real
    changes = np.diag(density - density_old).real
    total_decay = changes[changes < 0.0].sum()
    result = np.zeros(density.shape[0], dtype=float)

    if populations[i] >= 1.0e-12 and abs(total_decay) > 1.0e-12 and changes[i] < 0.0:
        for j in range(density.shape[0]):
            if j != i and changes[j] > 0.0:
                result[j] = max(
                    0.0, (changes[j] / populations[i]) * changes[i] / total_decay
                )
    result[i] = max(0.0, 1.0 - result.sum())
    return result


def hopping_probabilities_fssh2(params, density, density_old, active_state):
    """Compute the population-difference FSSH2 proposal probabilities.

    FSSH2 avoids explicit nonadiabatic couplings.  For initial state ``m``,

        P_m,out = max[0, -(ρ_mm(t+dt) - ρ_mm(t)) / ρ_mm(t)],

        q_mn = max[0, (ρ_nn(t+dt) - ρ_nn(t)) / ρ_mm(t)].

    With ``fssh2_revision == 0`` (the original prescription),
    ``g_mn = min(q_mn, P_m,out)`` independently for every target.  With
    revision 1, the positive ``q_mn`` values are normalized so their total is
    ``1 - P_m,out`` and the stored diagonal value is ``P_m,out``, matching the
    revised C++ implementation.  Zero initial population implies no outgoing
    proposals.
    """

    density = _square_matrix(density, "density")
    density_old = _density_matrix(density_old, density.shape[0])
    i = _state_index(active_state, density.shape[0])
    old_populations = np.diag(density_old).real
    changes = np.diag(density - density_old).real
    result = np.zeros(density.shape[0], dtype=float)
    old_active = old_populations[i]
    out_probability = max(0.0, -changes[i] / old_active) if old_active > 0.0 else 0.0

    if old_active > 0.0:
        result = np.maximum(0.0, changes / old_active)
        result[i] = 0.0
    revision = int(_param(params, "fssh2_revision", 0))
    if revision == 0:
        result = np.minimum(result, out_probability)
        result[i] = 1.0 - result.sum()
    elif revision == 1:
        incoming = result.sum()
        if incoming > 0.0:
            result *= (1.0 - out_probability) / incoming
            result[i] = out_probability
        else:
            result[i] = 1.0
    else:
        raise ValueError("fssh2_revision must be 0 or 1")
    return result


def hopping_probabilities_mssh(params, state, hvib, active_state=None):
    """Compute Markov-state surface-hopping (MSSH) probabilities.

    MSSH proposes target states directly in proportion to their electronic
    populations.  For ``N = <Ψ|Ψ> = Tr(ρ)``,

        g_ij = ρ_jj / N,   j != i,
        g_ii = 1 - sum(j != i, g_ij).

    Uphill off-diagonal terms may be Boltzmann-scaled.  ``state`` is an
    amplitude vector when ``active_state`` is omitted and a density matrix
    otherwise.

    Reference: A. V. Akimov, D. Trivedi, L. Wang, and O. V. Prezhdo,
    J. Phys. Soc. Jpn. 84, 094002 (2015).
    """

    hvib = _square_matrix(hvib, "hvib")
    if active_state is None:
        coeff = _amplitudes(state, hvib.shape[0])
        norm = np.vdot(coeff, coeff).real
        if norm <= 0.0:
            raise ValueError("amplitudes must have nonzero norm")
        populations = np.abs(coeff) ** 2 / norm
        result = np.empty((len(coeff), len(coeff)), dtype=float)
        for i in range(len(coeff)):
            result[i] = _mssh_vector(params, populations, hvib, i)
        return result

    density = _density_matrix(state, hvib.shape[0])
    norm = np.trace(density).real
    if norm <= 0.0:
        raise ValueError("density matrix must have positive trace")
    return _mssh_vector(params, np.diag(density).real / norm, hvib, active_state)


def hopping_probabilities_lz(
    ham, ham_prev, active_state, rep, momentum, inverse_mass
):
    """Compute Landau–Zener (LZ) hop-proposal probabilities.

    This is the NumPy analogue of the C++ ``nHamiltonian`` overload. ``ham``
    and ``ham_prev`` may be mappings or objects exposing these arrays:

    - ``ham_dia`` in every representation;
    - ``d1ham_dia`` for ``rep == 0``;
    - ``ham_adi`` and ``nac_adi`` for ``rep in (1, 2)``;
    - previous ``ham_dia`` for ``rep == 1`` or previous ``nac_adi`` for
      ``rep == 2``.

    For a diabatic calculation (``rep == 0``), a proposal is considered only
    when the diabatic gap changes sign between the current and previous steps:

        ΔH_ij(t) ΔH_ij(t-dt) < 0,
        ΔH_ij = H_ii - H_jj.

    At such a crossing, the Belyaev–Lebedev/Tully expression is

        g_ij = exp[-2π H_ij² / Σ_k |v_k(dH_ii/dR_k-dH_jj/dR_k)|],
        v_k = p_k M_k⁻¹.

    For the adiabatic variants, ``rep == 1`` locates a crossing by the same
    diabatic-gap sign change, whereas ``rep == 2`` locates it by a sign change
    of the time-derivative NAC.  The probability is then

        g_ij = exp[-π |E_i-E_j| / (4 |NAC_ij|)].

    Zero denominators or NACs give zero transition probability.  ``momentum``
    and ``inverse_mass`` are one-dimensional nuclear arrays; as in C++, they
    are used only by the diabatic formula.

    References
    ----------
    J. C. Tully, J. Chem. Phys. 93, 1061–1071 (1990); A. K. Belyaev and
    O. V. Lebedev, Phys. Rev. A 84, 014701 (2011).
    """

    rep = int(rep)
    if rep not in (0, 1, 2):
        raise ValueError("rep must be 0 (diabatic), 1 (adiabatic/gap), or 2 (adiabatic/NAC)")
    ham_dia = _square_matrix(_ham_value(ham, "ham_dia"), "ham_dia").real
    ham_dia_prev = _square_matrix(
        _ham_value(ham_prev, "ham_dia"), "ham_prev.ham_dia"
    ).real
    if ham_dia_prev.shape != ham_dia.shape:
        raise ValueError("current and previous diabatic Hamiltonians must match")
    nstates = ham_dia.shape[0]
    i = _state_index(active_state, nstates)
    result = np.zeros(nstates, dtype=float)

    if rep == 0:
        derivatives = np.asarray(_ham_value(ham, "d1ham_dia"), dtype=complex).real
        if derivatives.ndim != 3 or derivatives.shape[1:] != ham_dia.shape:
            raise ValueError("d1ham_dia must have shape (ndof, nstates, nstates)")
        velocity = _nuclear_vector(momentum, "momentum") * _nuclear_vector(
            inverse_mass, "inverse_mass"
        )
        if len(velocity) != derivatives.shape[0]:
            raise ValueError("nuclear arrays must match the d1ham_dia DOF dimension")
        for j in range(nstates):
            if j == i or not _gap_crossed(ham_dia, ham_dia_prev, i, j):
                continue
            denominator = np.sum(
                np.abs(velocity * (derivatives[:, i, i] - derivatives[:, j, j]))
            )
            coupling = ham_dia[i, j]
            if denominator > 0.0:
                result[j] = np.exp(-2.0 * np.pi * coupling * coupling / denominator)
    else:
        ham_adi = _square_matrix(_ham_value(ham, "ham_adi"), "ham_adi").real
        nac_adi = _square_matrix(_ham_value(ham, "nac_adi"), "nac_adi").real
        if ham_adi.shape != (nstates, nstates) or nac_adi.shape != ham_adi.shape:
            raise ValueError("adiabatic and diabatic state dimensions must match")
        nac_adi_prev = None
        if rep == 2:
            nac_adi_prev = _square_matrix(
                _ham_value(ham_prev, "nac_adi"), "ham_prev.nac_adi"
            ).real
        for j in range(nstates):
            if j == i:
                continue
            crossing = _gap_crossed(ham_dia, ham_dia_prev, i, j)
            if rep == 2:
                crossing = nac_adi[i, j] * nac_adi_prev[i, j] < 0.0
            nac = abs(nac_adi[i, j])
            if crossing and nac > 0.0:
                gap = abs(ham_adi[i, i] - ham_adi[j, j])
                result[j] = np.exp(-0.25 * np.pi * gap / nac)
    result[i] = 1.0 - result.sum()
    return result


def hopping_probabilities_zn(
    ham, ham_prev, active_state, rep, momentum, inverse_mass
):
    """Compute multidimensional Zhu–Nakamura (ZN) probabilities.

    This follows the Yu et al. multidimensional formula used in the C++ code.
    ``ham`` and ``ham_prev`` are mappings or objects containing ``ham_dia``.
    The current record must additionally provide either ``forces_adi`` with
    shape ``(ndof, nstates)`` or ``d1ham_adi`` with shape
    ``(ndof, nstates, nstates)``; in the latter case state forces are obtained
    as ``F_ki = -Re[dH_ii/dR_k]``.

    A proposal is evaluated only at a diabatic-gap sign change.  For active
    state ``i`` and target ``j``, define

        x = sqrt(|Σ_k F_ki M_k⁻¹ F_kj|),
        y = sqrt(|Σ_k (F_ki-F_kj)² M_k⁻¹|),
        h = |H_ij|,
        a² = x y / (16 h³),
        b² = y / (2 x h),
        s = sign(F_i · F_j).

    Then

        g_ij = exp{-[π/(4 sqrt(a²))]
                    sqrt[2/(b² + sqrt(b⁴ + s))]}.

    The C++ signature includes ``rep`` and ``momentum``, although its current
    ZN implementation does not use them; they remain in this API for parity.
    ``inverse_mass`` is required for the mass-weighted force contractions.

    The formulation corresponds to the multidimensional Zhu–Nakamura
    equations cited in the C++ source (Yu et al.; Eqs. 3.89–3.91 in the Libra
    dynamics notes).
    """

    del rep, momentum
    ham_dia = _square_matrix(_ham_value(ham, "ham_dia"), "ham_dia").real
    ham_dia_prev = _square_matrix(
        _ham_value(ham_prev, "ham_dia"), "ham_prev.ham_dia"
    ).real
    if ham_dia_prev.shape != ham_dia.shape:
        raise ValueError("current and previous diabatic Hamiltonians must match")
    nstates = ham_dia.shape[0]
    i = _state_index(active_state, nstates)
    inverse_mass = _nuclear_vector(inverse_mass, "inverse_mass")
    forces = _adiabatic_forces(ham, len(inverse_mass), nstates)
    result = np.zeros(nstates, dtype=float)

    force_i = forces[:, i]
    for j in range(nstates):
        if j == i or not _gap_crossed(ham_dia, ham_dia_prev, i, j):
            continue
        force_j = forces[:, j]
        sign = np.sign(np.dot(force_i, force_j))
        x = np.sqrt(abs(np.sum(force_i * inverse_mass * force_j)))
        delta_force = force_i - force_j
        y = np.sqrt(abs(np.sum(delta_force * inverse_mass * delta_force)))
        coupling = abs(ham_dia[i, j])
        if coupling <= 0.0 or x <= 0.0:
            continue
        a2 = 0.0625 * x * y / coupling**3
        b2 = 0.5 * y / (x * coupling)
        inner = b2 * b2 + sign
        denominator = b2 + np.sqrt(max(0.0, inner))
        if a2 > 0.0 and denominator > 0.0:
            result[j] = np.exp(
                -(0.25 * np.pi / np.sqrt(a2)) * np.sqrt(2.0 / denominator)
            )
    result[i] = 1.0 - result.sum()
    return result


def _mssh_vector(params, populations, hvib, active_state):
    i = _state_index(active_state, len(populations))
    result = np.zeros(len(populations), dtype=float)
    for j in range(len(populations)):
        if j != i:
            result[j] = populations[j] * _boltzmann_factor(
                params, hvib[j, j].real - hvib[i, i].real
            )
    result[i] = 1.0 - result.sum()
    return result


def hopping_probabilities_mash(params, density):
    """Compute mapped surface-hopping (MASH) proposal probabilities.

    The proposal is deterministic: the state with the largest diagonal density
    receives probability one and all other states receive zero.  NumPy's first
    maximum is used to break ties, matching the forward scan in C++.

    References
    ----------
    J. R. Mannouch and J. O. Richardson, J. Chem. Phys. 158, 104111 (2023);
    J. E. Runeson and D. E. Manolopoulos, J. Chem. Phys. 159, 094115 (2023);
    J. E. Runeson, T. P. Fay, and D. E. Manolopoulos, Phys. Chem. Chem. Phys.
    26, 4929–4938 (2024).
    """

    del params
    density = _square_matrix(density, "density")
    result = np.zeros(density.shape[0], dtype=float)
    result[int(np.argmax(np.diag(density).real))] = 1.0
    return result


def hop(*args):
    """Select a final state from a probability vector or transition matrix.

    The probability vector is renormalized before sampling.  The selected
    state is the interval containing ``ksi`` in the cumulative distribution.
    A zero-norm vector leaves the supplied initial state unchanged; in the
    two-argument form that default state is zero.  ``ksi >= 1`` selects the
    highest state, preserving the C++ rounding-error safeguard.

    Supported forms mirror the C++ overloads: ``hop(probabilities, ksi)`` and
    ``hop(initial_state, probabilities, ksi)``.  For a probability matrix, the
    row belonging to ``initial_state`` is sampled.
    """

    if len(args) == 2:
        initial_state, probabilities, ksi = 0, args[0], args[1]
    elif len(args) == 3:
        initial_state, probabilities, ksi = args
    else:
        raise TypeError("hop expects (probabilities, ksi) or (initial_state, probabilities, ksi)")

    probabilities = np.asarray(probabilities, dtype=float)
    if probabilities.ndim == 2:
        initial_state = _state_index(initial_state, probabilities.shape[0])
        probabilities = probabilities[initial_state]
    if probabilities.ndim != 1 or probabilities.size == 0:
        raise ValueError("probabilities must be a nonempty vector or square matrix")
    if np.any(probabilities < 0.0) or not np.all(np.isfinite(probabilities)):
        raise ValueError("probabilities must be finite and nonnegative")
    if ksi >= 1.0:
        return probabilities.size - 1
    if ksi < 0.0:
        raise ValueError("ksi must be in [0, 1]")
    norm = probabilities.sum()
    if norm <= 0.0:
        return int(initial_state)
    return min(int(np.searchsorted(np.cumsum(probabilities / norm), ksi, side="right")), probabilities.size - 1)


def propose_hops(probabilities, active_states, rng=None):
    """Draw proposed final states for a batch of trajectories.

    One uniform random number in ``[0, 1)`` is drawn per trajectory and passed
    to :func:`hop`.  ``rng`` may be any object implementing NumPy's
    ``uniform(low, high)`` interface; a fresh default generator is used when it
    is omitted.
    """

    active_states = np.asarray(active_states, dtype=int)
    if active_states.ndim != 1:
        raise ValueError("active_states must be one-dimensional")
    if len(probabilities) != len(active_states):
        raise ValueError("one probability array is required per active state")
    rng = rng or np.random.default_rng()
    return np.asarray(
        [hop(state, prob, rng.uniform(0.0, 1.0)) for prob, state in zip(probabilities, active_states)],
        dtype=int,
    )


def _boltzmann_factor(params, energy_gap):
    if not int(_param(params, "use_boltz_factor", 0)) or energy_gap <= 0.0:
        return 1.0
    temperature = float(_param(params, "Temperature", 300.0))
    if temperature <= 0.0:
        return 0.0
    exponent = -energy_gap / (KB_HARTREE_PER_K * temperature)
    return float(np.exp(exponent)) if exponent > -500.0 else 0.0


def _param(params: Any, name: str, default):
    return params.get(name, default) if isinstance(params, dict) else getattr(params, name, default)


def _square_matrix(value, name):
    array = np.asarray(value, dtype=complex)
    if array.ndim != 2 or array.shape[0] != array.shape[1]:
        raise ValueError(f"{name} must be a square matrix")
    return array


def _density_matrix(value, nstates):
    array = _square_matrix(value, "density")
    if array.shape != (nstates, nstates):
        raise ValueError("density and Hamiltonian dimensions must match")
    return array


def _amplitudes(value, nstates):
    array = np.asarray(value, dtype=complex)
    if array.ndim == 2 and 1 in array.shape:
        array = array.reshape(-1)
    if array.ndim != 1 or array.size != nstates:
        raise ValueError("amplitudes must be a vector matching the Hamiltonian dimension")
    return array


def _state_index(value, nstates):
    index = int(value)
    if index != value or not 0 <= index < nstates:
        raise ValueError(f"state index must be in [0, {nstates})")
    return index


def _ham_value(ham, name):
    if isinstance(ham, dict):
        if name not in ham:
            raise ValueError(f"Hamiltonian record is missing {name!r}")
        return ham[name]
    if not hasattr(ham, name):
        raise ValueError(f"Hamiltonian object is missing {name!r}")
    return getattr(ham, name)


def _nuclear_vector(value, name):
    array = np.asarray(value, dtype=float)
    if array.ndim == 2 and 1 in array.shape:
        array = array.reshape(-1)
    if array.ndim != 1:
        raise ValueError(f"{name} must be a one-dimensional nuclear vector")
    return array


def _gap_crossed(ham_dia, ham_dia_prev, state_i, state_j):
    gap = ham_dia[state_i, state_i] - ham_dia[state_j, state_j]
    gap_prev = ham_dia_prev[state_i, state_i] - ham_dia_prev[state_j, state_j]
    return gap * gap_prev < 0.0


def _adiabatic_forces(ham, ndof, nstates):
    try:
        forces = np.asarray(_ham_value(ham, "forces_adi"), dtype=float)
    except ValueError:
        derivatives = np.asarray(_ham_value(ham, "d1ham_adi"), dtype=complex)
        if derivatives.shape != (ndof, nstates, nstates):
            raise ValueError(
                "d1ham_adi must have shape (ndof, nstates, nstates)"
            ) from None
        forces = -np.diagonal(derivatives.real, axis1=1, axis2=2)
    if forces.shape == (nstates, ndof):
        forces = forces.T
    if forces.shape != (ndof, nstates):
        raise ValueError("forces_adi must have shape (ndof, nstates)")
    return forces


def _history_probabilities(
    function, ham, ham_prev, active_states, rep, momentum, inverse_mass
):
    states = np.asarray(active_states, dtype=int)
    if states.ndim == 0:
        return function(
            ham, ham_prev, int(states), rep, momentum, inverse_mass
        )
    if states.ndim != 1:
        raise ValueError("active_states must be a scalar or one-dimensional")
    if not isinstance(ham, (list, tuple)) or not isinstance(ham_prev, (list, tuple)):
        raise ValueError("batched LZ/ZN requires ham and ham_prev sequences")
    if len(ham) != len(states) or len(ham_prev) != len(states):
        raise ValueError("Hamiltonian-history batch sizes must match active_states")
    momenta = np.asarray(momentum, dtype=float)
    if momenta.ndim != 2:
        raise ValueError("batched momentum must have shape (ntraj, ndof) or (ndof, ntraj)")
    if momenta.shape[0] != len(states) and momenta.shape[1] == len(states):
        momenta = momenta.T
    if momenta.shape[0] != len(states):
        raise ValueError("momentum trajectory dimension must match active_states")
    return np.asarray(
        [
            function(current, previous, state, rep, p, inverse_mass)
            for current, previous, state, p in zip(ham, ham_prev, states, momenta)
        ]
    )
