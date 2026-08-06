"""Hop acceptance and nuclear momentum handling.

NumPy translation of ``src/dyn/dyn_hop_acceptance.cpp``. Energies and times
use atomic units, temperatures use kelvin, momenta are arranged as
``(ntraj, ndof)``, and Hamiltonian data as ``(ntraj, nstates, nstates)`` (a
single shared matrix is also accepted).
"""

from __future__ import annotations

from math import erf, exp, pi, sqrt

import numpy as np


KB_HARTREE_PER_K = 3.166811563e-6
HOP_ACCEPTANCE_ALGOS = {
    0: "accept all",
    10: "adiabatic energy",
    11: "diabatic energy",
    20: "derivative-coupling rescaling feasibility",
    21: "force-difference rescaling feasibility",
    31: "quantum Boltzmann ratio",
    32: "classical Maxwell-Boltzmann probability",
    33: "normalized quantum final-state probability",
    40: "TC-NBRA kinetic energy",
}
MOMENTA_RESCALING_ALGOS = {
    0: "none",
    40: "TC-NBRA kinetic energy",
    100: "adiabatic energy, no reversal",
    101: "adiabatic energy, reverse frustrated hop",
    110: "diabatic energy, no reversal",
    111: "diabatic energy, reverse frustrated hop",
    200: "derivative coupling, no reversal",
    201: "derivative coupling, reverse frustrated hop",
    210: "force difference, no reversal",
    211: "force difference, reverse frustrated hop",
}


def can_rescale_along_vector(
    energy_old, energy_new, momentum, inverse_mass, direction, which_dofs=None
):
    """Return whether momentum can be rescaled along ``direction``.

    Following Fabiano et al., solve ``a γ² - b γ - c = 0`` with

        a = 1/2 Σ_k t_k² M_k⁻¹,
        b = Σ_k p_k t_k M_k⁻¹,
        c = E_old - E_new.

    A real energy-conserving rescaling exists when
    ``D = b² + 4ac >= 0``. Only ``which_dofs`` enter these contractions.
    """

    _, _, _, discriminant = _rescaling_terms(
        energy_old, energy_new, momentum, inverse_mass, direction, which_dofs
    )
    return discriminant >= 0.0


def rescale_along_vector(
    energy_old,
    energy_new,
    momentum,
    inverse_mass,
    direction,
    do_reverse=False,
    which_dofs=None,
):
    """Rescale momentum along a vector while conserving total energy.

    For an allowed hop the quadratic root causing the smallest momentum change
    is selected. For a frustrated hop, ``do_reverse=True`` reverses the
    component along ``direction``. The supplied momentum array is updated in
    place and returned. As in C++, a direction with ``|a| <= 1e-12`` causes no
    change.
    """

    p, t, a, discriminant = _rescaling_terms(
        energy_old, energy_new, momentum, inverse_mass, direction, which_dofs
    )
    indices = _dof_indices(which_dofs, len(p))
    inv_mass = _vector(inverse_mass, "inverse_mass")
    b = float(np.sum(p[indices] * t[indices] * inv_mass[indices]))
    gamma = 0.0
    if discriminant < 0.0:
        if abs(a) > 1.0e-12 and do_reverse:
            gamma = b / a
    elif abs(a) > 1.0e-12:
        gamma = 0.5 * (b + sqrt(discriminant)) / a if b < 0.0 else 0.5 * (b - sqrt(discriminant)) / a
    p[indices] -= gamma * t[indices]
    return p


def Boltz_quant_prob(energies, temperature):
    """Normalized quantum Boltzmann probabilities ``exp(-E_i/kT)/Z``."""

    energies = _vector(energies, "energies")
    beta = _beta(temperature)
    weights = np.exp(-(energies - energies.min()) * beta)
    return weights / weights.sum()


def Boltz_cl_prob(energy, temperature):
    """Classical Maxwell distribution density used by the C++ routine.

    This retains the original expression
    ``(2β/sqrt(π)) sqrt(Eβ) exp[-sqrt(Eβ)]``. For an integrated probability of
    exceeding a threshold use :func:`Boltz_cl_prob_up`.
    """

    if energy < 0.0:
        raise ValueError("energy must be nonnegative")
    beta = _beta(temperature)
    x = sqrt(energy * beta)
    return (2.0 * beta / sqrt(pi)) * x * exp(-x)


def Boltz_cl_prob_up(energy, temperature):
    """Probability that classical kinetic energy exceeds ``energy``.

    ``P(E_k > E) = 1 - erf(x) + 2x exp(-x²)/sqrt(π)``, where
    ``x = sqrt(E/kT)``.
    """

    if energy < 0.0:
        raise ValueError("energy must be nonnegative")
    x = sqrt(energy * _beta(temperature))
    return 1.0 - erf(x) + sqrt(4.0 / pi) * x * exp(-x * x)


def HO_prob(energies, quantum_numbers, temperature, prob=None):
    """Probability of specified quantum states of independent oscillators.

    For oscillator ``i``, ``p_i(n_i) = ξ_i**n_i (1-ξ_i)`` with
    ``ξ_i = exp(-E_i/kT)``. Returns ``(product(p_i), p_i)``; if the optional
    mutable ``prob`` is supplied it is updated in place for C++ API parity.
    """

    energies, quantum_numbers = _oscillator_inputs(energies, quantum_numbers)
    xi = np.exp(-energies * _beta(temperature))
    probabilities = xi**quantum_numbers * (1.0 - xi)
    _store_vector(prob, probabilities)
    return float(np.prod(probabilities)), probabilities


def HO_prob_up(energies, quantum_numbers, temperature, prob=None):
    """Probability that each oscillator lies above its specified level.

    The public C++ implementation uses ``p_i = exp(-E_i/kT)`` independently of
    ``quantum_numbers``; this behavior is preserved. Returns
    ``(product(p_i), p_i)``.
    """

    energies, _ = _oscillator_inputs(energies, quantum_numbers)
    probabilities = np.exp(-energies * _beta(temperature))
    _store_vector(prob, probabilities)
    return float(np.prod(probabilities)), probabilities


def boltz_factor(energy_new, energy_old, temperature, boltz_opt):
    """Return an uphill-hop acceptance factor.

    Options match C++: 0 accepts all; 1 uses ``exp(-ΔE/kT)``; 2 uses the
    classical energy-tail probability; and 3 uses the normalized population of
    ``ΔE`` in a two-level quantum system. Downhill hops always return one.
    """

    option = int(boltz_opt)
    if option not in (0, 1, 2, 3):
        raise ValueError("boltz_opt must be 0, 1, 2, or 3")
    gap = float(energy_new) - float(energy_old)
    if option == 0 or gap <= 0.0:
        return 1.0
    if option == 1:
        exponent = gap * _beta(temperature)
        return 0.0 if exponent > 50.0 else exp(-exponent)
    if option == 2:
        return Boltz_cl_prob_up(gap, temperature)
    return float(Boltz_quant_prob([0.0, gap], temperature)[1])


def accept_hops(
    params,
    proposed_states,
    initial_states,
    energies_adi,
    *,
    momenta=None,
    inverse_mass=None,
    energies_dia=None,
    dc1_adi=None,
    d1ham_adi=None,
    tcnbra_ekin=None,
    rng=None,
    which_trajectories=None,
):
    """Accept or reject proposed hops using the C++ option numbering.

    Options 0, 10, 11, 20, 21, 31, 32, 33, and 40 are implemented. Matrix
    inputs may be shared across trajectories or batched. Entries not listed in
    ``which_trajectories`` are returned as ``-1``, matching the C++ overload.
    """

    proposed = np.asarray(proposed_states, dtype=int)
    initial = np.asarray(initial_states, dtype=int)
    if proposed.ndim != 1 or proposed.shape != initial.shape:
        raise ValueError("proposed_states and initial_states must be matching vectors")
    ntraj = len(initial)
    selected = np.arange(ntraj) if which_trajectories is None else np.asarray(which_trajectories, dtype=int)
    final = np.full(ntraj, -1, dtype=int)
    algorithm = int(_param(params, "hop_acceptance_algo", 0))
    if algorithm not in HOP_ACCEPTANCE_ALGOS:
        raise ValueError(f"unsupported hop_acceptance_algo {algorithm}")
    eadi = _energy_batch(energies_adi, ntraj, "energies_adi")
    edia = _energy_batch(energies_dia, ntraj, "energies_dia") if energies_dia is not None else None
    p = _momentum_batch(momenta, ntraj) if momenta is not None else None
    inv_mass = _vector(inverse_mass, "inverse_mass") if inverse_mass is not None else None
    dofs = _param(params, "quantum_dofs", None)
    rng = rng or np.random.default_rng()

    for traj in selected:
        old, new = int(initial[traj]), int(proposed[traj])
        accepted = True
        if old != new and algorithm in (10, 11):
            _require(p=p, inverse_mass=inv_mass)
            energies = eadi if algorithm == 10 else _require(energies_dia=edia)["energies_dia"]
            accepted = _kinetic_energy(p[traj], inv_mass, dofs) + energies[traj, old] - energies[traj, new] >= 0.0
        elif old != new and algorithm in (20, 21):
            _require(p=p, inverse_mass=inv_mass)
            if algorithm == 20:
                direction = _pair_vector(dc1_adi, traj, old, new, ntraj, "dc1_adi")
            else:
                deriv = _tensor_batch(d1ham_adi, ntraj, "d1ham_adi")[traj]
                direction = np.diagonal(deriv.real, axis1=1, axis2=2)[:, old] - np.diagonal(deriv.real, axis1=1, axis2=2)[:, new]
            accepted = can_rescale_along_vector(eadi[traj, old], eadi[traj, new], p[traj], inv_mass, direction, dofs)
        elif algorithm in (31, 32, 33):
            probability = boltz_factor(eadi[traj, new], eadi[traj, old], _param(params, "Temperature", 300.0), algorithm - 30)
            accepted = rng.uniform(0.0, 1.0) < probability
        elif algorithm == 40:
            kinetic = _require(tcnbra_ekin=tcnbra_ekin)["tcnbra_ekin"]
            accepted = eadi[traj, old] + kinetic[traj] - eadi[traj, new] >= 0.0
        final[traj] = new if accepted else old
    return final


def where_can_we_hop(traj, params, initial_states, energies_adi, **kwargs):
    """Return nontrivial target states accepted for one trajectory."""

    initial = np.asarray(initial_states, dtype=int)
    nstates = _energy_batch(energies_adi, len(initial), "energies_adi").shape[1]
    possible = []
    for target in range(nstates):
        if target == initial[traj]:
            continue
        proposed = initial.copy()
        proposed[traj] = target
        result = accept_hops(
            params, proposed, initial, energies_adi, which_trajectories=[traj], **kwargs
        )
        if result[traj] != initial[traj]:
            possible.append(int(result[traj]))
    return possible


def handle_hops_nuclear(
    params,
    momenta,
    inverse_mass,
    new_states,
    old_states,
    energies_adi,
    *,
    energies_dia=None,
    dc1_adi=None,
    d1ham_adi=None,
    tcnbra_ekin=None,
):
    """Apply C++ momentum-rescaling options after accepted/frustrated hops.

    Momentum is updated in place and returned. If option 40 is selected,
    ``tcnbra_ekin`` is updated in place instead.
    """

    p = np.asarray(momenta, dtype=float)
    if p.ndim == 1:
        p = p.reshape(1, -1)
    inv_mass = _vector(inverse_mass, "inverse_mass")
    old_states = np.asarray(old_states, dtype=int)
    new_states = np.asarray(new_states, dtype=int)
    ntraj = len(old_states)
    if p.shape[0] != ntraj or new_states.shape != old_states.shape:
        raise ValueError("trajectory dimensions must match")
    eadi = _energy_batch(energies_adi, ntraj, "energies_adi")
    edia = _energy_batch(energies_dia, ntraj, "energies_dia") if energies_dia is not None else None
    algorithm = int(_param(params, "momenta_rescaling_algo", 0))
    if algorithm not in MOMENTA_RESCALING_ALGOS:
        raise ValueError(f"unsupported momenta_rescaling_algo {algorithm}")
    dofs = _param(params, "quantum_dofs", None)

    for traj, (old, new) in enumerate(zip(old_states, new_states)):
        if old == new or algorithm == 0:
            continue
        if algorithm == 40:
            kinetic = _require(tcnbra_ekin=tcnbra_ekin)["tcnbra_ekin"]
            kinetic[traj] += eadi[traj, old] - eadi[traj, new]
        elif algorithm in (100, 101, 110, 111):
            energies = eadi if algorithm < 110 else _require(energies_dia=edia)["energies_dia"]
            kinetic_old = _kinetic_energy(p[traj], inv_mass, dofs)
            kinetic_new = kinetic_old + energies[traj, old] - energies[traj, new]
            scale = 1.0
            if kinetic_old > 0.0:
                scale = sqrt(kinetic_new / kinetic_old) if kinetic_new >= 0.0 else (-1.0 if algorithm in (101, 111) else 1.0)
            p[traj, _dof_indices(dofs, p.shape[1])] *= scale
        elif algorithm in (200, 201, 210, 211):
            if algorithm < 210:
                direction = _pair_vector(dc1_adi, traj, old, new, ntraj, "dc1_adi")
                reverse = algorithm == 201
                if reverse and int(_param(params, "use_Jasper_Truhlar_criterion", 1)):
                    deriv = _tensor_batch(d1ham_adi, ntraj, "d1ham_adi")[traj]
                    diag = np.diagonal(deriv.real, axis1=1, axis2=2)
                    force_old, force_new = -diag[:, old], -diag[:, new]
                    velocity = p[traj] * inv_mass
                    reverse = np.dot(force_old, direction) * np.dot(force_new, direction) < 0.0 and np.dot(velocity, direction) * np.dot(force_new, direction) < 0.0
            else:
                deriv = _tensor_batch(d1ham_adi, ntraj, "d1ham_adi")[traj]
                diag = np.diagonal(deriv.real, axis1=1, axis2=2)
                direction = diag[:, old] - diag[:, new]
                reverse = algorithm == 211
            rescale_along_vector(eadi[traj, old], eadi[traj, new], p[traj], inv_mass, direction, reverse, dofs)
    return momenta


def _rescaling_terms(e_old, e_new, momentum, inverse_mass, direction, which_dofs):
    p = _vector(momentum, "momentum")
    inv_mass = _vector(inverse_mass, "inverse_mass")
    t = _vector(direction, "direction")
    if not (len(p) == len(inv_mass) == len(t)):
        raise ValueError("momentum, inverse_mass, and direction dimensions must match")
    idx = _dof_indices(which_dofs, len(p))
    a = 0.5 * float(np.sum(t[idx] ** 2 * inv_mass[idx]))
    b = float(np.sum(p[idx] * t[idx] * inv_mass[idx]))
    return p, t, a, b * b + 4.0 * a * (float(e_old) - float(e_new))


def _beta(temperature):
    if temperature <= 0.0:
        raise ValueError("temperature must be positive")
    return 1.0 / (KB_HARTREE_PER_K * temperature)


def _vector(value, name):
    array = np.asarray(value, dtype=float).reshape(-1)
    if array.size == 0:
        raise ValueError(f"{name} must not be empty")
    return array


def _dof_indices(which_dofs, ndof):
    indices = np.arange(ndof) if which_dofs is None or len(which_dofs) == 0 else np.asarray(which_dofs, dtype=int)
    if np.any(indices < 0) or np.any(indices >= ndof):
        raise ValueError("which_dofs contains an invalid index")
    return indices


def _oscillator_inputs(energies, quantum_numbers):
    energies = _vector(energies, "energies")
    qn = np.asarray(quantum_numbers, dtype=int).reshape(-1)
    if qn.shape != energies.shape or np.any(qn < 0):
        raise ValueError("quantum_numbers must be nonnegative and match energies")
    return energies, qn


def _store_vector(target, values):
    if target is None:
        return
    if hasattr(target, "clear") and hasattr(target, "extend"):
        target.clear(); target.extend(values.tolist())
    else:
        target[...] = values


def _param(params, name, default):
    return params.get(name, default) if isinstance(params, dict) else getattr(params, name, default)


def _energy_batch(value, ntraj, name):
    array = np.asarray(value, dtype=complex).real
    if array.ndim == 2 and array.shape[0] == array.shape[1]:
        array = np.broadcast_to(np.diag(array), (ntraj, array.shape[0]))
    elif array.ndim == 3:
        array = np.diagonal(array, axis1=1, axis2=2)
    elif array.ndim == 1:
        array = np.broadcast_to(array, (ntraj, len(array)))
    if array.ndim != 2 or array.shape[0] != ntraj:
        raise ValueError(f"{name} must contain one state-energy vector per trajectory")
    return array


def _momentum_batch(value, ntraj):
    array = np.asarray(value, dtype=float)
    if array.ndim == 1:
        array = array.reshape(1, -1)
    if array.ndim != 2:
        raise ValueError("momenta must have shape (ntraj, ndof)")
    if array.shape[0] != ntraj and array.shape[1] == ntraj:
        array = array.T
    if array.shape[0] != ntraj:
        raise ValueError("momenta trajectory dimension must match states")
    return array


def _tensor_batch(value, ntraj, name):
    array = np.asarray(value, dtype=complex)
    if array.ndim == 3:
        array = np.broadcast_to(array, (ntraj,) + array.shape)
    if array.ndim != 4 or array.shape[0] != ntraj:
        raise ValueError(f"{name} must have shape (ntraj, ndof, nstates, nstates)")
    return array


def _pair_vector(value, traj, old, new, ntraj, name):
    tensor = _tensor_batch(value, ntraj, name)[traj]
    return tensor[:, old, new].real


def _kinetic_energy(momentum, inverse_mass, which_dofs):
    idx = _dof_indices(which_dofs, len(momentum))
    return 0.5 * float(np.sum(momentum[idx] ** 2 * inverse_mass[idx]))


def _require(**values):
    missing = [name for name, value in values.items() if value is None]
    if missing:
        raise ValueError(f"required input missing: {', '.join(missing)}")
    return values
