"""Compute decoherence rates/times and apply coefficient corrections."""

import numpy as np

from libra_py.dyn.decoherence import (
    coherence_intervals,
    collapse,
    decoherence_event,
    dephasing_informed_correction,
    edc_rates,
    gu_franco,
    project_out,
    schwartz_1,
    schwartz_2,
    sdm,
)


def main():
    """Print a two-trajectory, three-state decoherence demonstration."""

    # All batched Python functions are trajectory-first. Energies and kinetic
    # energies are in atomic units, so returned rates are inverse atomic time.
    energies = np.array([[0.00, 0.10, 0.25], [0.00, 0.12, 0.20]])
    hvib = np.array([np.diag(row) for row in energies], dtype=complex)
    kinetic_energy = np.array([0.50, 0.25])
    amplitudes = np.array(
        [
            [np.sqrt(0.70), 1j * np.sqrt(0.20), np.sqrt(0.10)],
            [np.sqrt(0.20), np.sqrt(0.50), 1j * np.sqrt(0.30)],
        ],
        dtype=complex,
    )

    edc = edc_rates(hvib, kinetic_energy, C_param=1.0, eps_param=0.1)
    intervals = coherence_intervals(amplitudes, edc)
    print("EDC rates [1/a.u.]:\n", edc)
    print("DISH coherence intervals [a.u.]:\n", intervals)

    average_gaps = np.array(
        [[0.05, 0.11, 0.24], [0.11, 0.05, 0.16], [0.24, 0.16, 0.05]]
    )
    corrected = dephasing_informed_correction(edc, hvib, average_gaps)
    print("Dephasing-informed EDC rates:\n", corrected)

    # State-resolved forces have shape (trajectory, state, nuclear DOF).
    state_forces = np.array(
        [
            [[0.10, -0.05], [-0.15, 0.02], [0.05, 0.12]],
            [[0.08, -0.04], [-0.12, 0.03], [0.03, 0.10]],
        ]
    )
    inverse_alpha = np.array([4.0, 2.0])
    print("Schwartz-1 diagonal rates:\n", schwartz_1(amplitudes, state_forces, inverse_alpha))
    print("Schwartz-2 pair rates:\n", schwartz_2(state_forces, inverse_alpha))
    print("Gu--Franco rates:\n", gu_franco(amplitudes, 0.10, 300.0))

    # SDM returns a corrected copy and conserves the norm trajectory by
    # trajectory. Trajectory 0 is active on state 0; trajectory 1 on state 1.
    decayed = sdm(amplitudes, dt=2.0, active_states=[0, 1], rates=edc)
    print("SDM amplitudes:\n", decayed)
    print("SDM norms:", np.sum(np.abs(decayed) ** 2, axis=1))

    projected = amplitudes.copy()
    project_out(projected, state=2, trajectory=0)
    collapsed = amplitudes.copy()
    collapse(collapsed, state=1, trajectory=1, collapse_option=0)
    print("Trajectory 0 after projecting out state 2:", projected[0])
    print("Trajectory 1 after phase-preserving collapse:", collapsed[1])

    # Direct comparison (option 0) makes this example deterministic: only
    # state 1 of trajectory 0 and state 2 of trajectory 1 have expired.
    clocks = np.array([[0.2, 3.0, 0.1], [0.1, 0.2, 4.0]])
    thresholds = np.ones_like(clocks)
    events = decoherence_event(
        clocks, thresholds, option=0, rng=np.random.default_rng(4)
    )
    print("Selected DISH event states:", events)


if __name__ == "__main__":
    main()
