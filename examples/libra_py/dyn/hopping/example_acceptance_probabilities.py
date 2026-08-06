"""Thermal probability functions used by stochastic hop acceptance."""

import numpy as np

from libra_py.dyn.hopping import (
    Boltz_cl_prob,
    Boltz_cl_prob_up,
    Boltz_quant_prob,
    HO_prob,
    HO_prob_up,
    boltz_factor,
)


def main():
    temperature = 300.0
    energy_gap = 0.01

    print("Two-level quantum Boltzmann distribution")
    print("energies [Ha]:", [0.0, energy_gap])
    print("probabilities:", Boltz_quant_prob([0.0, energy_gap], temperature))

    print("\nClassical Maxwell-Boltzmann functions")
    print("density-like C++ expression:", Boltz_cl_prob(energy_gap, temperature))
    print(
        "probability of kinetic energy above the gap:",
        Boltz_cl_prob_up(energy_gap, temperature),
    )

    print("\nUphill acceptance factors for a 0.01 Ha gap")
    descriptions = {
        0: "accept every hop",
        1: "ordinary Boltzmann ratio",
        2: "classical kinetic-energy tail",
        3: "normalized quantum final-state population",
    }
    for option, description in descriptions.items():
        probability = boltz_factor(energy_gap, 0.0, temperature, option)
        print(f"option {option}: {probability:.12g}  ({description})")
    print("downhill factor:", boltz_factor(0.0, energy_gap, temperature, 1))

    frequencies = [0.005, 0.01, 0.02]
    quantum_numbers = [0, 1, 2]
    total, individual = HO_prob(
        frequencies, quantum_numbers, temperature
    )
    print("\nIndependent harmonic oscillators")
    print("level spacings [Ha]:", frequencies)
    print("quantum numbers:    ", quantum_numbers)
    print("individual probabilities:", individual)
    print("joint probability:       ", total)

    total_up, individual_up = HO_prob_up(
        frequencies, quantum_numbers, temperature
    )
    print("upper-tail probabilities:", individual_up)
    print("joint upper-tail value:  ", total_up)
    np.testing.assert_allclose(total, np.prod(individual))


if __name__ == "__main__":
    main()
