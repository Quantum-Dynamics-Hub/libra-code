"""Hop acceptance decisions and post-hop momentum rescaling."""

import numpy as np

from libra_py.dyn.control_params import DynControlParams
from libra_py.dyn.hopping import (
    accept_hops,
    can_rescale_along_vector,
    handle_hops_nuclear,
    rescale_along_vector,
    where_can_we_hop,
)


def main():
    print("Energy-conserving rescaling along a coupling direction")
    direction = np.array([1.0])
    inverse_mass = np.array([1.0])
    for momentum in (0.5, 1.0, 2.0):
        possible = can_rescale_along_vector(
            0.0, 1.0, [momentum], inverse_mass, direction
        )
        print(f"p={momentum:3.1f}, 0 -> 1 Ha uphill: possible={possible}")

    momentum = np.array([0.0])
    rescale_along_vector(1.0, 0.0, momentum, inverse_mass, direction)
    print("downhill 1 Ha produces p=", momentum, "and KE=", 0.5 * momentum[0] ** 2)

    frustrated = np.array([1.0])
    rescale_along_vector(
        0.0, 1.0, frustrated, inverse_mass, direction, do_reverse=True
    )
    print("frustrated uphill hop with reversal gives p=", frustrated)

    energies = np.array([0.0, 0.2, 1.0])
    initial = np.array([0])
    proposed = np.array([2])
    energy_params = DynControlParams(hop_acceptance_algo=10)
    rejected = accept_hops(
        energy_params,
        proposed,
        initial,
        energies,
        momenta=[[1.0]],
        inverse_mass=inverse_mass,
    )
    print("\n0 -> 2 proposal with KE=0.5 Ha:", rejected, "(rejected)")
    possible = where_can_we_hop(
        0,
        energy_params,
        initial,
        energies,
        momenta=[[1.0]],
        inverse_mass=inverse_mass,
    )
    print("energy-accessible target states:", possible)

    dc1 = np.array([[[[0.0, 1.0], [-1.0, 0.0]]]])
    dc_params = DynControlParams(hop_acceptance_algo=20)
    for momentum_value in (1.0, 2.0):
        final = accept_hops(
            dc_params,
            [1],
            [0],
            [0.0, 1.0],
            momenta=[[momentum_value]],
            inverse_mass=inverse_mass,
            dc1_adi=dc1,
        )
        print(
            f"derivative-coupling acceptance with p={momentum_value}: final state {final[0]}"
        )

    momenta = np.array([[2.0]])
    handle_hops_nuclear(
        DynControlParams(momenta_rescaling_algo=100),
        momenta,
        inverse_mass,
        new_states=[1],
        old_states=[0],
        energies_adi=[0.0, 1.0],
    )
    print("\nuniform rescaling after accepted 0 -> 1 hop:", momenta)

    momenta = np.array([[2.0]])
    handle_hops_nuclear(
        DynControlParams(
            momenta_rescaling_algo=200, quantum_dofs=[0]
        ),
        momenta,
        inverse_mass,
        new_states=[1],
        old_states=[0],
        energies_adi=[0.0, 1.0],
        dc1_adi=dc1,
    )
    print("coupling-direction rescaling after the same hop:", momenta)

    tcnbra_kinetic = np.array([2.0])
    handle_hops_nuclear(
        DynControlParams(momenta_rescaling_algo=40),
        momenta,
        inverse_mass,
        new_states=[1],
        old_states=[0],
        energies_adi=[0.0, 0.5],
        tcnbra_ekin=tcnbra_kinetic,
    )
    print("TC-NBRA kinetic energy after 0.5 Ha uphill hop:", tcnbra_kinetic)


if __name__ == "__main__":
    main()
