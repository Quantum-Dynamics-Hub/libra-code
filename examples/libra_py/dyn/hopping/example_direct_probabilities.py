"""Concrete FSSH, GFSH, FSSH2/3, MSSH, and MASH probabilities."""

from _example_setup import use_repo_sources

use_repo_sources()

import numpy as np

from libra_py.dyn.control_params import DynControlParams
from libra_py.dyn.hopping import (
    hopping_probabilities_fssh,
    hopping_probabilities_fssh2,
    hopping_probabilities_fssh3,
    hopping_probabilities_gfsh,
    hopping_probabilities_gfsh_orig,
    hopping_probabilities_mash,
    hopping_probabilities_mssh,
)


def show(label, values):
    print(f"{label:34s} {np.asarray(values)}  sum={np.sum(values):.12f}")


def main():
    hvib = np.array([[0.0, -0.1j], [0.1j, 0.01]], dtype=complex)
    density = np.array([[0.5, 1.0], [0.0, 0.5]], dtype=complex)

    plain = DynControlParams(dt=0.25)
    boltzmann = DynControlParams(dt=41.0, Temperature=300.0, use_boltz_factor=1)
    show("FSSH, no Boltzmann scaling", hopping_probabilities_fssh(plain, density, hvib, 0))
    show(
        "FSSH, 300 K uphill scaling",
        hopping_probabilities_fssh(boltzmann, density, hvib, 0),
    )

    amplitudes = np.array([np.sqrt(0.8), np.sqrt(0.2)], dtype=complex)
    print("\nFSSH all-pairs matrix from amplitudes:")
    print(hopping_probabilities_fssh(DynControlParams(dt=0.1), amplitudes, hvib))

    gfsh_density = np.array([[0.5, 0.5], [0.5, 0.5]], dtype=complex)
    show("GFSH instantaneous flux", hopping_probabilities_gfsh(plain, gfsh_density, hvib, 0))

    old = np.diag([0.8, 0.2]).astype(complex)
    new = np.diag([0.6, 0.4]).astype(complex)
    show("Original GFSH, population step", hopping_probabilities_gfsh_orig({}, new, old, 0))
    show("FSSH2 revision 0", hopping_probabilities_fssh2({"fssh2_revision": 0}, new, old, 0))
    show("FSSH2 revision 1", hopping_probabilities_fssh2({"fssh2_revision": 1}, new, old, 0))

    fssh3_params = DynControlParams(
        dt=1.0,
        fssh3_dt=0.01,
        fssh3_max_steps=2000,
        fssh3_err_tol=1.0e-12,
    )
    errors = [0.0] * 5
    show("FSSH3", hopping_probabilities_fssh3(fssh3_params, new, old, 0, errors))
    print("FSSH3 errors [L, L0, L1, L2, L3]:", errors)

    population_density = np.diag([0.5, 0.3, 0.2]).astype(complex)
    energies = np.diag([0.0, 0.01, -0.005]).astype(complex)
    show("MSSH, three states", hopping_probabilities_mssh({}, population_density, energies, 0))
    show(
        "MSSH, Boltzmann scaled",
        hopping_probabilities_mssh(
            DynControlParams(use_boltz_factor=1, Temperature=300.0),
            population_density,
            energies,
            0,
        ),
    )
    show("MASH", hopping_probabilities_mash({}, np.diag([0.2, 0.7, 0.1])))


if __name__ == "__main__":
    main()
