"""SD and CSF overlaps for doublet, singlet, and triplet model spaces.

Positive orbital labels denote alpha spin and negative labels beta spin.
For every case, this example shows the generic ``configs_and_T_matrix``
construction used underneath ``sd_and_csf_overlaps`` and checks that both
routes produce the same SD-to-CSF result.
"""

import numpy as np

from libra_py.citools import interfaces
from libra_py.citools import slatdet as sd


def is_closed_shell(det):
    """Return True when every occupied spatial orbital is doubly occupied."""
    return len(det) % 2 == 0 and all(-orb in det for orb in det)


def unpaired_spins(det):
    """List singly occupied orbitals and their spin labels."""
    occupied = set(det)
    return [
        (abs(orb), "alpha" if orb > 0 else "beta")
        for orb in det
        if -orb not in occupied
    ]


def spin_projections(S):
    """Return all Ms components, from -S through +S."""
    return [(-2 * S + 2 * i) / 2 for i in range(int(2 * S) + 1)]


def run_component(label, nelec, norb, reference_det, excitations, S, Ms):
    orbital_space = list(range(1, norb + 1))
    raw_configs = [tuple(reference_det)]
    raw_configs.extend(
        tuple(sd.make_excitation(reference_det, occ, vir))
        for occ, vir in excitations
    )

    sd_basis, T = interfaces.configs_and_T_matrix(
        raw_configs,
        active_space=orbital_space,
        orbital_space=orbital_space,
        nelec=nelec,
        S=S,
        Ms=Ms,
    )
    T = T.toarray()

    # Identity one-electron overlaps make every normalized CSF overlap itself
    # by one. They also make the printed matrices especially transparent.
    spin_orbital_overlap = np.eye(2 * norb)
    csf_overlap, sd_overlap = interfaces.sd_and_csf_overlaps(
        spin_orbital_overlap,
        lowest_orbital=1,
        highest_orbital=norb,
        nelec=nelec,
        homo_indx=max(abs(orb) for orb in reference_det),
        common_sd_basis=excitations,
        _active_space=orbital_space,
        S=S,
        Ms=Ms,
        reference_det=reference_det,
    )

    reconstructed_csf_overlap = T.conj().T @ sd_overlap @ T
    assert np.allclose(csf_overlap, reconstructed_csf_overlap)
    assert np.allclose(csf_overlap, np.eye(csf_overlap.shape[0]))

    open_shell_sector = S > 0
    if open_shell_sector:
        assert not any(is_closed_shell(det) for det in sd_basis)

    print("=" * 72)
    print(f"{label}: N={nelec}, orbitals={norb}, S={S}, Ms={Ms}")
    print(f"Reference determinant: {tuple(reference_det)}")
    print("SD basis:")
    for i, det in enumerate(sd_basis):
        print(f"  SD {i}: {det}")
        if S > 0:
            spins = ", ".join(
                f"orbital {orbital}: {spin}" for orbital, spin in unpaired_spins(det)
            )
            print(f"        unpaired spins: {spins}")
    print("SD-to-CSF transformation T (columns are CSFs):")
    print(T.real)
    print("SD overlap:")
    print(sd_overlap)
    print("CSF overlap = T^dagger S_SD T:")
    print(csf_overlap.real)
    if open_shell_sector:
        print("Closed-shell determinant absent: True")


MULTIPLETS = [
    # Doublets use explicit restricted-open-shell references.
    ("Doublet: 1 electron in 2 orbitals", 1, 2, [1], [[1, 2]], 0.5),
    ("Doublet: 3 electrons in 3 orbitals", 3, 3, [1, -1, 2], [[2, 3]], 0.5),

    # Singlets start with the closed-shell determinant (1, -1).
    ("Singlet: 2 electrons in 2 orbitals", 2, 2, [1, -1], [[1, 2]], 0),
    ("Singlet: 2 electrons in 3 orbitals", 2, 3, [1, -1], [[1, 2], [1, 3]], 0),

    # Triplets are open-shell references; no closed-shell state is prepended.
    ("Triplet: 2 electrons in 2 orbitals", 2, 2, [1, 2], [], 1),
    ("Triplet: 2 electrons in 3 orbitals", 2, 3, [1, 2], [[2, 3]], 1),
]


for multiplet in MULTIPLETS:
    label, nelec, norb, reference_det, excitations, S = multiplet
    for Ms in spin_projections(S):
        run_component(label, nelec, norb, reference_det, excitations, S, Ms)
