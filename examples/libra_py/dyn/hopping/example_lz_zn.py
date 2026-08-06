"""Landau–Zener and Zhu–Nakamura crossing examples."""

import numpy as np

from libra_py.dyn.hopping import hopping_probabilities_lz, hopping_probabilities_zn


def show(label, values):
    print(f"{label:35s} {np.asarray(values)}")


def main():
    previous = {"ham_dia": np.diag([-0.1, 0.1])}
    current_dia = {
        "ham_dia": np.array([[0.1, 0.01], [0.01, -0.1]]),
        "d1ham_dia": np.array([[[1.0, 0.0], [0.0, -1.0]]]),
    }
    show(
        "LZ diabatic, velocity 1",
        hopping_probabilities_lz(current_dia, previous, 0, 0, [1.0], [1.0]),
    )
    show(
        "LZ diabatic, velocity 4",
        hopping_probabilities_lz(current_dia, previous, 0, 0, [4.0], [1.0]),
    )

    no_crossing = {
        "ham_dia": np.array([[0.2, 0.01], [0.01, 0.0]]),
        "d1ham_dia": current_dia["d1ham_dia"],
    }
    show(
        "LZ without a crossing",
        hopping_probabilities_lz(
            no_crossing, {"ham_dia": np.diag([0.1, 0.0])}, 0, 0, [1.0], [1.0]
        ),
    )

    current_adi = {
        "ham_dia": np.diag([0.1, -0.1]),
        "ham_adi": np.diag([0.0, 0.02]),
        "nac_adi": np.array([[0.0, 0.1], [-0.1, 0.0]]),
    }
    previous_adi = {
        "ham_dia": np.diag([-0.1, 0.1]),
        "nac_adi": np.array([[0.0, -0.1], [0.1, 0.0]]),
    }
    show(
        "LZ adiabatic, gap crossing",
        hopping_probabilities_lz(current_adi, previous_adi, 0, 1, [1.0], [1.0]),
    )
    show(
        "LZ adiabatic, NAC crossing",
        hopping_probabilities_lz(current_adi, previous_adi, 0, 2, [1.0], [1.0]),
    )

    current_zn = {
        "ham_dia": np.array([[0.1, 0.1], [0.1, -0.1]]),
        "forces_adi": np.array([[1.0, -1.0]]),
    }
    show(
        "ZN, opposing forces",
        hopping_probabilities_zn(current_zn, previous, 0, 0, [0.0], [1.0]),
    )

    current_zn_2d = {
        "ham_dia": np.array([[0.1, 0.08], [0.08, -0.1]]),
        "forces_adi": np.array([[1.0, -0.5], [0.4, -1.2]]),
    }
    show(
        "ZN, two nuclear DOFs",
        hopping_probabilities_zn(
            current_zn_2d, previous, 0, 0, [0.0, 0.0], [1.0, 0.5]
        ),
    )


if __name__ == "__main__":
    main()
