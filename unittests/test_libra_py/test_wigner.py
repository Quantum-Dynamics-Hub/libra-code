from pathlib import Path

import numpy as np
import pytest

from libra_py import units
from libra_py.wigner import (
    build_mass_weighted_eigenvectors,
    build_modes_from_hessian,
    generate_wigner_from_hessian,
    generate_wigner_ics,
    prepare_wigner_from_modes,
    read_modes_i_range,
)


REFERENCE_DIR = (
    Path(__file__).resolve().parents[2]
    / "examples" / "libra_py" / "wigner_reference"
)
WATER_LABELS = ["O", "H", "H"]
WATER_MASSES = {"O": 15.999, "H": 1.00784}
WATER_FREQUENCIES_CM = np.array([1594.75, 3657.05, 3755.93])
WATER_Q_EQ_BOHR = np.array([
    0.0, 0.0, 0.0,
    0.75716, 0.58626, 0.0,
    -0.75716, 0.58626, 0.0,
]) * units.Angst


def test_read_realistic_water_mode_files():
    natoms, frequencies, displacements = read_modes_i_range(
        1, 3, str(REFERENCE_DIR / "mode_{}.xyz")
    )

    assert natoms == 3
    np.testing.assert_allclose(frequencies, WATER_FREQUENCIES_CM)
    assert displacements.shape == (9, 3)
    assert np.all(np.linalg.norm(displacements, axis=0) > 1.0)


def test_water_mode_transforms_are_mass_orthonormal():
    _, _, displacements = read_modes_i_range(
        1, 3, str(REFERENCE_DIR / "mode_{}.xyz")
    )
    d_cart, d_p, mass_cart = build_mass_weighted_eigenvectors(
        displacements, WATER_LABELS, WATER_MASSES
    )

    np.testing.assert_allclose(d_cart.T @ d_p, np.eye(3), atol=2.0e-10)
    np.testing.assert_allclose(
        mass_cart[::3],
        np.array([15.999, 1.00784, 1.00784]) * 1822.888,
        rtol=1.0e-12,
    )


def test_prepare_water_from_mode_files_is_reproducible():
    kwargs = dict(
        labels=WATER_LABELS,
        q_eq=WATER_Q_EQ_BOHR,
        mode_start=1,
        mode_end=3,
        mode_file_pattern=str(REFERENCE_DIR / "mode_{}.xyz"),
        mass_map=WATER_MASSES,
        temperature=300.0,
        ntraj=4,
        seed=81,
    )
    first = prepare_wigner_from_modes(**kwargs)
    second = prepare_wigner_from_modes(**kwargs)

    np.testing.assert_allclose(first["freqs_cm"], WATER_FREQUENCIES_CM)
    for ic_first, ic_second in zip(first["ics"], second["ics"]):
        np.testing.assert_array_equal(ic_first["q"], ic_second["q"])
        np.testing.assert_array_equal(ic_first["p"], ic_second["p"])


def test_build_modes_from_diagonal_hessian():
    masses_au = np.array([2.0, 8.0, 18.0])
    expected_omega = np.array([0.1, 0.2, 0.3])
    hessian = np.diag(masses_au * expected_omega**2)

    omega, d_cart, d_p, mass_cart = build_modes_from_hessian(
        hessian, masses_au, amu_to_au=1.0
    )

    np.testing.assert_allclose(omega, expected_omega)
    np.testing.assert_allclose(mass_cart, masses_au)
    np.testing.assert_allclose(d_cart.T @ d_p, np.eye(3), atol=1.0e-14)
    np.testing.assert_allclose(
        d_cart.T @ hessian @ d_cart, np.diag(expected_omega**2),
        atol=1.0e-14
    )


def test_generate_from_hessian_has_ground_state_variances():
    mass = 4.0
    omega = 0.25
    result = generate_wigner_from_hessian(
        np.zeros(1), np.array([[mass * omega**2]]), [mass], 0.0,
        ntraj=120000, seed=19, amu_to_au=1.0
    )
    q = np.array([ic["q"][0] for ic in result["ics"]])
    p = np.array([ic["p"][0] for ic in result["ics"]])

    assert np.var(q) == pytest.approx(1.0 / (2.0 * mass * omega), rel=0.02)
    assert np.var(p) == pytest.approx(mass * omega / 2.0, rel=0.02)


def test_realistic_water_hessian_recovers_reference_frequencies():
    hessian = np.loadtxt(REFERENCE_DIR / "water_hessian_au.txt")
    omega, d_cart, d_p, mass_cart = build_modes_from_hessian(
        hessian, [15.999, 1.00784, 1.00784]
    )

    assert hessian.shape == (9, 9)
    np.testing.assert_allclose(hessian, hessian.T, atol=1.0e-14)
    np.testing.assert_allclose(
        omega[omega > 0.0] / units.inv_cm2Ha,
        WATER_FREQUENCIES_CM,
        atol=0.06,
    )
    np.testing.assert_allclose(d_cart.T @ d_p, np.eye(9), atol=1.0e-13)
    np.testing.assert_allclose(
        d_cart.T @ hessian @ d_cart, np.diag(omega**2), atol=1.0e-14
    )
    assert mass_cart.shape == (9,)


def test_water_hessian_sampling_has_expected_thermal_modal_variances():
    hessian = np.loadtxt(REFERENCE_DIR / "water_hessian_au.txt")
    result = generate_wigner_from_hessian(
        WATER_Q_EQ_BOHR, hessian, [15.999, 1.00784, 1.00784], 300.0,
        ntraj=40000, seed=2026
    )
    active = result["omega_au"] > 0.0
    omega = result["omega_au"][active]
    d_p = result["D_p"][:, active]
    d_cart = result["D_cart"][:, active]
    q = np.stack([ic["q"] for ic in result["ics"]]) - WATER_Q_EQ_BOHR
    p = np.stack([ic["p"] for ic in result["ics"]])
    q_normal = q @ d_p
    p_normal = p @ d_cart
    coth = 1.0 / np.tanh(omega / (2.0 * units.kB * 300.0))

    np.testing.assert_allclose(
        np.var(q_normal, axis=0), coth / (2.0 * omega), rtol=0.025
    )
    np.testing.assert_allclose(
        np.var(p_normal, axis=0), omega * coth / 2.0, rtol=0.025
    )


def test_generate_wigner_ics_zero_temperature_seed_and_shape():
    d_cart = np.eye(2)
    d_p = np.eye(2)
    first = generate_wigner_ics(
        np.array([2.0, -1.0]), d_cart, d_p, [0.2, 0.4], 0.0,
        ntraj=2, seed=5
    )
    second = generate_wigner_ics(
        np.array([2.0, -1.0]), d_cart, d_p, [0.2, 0.4], 0.0,
        ntraj=2, seed=5
    )

    assert [ic["traj"] for ic in first] == [0, 1]
    assert all(ic["q"].shape == (2,) and ic["p"].shape == (2,) for ic in first)
    np.testing.assert_array_equal(first[0]["q"], second[0]["q"])
    np.testing.assert_array_equal(first[0]["p"], second[0]["p"])


def test_hessian_zero_mode_is_frozen_and_imaginary_mode_is_rejected():
    result = generate_wigner_from_hessian(
        np.array([1.5]), np.zeros((1, 1)), [2.0], 300.0,
        ntraj=3, seed=7, amu_to_au=1.0
    )
    assert all(ic["q"][0] == 1.5 and ic["p"][0] == 0.0 for ic in result["ics"])

    with pytest.raises(ValueError, match="imaginary mode"):
        build_modes_from_hessian(np.array([[-0.1]]), [1.0], amu_to_au=1.0)


@pytest.mark.parametrize(
    "hessian,masses,message",
    [
        (np.zeros((2, 3)), [1.0, 1.0], "square"),
        (np.array([[1.0, 0.1], [0.0, 1.0]]), [1.0, 1.0], "symmetric"),
        (np.eye(2), [1.0, -1.0], "positive"),
    ],
)
def test_build_modes_from_hessian_rejects_invalid_inputs(hessian, masses, message):
    with pytest.raises(ValueError, match=message):
        build_modes_from_hessian(hessian, masses, amu_to_au=1.0)
