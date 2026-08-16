import math

import pytest

from liblibra_core import MATRIX, nac_npi


def rotation(angle):
    overlap = MATRIX(2, 2)
    overlap.set(0, 0, math.cos(angle))
    overlap.set(0, 1, -math.sin(angle))
    overlap.set(1, 0, math.sin(angle))
    overlap.set(1, 1, math.cos(angle))
    return overlap


def three_state_rotation(angle):
    overlap = MATRIX(3, 3)
    overlap.set(0, 0, math.cos(angle))
    overlap.set(0, 1, -math.sin(angle))
    overlap.set(1, 0, math.sin(angle))
    overlap.set(1, 1, math.cos(angle))
    overlap.set(2, 2, 1.0)
    return overlap


def test_npi_removable_singularity_has_finite_two_state_limit():
    angle = 0.4572208409065037
    # These values come from exponentiating an antisymmetric generator. Their
    # last-bit asymmetry makes the old E expression evaluate as 0/0.
    overlap = MATRIX(2, 2)
    overlap.set(0, 0, 0.8972828379263506)
    overlap.set(0, 1, -0.44145612325896505)
    overlap.set(1, 0, 0.441456123258965)
    overlap.set(1, 1, 0.8972828379263506)
    nac = nac_npi(overlap, 1.0)

    assert math.isfinite(nac.get(1, 0))
    assert nac.get(1, 0) == pytest.approx(angle, abs=1.0e-12)
    assert nac.get(0, 1) == pytest.approx(-angle, abs=1.0e-12)


def test_npi_rejects_non_square_overlap():
    with pytest.raises(ValueError, match="non-empty square matrix"):
        nac_npi(MATRIX(2, 3), 1.0)


def test_npi_rejects_non_positive_timestep():
    with pytest.raises(ValueError, match="dt must be finite and positive"):
        nac_npi(rotation(0.1), 0.0)


def test_npi_rejects_non_finite_overlap():
    overlap = rotation(0.1)
    overlap.set(0, 0, float("nan"))
    with pytest.raises(ValueError, match="overlap element .* is not finite"):
        nac_npi(overlap, 1.0)


def test_npi_rejects_unmatched_state_phase():
    overlap = MATRIX(2, 2)
    overlap.set(0, 0, -1.0)
    overlap.set(1, 1, -1.0)
    with pytest.raises(ValueError, match="discontinuous electronic-state phase"):
        nac_npi(overlap, 1.0)


def test_npi_rejects_non_orthogonal_overlap():
    overlap = MATRIX(2, 2)
    overlap.set(0, 0, 0.9)
    overlap.set(1, 1, 1.0)
    with pytest.raises(ValueError, match="not orthogonal"):
        nac_npi(overlap, 1.0)


def test_npi_rejects_improper_rotation():
    overlap = MATRIX(2, 2)
    overlap.set(0, 1, 1.0)
    overlap.set(1, 0, 1.0)
    with pytest.raises(ValueError, match="proper rotation"):
        nac_npi(overlap, 1.0)


def test_npi_accepts_valid_three_state_rotation():
    angle = 0.2
    dt = 0.5
    nac = nac_npi(three_state_rotation(angle), dt)

    assert nac.get(1, 0) == pytest.approx(angle / dt, abs=1.0e-12)
    assert nac.get(0, 1) == pytest.approx(-angle / dt, abs=1.0e-12)
    for i, j in ((0, 2), (2, 0), (1, 2), (2, 1)):
        assert nac.get(i, j) == pytest.approx(0.0, abs=1.0e-12)


def test_npi_rejects_non_orthogonal_three_state_overlap():
    overlap = three_state_rotation(0.2)
    overlap.set(2, 2, 0.95)
    with pytest.raises(ValueError, match=r"max\|S\^T S - I\|"):
        nac_npi(overlap, 1.0)


def test_npi_rejects_unmatched_phase_in_three_state_overlap():
    overlap = MATRIX(3, 3)
    overlap.set(0, 0, -1.0)
    overlap.set(1, 1, -1.0)
    overlap.set(2, 2, 1.0)
    with pytest.raises(ValueError, match=r"diagonal overlap S\(0, 0\).+negative"):
        nac_npi(overlap, 1.0)


def test_npi_rejects_improper_three_state_rotation():
    overlap = MATRIX(3, 3)
    overlap.set(0, 1, 1.0)
    overlap.set(1, 0, 1.0)
    overlap.set(2, 2, 1.0)
    with pytest.raises(ValueError, match=r"det\(S\).+positive"):
        nac_npi(overlap, 1.0)
