import math

import numpy as np
import pytest

from liblibra_core import MATRIX, nac_npi as cpp_nac_npi
from libra_py.dyn.algorithms.npi import nac_npi
from libra_py.dyn.core.storage import TensorStorage


def rotation2(angle):
    c, s = math.cos(angle), math.sin(angle)
    return np.array([[c, -s], [s, c]], dtype=float)


def rotation3(angle01, angle12=0.0):
    c01, s01 = math.cos(angle01), math.sin(angle01)
    c12, s12 = math.cos(angle12), math.sin(angle12)
    r01 = np.array([[c01, -s01, 0.0], [s01, c01, 0.0], [0.0, 0.0, 1.0]])
    r12 = np.array([[1.0, 0.0, 0.0], [0.0, c12, -s12], [0.0, s12, c12]])
    return r01 @ r12


def to_cpp(array):
    matrix = MATRIX(*array.shape)
    for i in range(array.shape[0]):
        for j in range(array.shape[1]):
            matrix.set(i, j, float(array[i, j]))
    return matrix


def from_cpp(matrix, size):
    return np.array(
        [[matrix.get(i, j) for j in range(size)] for i in range(size)]
    )


def test_numpy_two_state_removable_singularity():
    overlap = np.array(
        [
            [0.8972828379263506, -0.44145612325896505],
            [0.441456123258965, 0.8972828379263506],
        ]
    )
    result = nac_npi(overlap, 1.0)

    assert np.all(np.isfinite(result))
    assert result[1, 0] == pytest.approx(0.4572208409065037, abs=1.0e-12)


def test_numpy_accepts_tensor_storage_batch_layout_and_complex_dtype():
    storage = TensorStorage(
        backend=np, ntraj=2, ndof=1, nstates=3, ntbf_initial=1, ntbf_capacity=3
    )
    overlaps = storage.time_overlap_adi
    for trajectory in range(2):
        for tbf in range(3):
            overlaps[trajectory, tbf] = rotation3(0.05 * (trajectory + tbf + 1))

    result = nac_npi(overlaps, 0.5)

    assert result.shape == overlaps.shape
    assert result.dtype == overlaps.dtype
    np.testing.assert_allclose(result + result.swapaxes(-1, -2), 0.0, atol=1.0e-12)
    np.testing.assert_allclose(result.imag, 0.0)


def test_numpy_validates_every_matrix_in_a_batch():
    overlaps = np.stack([rotation3(0.1), rotation3(0.2)])
    overlaps[1, 2, 2] = 0.9
    with pytest.raises(ValueError, match="not orthogonal"):
        nac_npi(overlaps, 1.0)


def test_complex_overlap_requires_real_gauge():
    overlap = rotation2(0.1).astype(complex)
    overlap[0, 1] += 1.0e-3j
    with pytest.raises(ValueError, match=r"requires real overlaps.+max\|Im\(S\)\|"):
        nac_npi(overlap, 1.0)


@pytest.mark.parametrize(
    "overlap, dt",
    [
        (rotation2(0.31), 0.7),
        (rotation3(0.18), 0.4),
        (rotation3(0.12, -0.08), 0.25),
    ],
)
def test_numpy_matches_cpp(overlap, dt):
    expected = from_cpp(cpp_nac_npi(to_cpp(overlap), dt), overlap.shape[0])
    np.testing.assert_allclose(nac_npi(overlap, dt), expected, atol=2.0e-12, rtol=0.0)


def test_torch_preserves_backend_dtype_device_and_batch_shape():
    torch = pytest.importorskip("torch")
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    storage = TensorStorage(
        backend=torch,
        ntraj=2,
        ndof=1,
        nstates=3,
        ntbf_initial=1,
        ntbf_capacity=1,
        device=str(device),
    )
    overlaps = storage.time_overlap_adi.to(device=device)
    overlaps[:] = torch.as_tensor(
        np.stack([rotation3(0.1), rotation3(0.2)])[:, None],
        dtype=overlaps.dtype,
        device=device,
    )

    result = nac_npi(overlaps, 0.5)

    assert isinstance(result, torch.Tensor)
    assert result.shape == overlaps.shape
    assert result.dtype == overlaps.dtype
    assert result.device == overlaps.device
    np.testing.assert_allclose(
        (result + result.transpose(-1, -2)).cpu().numpy(), 0.0, atol=1.0e-12
    )


def test_torch_matches_numpy_and_cpp():
    torch = pytest.importorskip("torch")
    overlap = rotation3(0.12, -0.08)
    torch_result = nac_npi(torch.as_tensor(overlap, dtype=torch.float64), 0.25)
    numpy_result = nac_npi(overlap, 0.25)
    cpp_result = from_cpp(cpp_nac_npi(to_cpp(overlap), 0.25), overlap.shape[0])

    np.testing.assert_allclose(torch_result.cpu().numpy(), numpy_result, atol=2.0e-12)
    np.testing.assert_allclose(torch_result.cpu().numpy(), cpp_result, atol=2.0e-12)
