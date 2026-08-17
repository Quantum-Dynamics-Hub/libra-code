from pathlib import Path

import numpy as np
import pytest
import util


# Source-tree test runs obtain the compiled extensions from the local build
# tree while importing the Python package under test from src.
_REPO_ROOT = Path(__file__).resolve().parents[4]
_BUILD_UTIL = _REPO_ROOT / "_build" / "src" / "util"
if _BUILD_UTIL.is_dir():
    util.__path__.append(str(_BUILD_UTIL))

from libra_py.packages.dftbplus import methods


def _mock_dftb_readers(monkeypatch, vector):
    # alpha 2->3, unlabelled/restricted 1->2, beta 2->3
    lookup = np.array([[1, 2, 0], [0, 1, 0], [1, 2, 1]], dtype=int)
    monkeypatch.setattr(
        methods,
        "read_spx_mappings",
        lambda filename: (np.empty((3, 3, 2), dtype=int), lookup, 1),
    )
    monkeypatch.setattr(
        methods,
        "read_mo_matrix",
        lambda filename, ndim, spin_polarized=False: np.array([np.eye(ndim)]),
    )
    monkeypatch.setattr(
        methods,
        "read_xplusy_ascii",
        lambda filename: (3, 1, [{"vector": vector, "energy": 0.1}]),
    )


def test_read_dftb_orbital_info_builds_signed_doublet_configurations(
    monkeypatch, tmp_path
):
    monkeypatch.chdir(tmp_path)
    _mock_dftb_readers(monkeypatch, [1.0, 0.5, 0.25])

    info, mos, data = methods.read_dftb_orbital_info(
        {
            "source_directory": "unused",
            "nstates": 2,
            "ci_threshold": 0.0,
            "multiplicity": 2,
            "spin_projection": -0.5,
        }
    )

    assert info["spin"] == pytest.approx(0.5)
    assert info["spin_projection"] == pytest.approx(-0.5)
    assert info["reference_det"] == [1, -1, 2]
    assert info["nelec"] == 3
    assert info["nocc"] == 2
    assert np.array_equal(mos, np.eye(3))
    assert data[1][0].tolist() == [[2, 3], [-1, -2], [-2, -3]]
    assert data[2][0].tolist() == pytest.approx([1.0, 0.5, 0.25])
    assert list(tmp_path.iterdir()) == []


def test_dftb_reference_does_not_depend_on_ci_truncation(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    _mock_dftb_readers(monkeypatch, [0.0, 1.0, 0.0])

    info, _, data = methods.read_dftb_orbital_info(
        {
            "source_directory": "unused",
            "nstates": 2,
            "ci_threshold": 0.5,
            "multiplicity": 2,
            "spin_projection": 0.5,
        }
    )

    # Only 1->2 survives, but the full SPX map still establishes orbital 2 as
    # the SOMO and three electrons as the doublet reference.
    assert data[1][0].tolist() == [[-1, -2]]
    assert info["nocc"] == 2
    assert info["nelec"] == 3
    assert info["reference_det"] == [1, -1, 2]
    assert list(tmp_path.iterdir()) == []

