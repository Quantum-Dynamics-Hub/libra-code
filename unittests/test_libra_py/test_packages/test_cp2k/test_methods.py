from pathlib import Path

import numpy as np
import pytest
import scipy.sparse as sp
import util


# Source-tree test runs obtain the compiled extensions from the local build
# tree while importing the Python package under test from src.
_REPO_ROOT = Path(__file__).resolve().parents[4]
_BUILD_UTIL = _REPO_ROOT / "_build" / "src" / "util"
if _BUILD_UTIL.is_dir():
    util.__path__.append(str(_BUILD_UTIL))

from liblibra_core import MATRIX, Py2Cpp_int
from libra_py.packages.cp2k import methods


def _doublet_tddfpt_data():
    info = {
        "nelec": 3,
        "nocc": 2,
        "min_occ": 1,
        "max_occ": 1,
        "min_vir": 2,
        "max_vir": 3,
    }
    # The excitation enters the alpha-only SOMO.  In a doublet reference it
    # must therefore be represented as a beta excitation.
    data = [np.array([0.1]), [[[1, 2]]], [[1.0]], [["alp"]]]
    return info, data


def test_cp2k_compute_adi_passes_open_shell_reference_to_citools(
    monkeypatch, tmp_path
):
    captured = []

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(
        methods, "read_cp2k_tddfpt_log_file", lambda params: _doublet_tddfpt_data()
    )
    monkeypatch.setattr(methods.sp, "load_npz", lambda filename: sp.eye(6, format="csr"))

    def fake_overlap(mo_overlap, data1, data2, params):
        captured.append((data1, data2, dict(params)))
        return np.eye(params["nstates"])

    monkeypatch.setattr(methods.ci, "overlap", fake_overlap)
    monkeypatch.setattr(methods.ortho, "lowdin_inverse_sqrt", np.linalg.inv)

    q = MATRIX(3, 1)
    params = {
        "atom_labels": ["H"],
        "nstates": 2,
        "logfile_name": "unused.log",
        "time_overlap_filename": "unused-st.npz",
        "overlap_filename": "unused-s.npz",
        "lowest_orbital": 1,
        "multiplicity": 2,
        "spin_projection": -0.5,
        "nelec_act_space": 3,
    }

    result = methods.cp2k_compute_adi(q, params, Py2Cpp_int([0, 0]))

    assert len(captured) == 2
    for data1, data2, overlap_params in captured:
        assert overlap_params["spin"] == pytest.approx(0.5)
        assert overlap_params["spin_projection"] == pytest.approx(-0.5)
        assert overlap_params["reference_det"] == [1, -1, 2]
        assert overlap_params["nelec"] == 3
        assert overlap_params["homo_indx"] == 2
        assert overlap_params["active_space"] == [1, 2, 3]
        assert data1[1][0][0] == [-1, -2]
        assert data2[1][0][0] == [-1, -2]

    assert result.ham_adi.get(1, 1).real == pytest.approx(0.1)
    assert result.time_overlap_adi.get(0, 0).real == pytest.approx(1.0)
    assert list(tmp_path.iterdir()) == []


def test_cp2k_compute_adi_rejects_invalid_doublet_spin_projection(tmp_path):
    q = MATRIX(3, 1)
    params = {
        "atom_labels": ["H"],
        "nstates": 2,
        "logfile_name": "unused.log",
        "time_overlap_filename": "unused-st.npz",
        "overlap_filename": "unused-s.npz",
        "lowest_orbital": 1,
        "multiplicity": 2,
        "spin_projection": 0,
    }

    with pytest.raises(ValueError, match="not a component"):
        methods.cp2k_compute_adi(q, params, Py2Cpp_int([0, 0]))

