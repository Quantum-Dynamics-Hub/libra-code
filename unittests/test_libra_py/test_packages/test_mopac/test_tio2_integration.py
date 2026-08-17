"""Numerical MOPAC regression tests for neutral and anionic TiO2."""

import itertools
import os
from pathlib import Path

import numpy as np
import pytest
import util


REPO_ROOT = Path(__file__).resolve().parents[4]
BUILD_UTIL = REPO_ROOT / "_build" / "src" / "util"
if BUILD_UTIL.is_dir():
    util.__path__.append(str(BUILD_UTIL))

from liblibra_core import Py2Cpp_int  # noqa: E402
from libra_py.packages.cp2k import methods as cp2k  # noqa: E402
from libra_py.packages.mopac import methods as mopac  # noqa: E402
from libra_py import units  # noqa: E402


MOPAC_EXE = Path(os.environ.get(
    "MOPAC_EXE", "/home/alexvakimov/SOFTWARE/mopac/_build/mopac"
))
EXAMPLE_ROOT = REPO_ROOT / "examples" / "libra_py" / "packages" / "mopac"


REFERENCES = {
    "neutral_tio2": {
        "nstates": 5,
        "multiplicity": 1,
        "spin_projection": 0,
        "nelec_act_space": 6,
        "run_params": (
            "INDO C.I.=(6,3) CHARGE=0 RELSCF=0.000001 "
            "ALLVEC WRTCONF=0.00 WRTCI=5"
        ),
        "energies": np.array([
            [0.0, 0.16222851, 0.18798648, 0.19376355, 0.20229319],
            [0.0, 0.16201904, 0.18747933, 0.19335563, 0.20143324],
            [0.0, 0.16183161, 0.18710080, 0.19303590, 0.20065415],
            [0.0, 0.16175076, 0.18688398, 0.19284848, 0.20011392],
            [0.0, 0.16188306, 0.18701261, 0.19294036, 0.20008452],
        ]),
        "overlaps": np.array([
            [[1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [0, 0, 1, 0, 0],
             [0, 0, 0, 1, 0], [0, 0, 0, 0, 1]],
            [[.99794567, .00084375332, 0, .00075930556, 0],
             [-.00084912071, .99796231, 0, .0023739378, 0],
             [0, 0, -.99794420, 0, .00010283060],
             [-.00075865193, -.0023836136, 0, .99792351, 0],
             [0, 0, -.000099759003, 0, -.99877980]],
            [[.99765734, -.00086388830, 0, .00089736066, 0],
             [-.00086818871, -.99767564, 0, .0014178572, 0],
             [0, 0, -.99765485, 0, -.00016617806],
             [-.00089688867, .0014291709, 0, .99763215, 0],
             [0, 0, -.00016493670, 0, .99861178]],
            [[.99886999, -.00068526781, 0, .00063708613, 0],
             [.00068686604, .99888193, 0, -.00083007155, 0],
             [0, 0, -.99887065, 0, .00014006627],
             [-.00063681549, .00083597676, 0, .99885946, 0],
             [0, 0, .00013931512, 0, .99931851]],
            [[.99968783, .00037229400, 0, .00035532963, 0],
             [.00037182454, -.99969550, 0, .00032901486, 0],
             [0, 0, .99969243, 0, -.000062855838],
             [-.00035544264, .00032702680, 0, .99969051, 0],
             [0, 0, .000062483783, 0, .99980277]],
        ]),
    },
    "anion_tio2": {
        "nstates": 4,
        "multiplicity": 2,
        "spin_projection": 0.5,
        "nelec_act_space": 7,
        "run_params": (
            "INDO C.I.=(6,3) CHARGE=-1 RELSCF=0.000001 "
            "ALLVEC WRTCONF=0.10 WRTCI=4"
        ),
        "energies": np.array([
            [0.0, .041604498, .041711073, .15741428],
            [0.0, .041332549, .041431774, .15724156],
            [0.0, .041167175, .041259050, .15708721],
            [0.0, .041056926, .041148800, .15701738],
            [0.0, .041064275, .041159825, .15712763],
        ]),
        "overlaps": np.array([
            [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]],
            [[.99742084, -.0013211228, 0, -.00072955484],
             [.0013222566, .99741234, 0, -.00000096692052],
             [0, 0, -.99741245, 0],
             [.00073270943, -.00000097124903, 0, .99743945]],
            [[.99713108, -.00093228517, 0, .00075471387],
             [.00093220214, .99712257, 0, .00000070565100],
             [0, 0, .99712519, 0],
             [.00075747105, -.00000070820811, 0, -.99715021]],
            [[.99866787, .00044983797, 0, -.00059710023],
             [.00044988448, -.99866293, 0, -.00000026921077],
             [0, 0, .99866466, 0],
             [-.00059799388, -.00000026658227, 0, -.99868016]],
            [[.99965403, -.00052279759, 0, -.00031410303],
             [.00052294673, .99965267, 0, -.00000016559835],
             [0, 0, .99965348, 0],
             [.00031359345, -.00000016353202, 0, .99966184]],
        ]),
    },
}


def _matrix_to_numpy(matrix):
    return np.array([
        [matrix.get(i, j) for j in range(matrix.num_of_cols)]
        for i in range(matrix.num_of_rows)
    ])


def _sign_vectors(nstates):
    return [np.asarray(signs) for signs in itertools.product((-1.0, 1.0), repeat=nstates)]


def _assert_phase_equivalent_sequence(actual, reference, atol=2.0e-5):
    """Compare overlaps using one consistent state phase at each geometry."""
    nstates = reference.shape[1]
    gauges = _sign_vectors(nstates)

    # At step zero bra and ket are the same geometry, hence they use the same
    # gauge. Retain every gauge compatible with the same-time overlap.
    possible = {
        tuple(gauge)
        for gauge in gauges
        if np.allclose(actual[0], gauge[:, None] * reference[0] * gauge[None, :], atol=atol)
    }
    assert possible, "No phase gauge reproduces the initial same-time overlap"

    # For t>0, a state's phase at geometry t is shared by overlap(t-1,t)
    # and overlap(t,t+1). Dynamic programming enforces that shared gauge.
    for step in range(1, len(reference)):
        next_possible = set()
        for previous_tuple in possible:
            previous = np.asarray(previous_tuple)
            for current in gauges:
                gauged = previous[:, None] * reference[step] * current[None, :]
                if np.allclose(actual[step], gauged, atol=atol):
                    next_possible.add(tuple(current))
        assert next_possible, f"No consistent state-phase gauge at overlap step {step}"
        possible = next_possible


@pytest.mark.parametrize("case_name", ["neutral_tio2", "anion_tio2"])
def test_tio2_energies_and_time_overlaps(case_name, tmp_path):
    if not MOPAC_EXE.is_file():
        pytest.skip(f"MOPAC executable not found: {MOPAC_EXE}")

    reference = REFERENCES[case_name]
    trajectory = EXAMPLE_ROOT / case_name / "TiO2-aligned.xyz"
    coordinate_dir = tmp_path / "coordinates"
    labels, _ = cp2k.read_trajectory_xyz_file(
        str(trajectory), 0, coordinate_dir
    )
    params = {
        "atom_labels": labels,
        "exe": str(MOPAC_EXE),
        "mopac_run_params": reference["run_params"],
        "multiplicity": reference["multiplicity"],
        "spin_projection": reference["spin_projection"],
        "nelec_act_space": reference["nelec_act_space"],
        "working_directory_prefix": str(tmp_path / "mopac_wd"),
        "mopac_input_prefix": "tio2_",
        "mopac_output_prefix": "output_",
        "nstates": reference["nstates"],
        "dt": units.fs2au,
        "do_Lowdin": True,
        "is_first_time": {0: True},
        "act_state": {0: 1},
    }
    full_id = Py2Cpp_int([0, 0])
    energies = []
    overlaps = []

    for step in range(5):
        _, q = cp2k.read_trajectory_xyz_file(
            str(trajectory), step, coordinate_dir
        )
        params["timestep"] = step
        result = mopac.mopac_compute_adi(q, params, full_id)
        energies.append(np.diag(_matrix_to_numpy(result.ham_adi)).real)
        overlaps.append(_matrix_to_numpy(result.time_overlap_adi).real)

    np.testing.assert_allclose(energies, reference["energies"], atol=2.0e-6, rtol=0)
    _assert_phase_equivalent_sequence(np.asarray(overlaps), reference["overlaps"])

    assert sorted(path.name for path in coordinate_dir.glob("coord-*.xyz")) == [
        f"coord-{step}.xyz" for step in range(5)
    ]
    assert (tmp_path / "mopac_wd_itraj0").is_dir()
    assert not list(REPO_ROOT.glob("coord-*.xyz"))
    assert not list(REPO_ROOT.glob("mopac_wd*"))
    assert not list(REPO_ROOT.glob("workflow_wd*"))
