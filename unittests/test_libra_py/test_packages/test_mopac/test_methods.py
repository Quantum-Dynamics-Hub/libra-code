from pathlib import Path

import pytest
import util

# Source-tree test runs obtain the compiled ``util.libutil`` extension from
# the local build tree while importing the Python package under test from src.
_BUILD_UTIL = Path(__file__).resolve().parents[4] / "_build" / "src" / "util"
if _BUILD_UTIL.is_dir():
    util.__path__.append(str(_BUILD_UTIL))

from libra_py.packages.mopac import methods


@pytest.mark.parametrize(
    "n_doubly,multiplicity,expected",
    [
        (0, 2, [1]),
        (1, 2, [1, -1, 2]),
        (1, 4, [1, -1, 2, 3, 4]),
    ],
)
def test_make_open_shell_reference(n_doubly, multiplicity, expected):
    assert methods.make_open_shell_reference(n_doubly, multiplicity) == expected


@pytest.mark.parametrize(
    "multiplicity,components",
    [
        (1, [0]),
        (2, [-0.5, 0.5]),
        (3, [-1, 0, 1]),
        (4, [-1.5, -0.5, 0.5, 1.5]),
    ],
)
def test_spin_quantum_numbers_accepts_every_ms_component(multiplicity, components):
    for component in components:
        spin, projection = methods.spin_quantum_numbers(multiplicity, component)
        assert spin == 0.5 * (multiplicity - 1)
        assert projection == component


def test_spin_quantum_numbers_rejects_wrong_projection_parity():
    with pytest.raises(ValueError, match="not a component"):
        methods.spin_quantum_numbers(2, 0)


def test_add_mopac_spin_keyword_preserves_matching_explicit_selection():
    base = "INDO C.I.=(6,3)"
    assert methods.add_mopac_spin_keyword(base, 4).endswith("QUARTET")
    assert methods.add_mopac_spin_keyword(base + " QUARTET", 4).endswith("QUARTET")


def test_add_mopac_spin_keyword_rejects_conflicting_selection():
    with pytest.raises(ValueError, match="selects multiplicity 2"):
        methods.add_mopac_spin_keyword("INDO C.I.=(6,3) DOUBLET", 4)
