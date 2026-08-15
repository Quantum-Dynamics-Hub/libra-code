from pathlib import Path

import numpy as np
import pytest

from libra_py.packages.openmolcas.methods import (
    make_molcas_input,
    read_ao_overlap,
    read_alaska_properties,
    read_alaska_vectors,
    read_ci_vectors,
    read_molcas_orbital_info,
    read_rasorb,
    read_rasscf_energies,
)


REFERENCE_DIR = Path(__file__).parent / "reference"
REFERENCE_STEM = REFERENCE_DIR / "input__timestep_0_traj_0"


def _reference_file(suffix):
    return Path(f"{REFERENCE_STEM}{suffix}")


def test_read_ci_vectors_does_not_parse_sguga_indices_as_occupations(tmp_path):
    output = tmp_path / "molcas.out"
    output.write_text(
        """
printout of CI-coefficients larger than  0.00 for root  1
energy=    -262.743158
Conf   SGUGA info  Occupation       Coef       Weight
 1 ( 1:1:  1/  1)  222000       -0.96643      0.93398
printout of CI-coefficients larger than  0.00 for root  2
energy=    -262.513349
 2 ( 2:1:  1/  1)  22ud00        0.92662      0.85862
"""
    )

    confs, coefficients = read_ci_vectors(str(output), expected_states=2)

    assert confs == [[(2, 2, 2, 0, 0, 0)], [(2, 2, 1, -1, 0, 0)]]
    assert coefficients == [[-0.96643], [0.92662]]


def test_reference_ci_vectors():
    confs, coefficients = read_ci_vectors(str(_reference_file(".out")), expected_states=3)

    assert [len(state) for state in confs] == [80, 65, 67]
    assert all(len(conf) == 6 for state in confs for conf in state)
    assert confs[0][0] == (2, 2, 2, 0, 0, 0)
    assert confs[1][0] == (2, 2, 1, -1, 0, 0)
    assert coefficients[0][0] == pytest.approx(-0.96643)


def test_reference_rasscf_energies():
    energies = read_rasscf_energies(str(_reference_file(".rasscf.h5")))

    assert energies == pytest.approx(
        [-262.7431577786203, -262.5133486129196, -262.4794630257731]
    )


def test_reference_rasorb_coefficients():
    mos = read_rasorb(str(_reference_file(".RasOrb")))

    assert mos.shape == (104, 104)
    assert mos[0, 0] == pytest.approx(-0.708242867162235)
    assert mos[-1, -1] == pytest.approx(-6.45024439706388e-13)
    assert np.linalg.norm(mos) == pytest.approx(48.2002240142368)


def test_reference_ao_overlap():
    overlap = read_ao_overlap(str(_reference_file(".out")), 104)

    assert overlap.shape == (104, 104)
    assert np.allclose(overlap, overlap.T)
    assert np.allclose(np.diag(overlap), 1.0, atol=1e-3)


def test_reference_orbital_info():
    info, mos, data = read_molcas_orbital_info(
        {
            "filename": str(_reference_file(".out")),
            "rasorb_file": str(_reference_file(".RasOrb")),
            "h5_path": str(_reference_file(".rasscf.h5")),
        }
    )

    assert info["nao"] == 104
    assert info["nocc"] == 18
    assert info["nact"] == 6
    assert info["actual_orbital_space"] == [19, 20, 21, 22, 23, 24]
    assert mos.shape == (104, 104)
    assert len(data["energies"]) == len(data["confs"]) == len(data["CI"]) == 3


def test_make_molcas_input_requests_gradients_and_nacs(tmp_path):
    input_file = tmp_path / "molcas.in"
    params = {
        "ciroot": "3 3 1",
        "nstates": 3,
        "gradient_states": "all",
        "nac_pairs": [(0, 2)],
        "nac_nocsf": True,
    }

    make_molcas_input(input_file, params, ["H", "H"], {})
    text = input_file.read_text()

    assert text.count("&ALASKA") == 4
    assert "ROOT=1" in text
    assert "ROOT=2" in text
    assert "ROOT=3" in text
    assert "NAC=1 3\nNOCSF" in text


def test_read_alaska_gradients_and_nac_vectors(tmp_path):
    output = tmp_path / "molcas.out"
    output.write_text(
        """
                       Molecular gradients
 Irreducible representation: a
 ------------------------------------------------------------------------------------------
                              X                      Y                      Z
 ------------------------------------------------------------------------------------------
  N1      1.00000000000000D-03   2.00000000000000D-03   3.00000000000000D-03
  H2     -4.00000000000000D-03   5.00000000000000D-03  -6.00000000000000D-03
 ------------------------------------------------------------------------------------------
                       Total derivative coupling
 Irreducible representation: a
 ------------------------------------------------------------------------------------------
                              X                      Y                      Z
 ------------------------------------------------------------------------------------------
  N1      1.00000000000000D-02  -2.00000000000000D-02   3.00000000000000D-02
  H2      4.00000000000000D-02   5.00000000000000D-02  -6.00000000000000D-02
 ------------------------------------------------------------------------------------------
"""
    )

    blocks = read_alaska_vectors(output)
    gradients, nacs = read_alaska_properties(
        output, natoms=2, gradient_states=[1], nac_pairs=[(0, 2)], nstates=3
    )

    assert [kind for kind, _ in blocks] == ["gradient", "nac"]
    assert np.allclose(gradients[1][0], [1.0e-3, 2.0e-3, 3.0e-3])
    assert np.allclose(nacs[(0, 2)][1], [4.0e-2, 5.0e-2, -6.0e-2])


def test_read_alaska_properties_checks_requested_table_count(tmp_path):
    output = tmp_path / "molcas.out"
    output.write_text("OpenMolcas completed without an ALASKA request\n")

    with pytest.raises(ValueError, match="0 ALASKA gradient tables, expected 1"):
        read_alaska_properties(output, natoms=2, gradient_states=[0], nstates=2)
