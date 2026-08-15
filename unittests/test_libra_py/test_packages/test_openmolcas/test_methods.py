from pathlib import Path

import numpy as np
import pytest

import libra_py.packages.openmolcas.methods as molcas_methods
from libra_py.packages.openmolcas.methods import (
    ci_overlap_general,
    expand_ms_projections,
    make_molcas_input,
    read_ao_overlap,
    read_alaska_properties,
    read_alaska_vectors,
    read_ci_vectors,
    read_molcas_orbital_info,
    read_rasorb,
    read_rasscf_energies,
    read_vecdet,
)


REFERENCE_DIR = Path(__file__).parent / "reference"
REFERENCE_STEM = REFERENCE_DIR / "input__timestep_0_traj_0"
CSF_REFERENCE_DIR = REFERENCE_DIR / "csf_expansions"


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
+ sqrt(     1/1     )  |222000|
printout of CI-coefficients larger than  0.00 for root  2
energy=    -262.513349
 2 ( 2:1:  1/  1)  22ud00        0.92662      0.85862
+ sqrt(     1/2     )  |22ab00|
- sqrt(     1/2     )  |22ba00|
"""
    )

    confs, coefficients = read_ci_vectors(str(output), expected_states=2)

    factor = 0.92662 / np.sqrt(2.0)
    assert confs == [
        [(2, 2, 2, 0, 0, 0)],
        [(2, 2, 1, -1, 0, 0), (2, 2, -1, 1, 0, 0)],
    ]
    assert coefficients[0] == pytest.approx([-0.96643])
    assert coefficients[1] == pytest.approx([factor, -factor])


def test_reference_ci_vectors_expand_csfs_to_determinants():
    confs, coefficients = read_ci_vectors(
        str(_reference_file(".out")),
        expected_states=3,
        h5_path=str(_reference_file(".rasscf.h5")),
    )

    assert [len(state) for state in confs] == [172, 156, 164]
    assert all(len(conf) == 6 for state in confs for conf in state)
    assert confs[0][0] == (2, 2, 2, 0, 0, 0)
    assert confs[1][0] == (2, 2, 1, -1, 0, 0)
    assert coefficients[0][0] == pytest.approx(-0.9664253832196845)
    assert sum(value * value for value in coefficients[0]) == pytest.approx(1.0)


@pytest.mark.parametrize(
    "case,nstates,nactive",
    [("singlet", 2, 2), ("doublet", 2, 3), ("triplet", 1, 2)],
)
def test_generated_casscf_vecdet_states_are_orthonormal(case, nstates, nactive):
    stem = CSF_REFERENCE_DIR / case / case
    determinants, coefficients = read_ci_vectors(
        f"{stem}.out",
        expected_states=nstates,
        h5_path=f"{stem}.rasscf.h5",
    )
    data = ([], determinants, coefficients)

    overlap = ci_overlap_general(
        data,
        data,
        np.eye(nactive),
        active_space=list(range(1, nactive + 1)),
        inactive_orbs=[],
        nstates=nstates,
        coeff_thresh=0.0,
    )

    assert np.allclose(overlap, np.eye(nstates), atol=2e-8)


def test_generated_singlet_csf_expands_to_two_spin_determinants():
    stem = CSF_REFERENCE_DIR / "singlet" / "singlet"
    determinants, coefficients = read_ci_vectors(f"{stem}.out", expected_states=2)

    assert determinants[1] == [(2, 0), (1, -1), (-1, 1), (0, 2)]
    assert coefficients[1][1:3] == pytest.approx(
        [1.0 / np.sqrt(2.0), -1.0 / np.sqrt(2.0)], abs=1e-8
    )


def test_generated_doublet_combines_repeated_determinants():
    vecdet = CSF_REFERENCE_DIR / "doublet" / "doublet.VecDet.1"
    determinants, coefficients = read_vecdet(vecdet)

    coefficient_by_det = dict(zip(determinants, coefficients))
    assert coefficient_by_det[(1, -1, 1)] == pytest.approx(-0.0913569781)
    assert coefficient_by_det[(-1, 1, 1)] == pytest.approx(0.0610093309)


def test_ci_threshold_is_applied_after_determinants_are_combined():
    vecdet = CSF_REFERENCE_DIR / "doublet" / "doublet.VecDet.1"
    determinants, coefficients = read_vecdet(vecdet, ci_threshold=0.005)

    coefficient_by_det = dict(zip(determinants, coefficients))
    assert (1, -1, 1) in coefficient_by_det
    assert (-1, 1, 1) not in coefficient_by_det


def test_triplet_ms_projections_are_distinct_and_orthonormal():
    stem = CSF_REFERENCE_DIR / "triplet" / "triplet"
    determinants, coefficients = read_ci_vectors(f"{stem}.out", expected_states=1)
    determinants, coefficients, labels = expand_ms_projections(
        determinants, coefficients, spin_multiplicity=3
    )

    overlap = ci_overlap_general(
        ([], determinants, coefficients),
        ([], determinants, coefficients),
        np.eye(2),
        active_space=[1, 2],
        inactive_orbs=[],
        nstates=3,
        coeff_thresh=0.0,
    )

    assert labels == [(0, 2), (0, 0), (0, -2)]
    assert np.allclose(overlap, np.eye(3), atol=2e-8)


def test_quartet_ms_projection_generalization():
    determinants, coefficients, labels = expand_ms_projections(
        [[(1, 1, 1)]], [[1.0]], spin_multiplicity=4
    )

    overlap = ci_overlap_general(
        ([], determinants, coefficients),
        ([], determinants, coefficients),
        np.eye(3),
        active_space=[1, 2, 3],
        inactive_orbs=[],
        nstates=4,
        coeff_thresh=0.0,
    )

    assert labels == [(0, 3), (0, 1), (0, -1), (0, -3)]
    assert np.allclose(overlap, np.eye(4), atol=2e-8)


def test_mixed_spin_manifolds_are_assembled_as_zero_coupled_blocks(monkeypatch):
    results = [
        {
            "multiplicity": 1, "nroots": 1, "nstates": 1,
            "energies": [-1.0], "ms_labels": [(0, 0)],
            "time_overlap": np.array([[0.9]]),
            "overlap": np.array([[1.0]]),
            "gradients": {}, "nac_vectors": {},
        },
        {
            "multiplicity": 3, "nroots": 1, "nstates": 3,
            "energies": [-0.8, -0.8, -0.8],
            "ms_labels": [(0, 2), (0, 0), (0, -2)],
            "time_overlap": 0.8 * np.eye(3),
            "overlap": np.eye(3),
            "gradients": {}, "nac_vectors": {},
        },
    ]
    calls = []

    def fake_compute(coords, params, itraj, manifold, index, multiple):
        calls.append((index, multiple, manifold["spin"]))
        return results[index]

    monkeypatch.setattr(molcas_methods, "Cpp2Py", lambda full_id: [0])
    monkeypatch.setattr(molcas_methods, "_compute_spin_manifold", fake_compute)
    q = type("Coordinates", (), {"col": lambda self, index: object()})()
    obj = molcas_methods.molcas_compute_adi(q, {
        "atom_labels": ["H"],
        "spin_manifolds": [{"spin": 1}, {"spin": 3}],
    }, 0)

    assert calls == [(0, True, 1), (1, True, 3)]
    assert obj.spin_labels == [(1, 0, 0), (3, 0, 2), (3, 0, 0), (3, 0, -2)]
    assert [obj.ham_adi.get(i, i).real for i in range(4)] == pytest.approx(
        [-1.0, -0.8, -0.8, -0.8]
    )
    assert all(obj.time_overlap_adi.get(0, j) == 0.0j for j in range(1, 4))
    assert all(obj.time_overlap_adi.get(j, 0) == 0.0j for j in range(1, 4))
    assert all(obj.ham_adi.get(0, j) == 0.0j for j in range(1, 4))


def test_single_spin_driver_uses_legacy_special_case(monkeypatch):
    calls = []
    result = {
        "multiplicity": 2, "nroots": 1, "nstates": 1,
        "energies": [-0.5], "ms_labels": [(0, 1)],
        "time_overlap": np.eye(1), "overlap": np.eye(1),
        "gradients": {}, "nac_vectors": {},
    }

    def fake_compute(coords, params, itraj, manifold, index, multiple):
        calls.append((manifold, index, multiple))
        return result

    monkeypatch.setattr(molcas_methods, "Cpp2Py", lambda full_id: [0])
    monkeypatch.setattr(molcas_methods, "_compute_spin_manifold", fake_compute)
    q = type("Coordinates", (), {"col": lambda self, index: object()})()
    obj = molcas_methods.molcas_compute_adi(q, {
        "atom_labels": ["H"], "molcas_run_params": {"spin": 2},
    }, 0)

    assert calls == [({}, 0, False)]
    assert obj.ms_labels == [(0, 1)]
    assert obj.spin_labels == [(2, 0, 1)]


@pytest.mark.parametrize("case,nstates", [("singlet", 2), ("doublet", 2), ("triplet", 1)])
def test_text_and_hdf5_fallback_matches_vecdet(case, nstates):
    stem = CSF_REFERENCE_DIR / case / case
    vecdet_data = read_ci_vectors(f"{stem}.out", expected_states=nstates)
    parsed_data = read_ci_vectors(
        f"{stem}.out",
        expected_states=nstates,
        h5_path=f"{stem}.rasscf.h5",
        vecdet_prefix=f"{stem}.missing",
    )

    assert parsed_data[0] == vecdet_data[0]
    for parsed, vecdet in zip(parsed_data[1], vecdet_data[1]):
        assert parsed == pytest.approx(vecdet, abs=2e-8)


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


def test_make_molcas_input_uses_uhf_for_open_shell(tmp_path):
    input_file = tmp_path / "open_shell.in"
    params = {
        "charge": 0,
        "spin": 2,
        "nactel": "3 0 0",
        "inactive": 18,
        "ras2": 6,
        "ciroot": "3 3 1",
        "scf_method": "auto",
        "uhf_orbital_set": "beta",
    }

    make_molcas_input(input_file, params, ["Al", "Al", "Al"], {})
    text = input_file.read_text()

    assert text.index("&SCF\nUHF\n") < text.index("Spin=2")
    assert "FILEORB=$Project.UhfOrb" in text
    assert "AlphaOrBeta=-1" in text


def test_make_molcas_input_rohf_uses_spin_adapted_rasscf(tmp_path):
    input_file = tmp_path / "rohf.in"
    params = {"spin": 2, "nactel": 3, "scf_method": "rohf"}

    make_molcas_input(input_file, params, ["H"], {})
    text = input_file.read_text()

    assert "&SCF" not in text
    assert "&RASSCF" in text
    assert "FILEORB=$Project.GssOrb" in text


def test_make_molcas_input_rejects_incompatible_spin_and_electron_parity(tmp_path):
    with pytest.raises(ValueError, match="incompatible electron parity"):
        make_molcas_input(
            tmp_path / "invalid.in",
            {"spin": 2, "nactel": "6 0 0", "scf_method": "uhf"},
            ["Al", "Al", "Al"],
            {},
        )


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
