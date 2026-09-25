# Small overlap checks

Checked with PySCF 2.13.1 using the repository's Python implementation and one thread. The PySCF interface unit-test directory passed all 14 tests, including six focused TDDFT/TDA determinant-overlap tests. The combined run emitted existing duplicate Boost/Python converter registration warnings.

All ten independent-geometry examples converged. The table lists absolute matrix elements because independently chosen electronic signs can differ. Geometry: Li–H 3.0 → 3.1 Bohr; basis: STO-3G; two roots. These are illustrative method-dependent overlaps, not accuracy benchmarks or comparable spectra.

| System/method | abs(U00) | abs(U01) | abs(U10) | abs(U11) | Current self-overlap error | MO metric error |
|---|---:|---:|---:|---:|---:|---:|
| lih_casscf_hf | 0.99848616 | 0.00419523 | 0.00477009 | 0.99910087 | 2.05e-15 | 3.54e-15 |
| lih_casscf | 0.99801843 | 0.00781091 | 0.00854537 | 0.99879477 | 2.13e-15 | 3.54e-15 |
| lih_cisd | 0.99833398 | 0.00661548 | 0.00664398 | 0.99816070 | 1.90e-15 | 3.41e-15 |
| lih_tda | 0.99846883 | 0.00655467 | 0.00678468 | 0.99918969 | 2.71e-15 | 3.32e-15 |
| lih_tddft | 0.99846883 | 0.00618606 | 0.00652725 | 0.99918359 | 2.81e-15 | 3.32e-15 |
| lih_plus_casscf_hf | 0.99834954 | 0.00231556 | 0.00173750 | 0.99994342 | 2.35e-15 | 3.32e-15 |
| lih_plus_casscf | 0.99834917 | 0.00629382 | 0.00458160 | 0.99938260 | 2.73e-15 | 3.32e-15 |
| lih_plus_cisd | 0.99834917 | 0.00626820 | 0.00456098 | 0.99938456 | 1.71e-15 | 3.61e-15 |
| lih_plus_tda | 0.99842802 | 0.00293289 | 0.00203494 | 0.99984711 | 1.81e-15 | 3.87e-15 |
| lih_plus_tddft | 0.99842802 | 0.00295125 | 0.00204742 | 0.99984472 | 2.18e-15 | 3.87e-15 |

The additional sequential LiH/CASSCF run reproduced a current MO metric error of 0.04823816 and a self-overlap identity error of 0.01894549. It completed and the direct and request/result overlap evaluations agreed; its nonorthonormal orbitals are a known limitation of the existing warm-start energy path, not a passing physical-normalization check. The default independent calculations isolate the overlap representation from this issue.

CASSCF overlaps in both modes omit inactive-core contributions. Restricted TDDFT/TDA overlaps now use normalized singlet determinant expansions; unrestricted overlaps use their alpha/beta determinant expansion. Full TDDFT still omits Y amplitudes from the overlap pseudo-wavefunction. The example README explains these limits and reproduction commands. No molecular dynamics or production ensembles were run.
