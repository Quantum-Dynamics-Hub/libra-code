# PySCF electronic-state overlaps

This example computes overlaps between LiH geometries at 3.0 and 3.1 Bohr with STO-3G. It also runs LiH⁺ to exercise the open-shell implementations. Each calculation retains two electronic states. These are inexpensive electronic-structure calculations, with no trajectories or gradients.

Run from the Libra repository root in an environment containing PySCF, NumPy, and Libra:

```bash
python examples/libra_py/packages/pyscf/02_example_overlaps/overlaps.py --source-tree --output /tmp/pyscf_overlap_comparison
```

You can also run directly from the example folder:

```bash
cd examples/libra_py/packages/pyscf/02_example_overlaps
python overlaps.py --source-tree --output /tmp/pyscf_overlap_local
```

The source-tree location is resolved relative to the script, so both invocation locations work. Relative output paths are resolved against your current working directory.

`--source-tree` uses this checkout's Python implementation alongside the installed Libra package. Omit it to use the installed implementation. Output directories must be new. All calculations use one PySCF thread.

| `--method` | Overlap representation |
|---|---|
| `casscf_hf` | CASSCF strategy, historical auxiliary CASCI roots in HF orbitals |
| `casscf` | Optimized CASSCF active orbitals and their stored CI roots |
| `cisd` | Restricted CISD for LiH; UCISD for LiH⁺; no frozen orbitals in this example |
| `tda` | RKS PBE0 TDA: normalized spin-adapted singlet determinant expansion; UKS: alpha/beta determinant expansion |
| `tddft` | Same determinant construction using normalized X amplitudes from full TDDFT; Y amplitudes are not represented |

`--method all` and `--system all` are the defaults: ten comparisons in total. Other system choices are `lih` and `lih_plus`. CASSCF uses CAS(2 electrons, 2 orbitals) for LiH and CAS(1 electron, 2 orbitals) for LiH⁺, both with one doubly occupied inactive orbital. The open-shell calculation targets a doublet. TDDFT/TDA retain the reference and one response root; CASSCF/CISD retain two CI roots. These method-dependent roots need not describe identical physical states.

## Selecting the CASSCF overlap basis

```python
from libra_py.packages.pyscf.implementations import CASSCF

# Default remains auxiliary CASCI on HF orbitals.
hf_overlaps = CASSCF(norbcas=2, nelecas=2, nroots=2)
# Equivalent explicit selection:
hf_overlaps = CASSCF(norbcas=2, nelecas=2, nroots=2,
                    overlap_orbitals="hf")

# Use the optimized CASSCF orbitals AND the matching stored CI coefficients.
optimized_overlaps = CASSCF(norbcas=2, nelecas=2, nroots=2,
                           overlap_orbitals="casscf")
```

The new setting affects only the time-overlap representation. CASSCF energies, gradients, state averaging, spin handling, and the existing column-sign convention are unchanged. The setting survives `copy()`, including the copies created by the dynamics adapter. Invalid values raise `ValueError`. The optimized option requires CASSCF orbitals and CI roots in both snapshots; it never silently falls back to HF orbitals.

For independently calculated snapshots, the call is:

```python
# Energies have already been calculated for both geometries.
U = optimized_overlaps.compute_time_overlap(current_state, previous_state)
```

Despite the argument order, `U[i, j] = <previous root i | current root j>`. In the ordinary sequential interface, request `ES_Request(n_singlets=2, time_overlap=True)` and read `result.time_overlap` after the second geometry. The direct strategy returns `None` on the first geometry because there is no previous snapshot; the dynamics adapter supplies its initial identity matrix.

## Restricted TDDFT and TDA determinant overlaps

For an RKS reference, each normalized X-only response root is represented as

```text
|Psi_I> = sum_ia X^I_ia (|Phi_i^a alpha> + |Phi_i^a beta>)/sqrt(2),
sum_ia |X^I_ia|^2 = 1.
```

The implementation evaluates the exact alpha/beta determinant overlap of the selected pseudo-states with the full cross-geometry MO overlap matrix. For a nonsingular occupied block it uses the algebraically equivalent determinant-minor formulas, avoiding an explicit quadratic sum over configurations; an explicit determinant expansion is the fallback for a singular or ill-conditioned occupied block. Consequently, the closed-shell reference overlap is `det(S_occ)^2`; reference–excited elements contain the singlet `sqrt(2)` factor and determinant cofactors; and excited–excited elements retain occupied–virtual cross terms. This replaces the former occupied/virtual contraction.

The public calculation is unchanged:

```python
from libra_py.packages.pyscf.implementations import TDDFT
from libra_py.packages.pyscf.interfaces import ES_Request

method = TDDFT(atom_labels=("Li", "H"), nexc=1, basis="sto-3g",
               xc="pbe0", use_tda=True)
request = ES_Request(n_singlets=2, time_overlap=True)

first = method.compute_result(geometry_t, request)
second = method.compute_result(geometry_t_plus_dt, request)
U = second.time_overlap       # U[I,J] = <Psi_I(t)|Psi_J(t+dt)>
```

Run only the closed-shell TDA or TDDFT examples with:

```bash
python examples/libra_py/packages/pyscf/02_example_overlaps/overlaps.py \
  --source-tree --system lih --method tda --output /tmp/lih_tda_overlaps
python examples/libra_py/packages/pyscf/02_example_overlaps/overlaps.py \
  --source-tree --system lih --method tddft --output /tmp/lih_tddft_overlaps
```

## Outputs and interpretation

Each comparison writes a JSON record and an overlap CSV; `summary.json` collects all records. The terminal prints the matrices. JSON includes both sets of energies, the current-state self-overlap, singular values, the overlap orthogonality defect, and the full-MO metric error `||C† S_AO C − I||F` (maximum over spin channels for unrestricted references). All matrices use the interface's column-phase alignment. Off-diagonal signs can vary between independent electronic calculations; comparisons of absolute matrix elements remove that sign ambiguity but do not solve root tracking.

Both CASSCF options retain the existing **active-space-only contraction**. Neither includes inactive-core determinant factors or core–active cross terms. The optimized option therefore improves orbital/CI consistency without claiming a complete all-electron CASSCF overlap. CISD uses the existing PySCF determinant overlap routines. Restricted TDDFT/TDA now uses exact determinant-minor expressions, with a determinant-sum fallback, for the selected normalized singlet X-only pseudo-states; the unrestricted path uses determinant sums of its selected X-only pseudo-states. Full TDDFT's Y amplitudes remain omitted, so these are pseudo-wavefunction overlaps rather than overlaps of the complete linear-response object.

An overlap between different geometries need not be unitary in a truncated state space. Even identical-geometry X-only TDDFT pseudo-states need not be mutually orthogonal when several response roots are retained. Report diagnostics rather than automatically replacing the matrices by identities.

## Sequential CASSCF warm-start diagnostic

The default comparison solves each geometry independently. To exercise the current sequential strategy instead:

```bash
python examples/libra_py/packages/pyscf/02_example_overlaps/overlaps.py --source-tree --system lih --method casscf --sequential --output /tmp/pyscf_overlap_sequential
```

The existing CASSCF energy routine passes the previous geometry's MO coefficients to the next optimization without reorthonormalizing them in the new AO metric. In the tested 0.1-Bohr displacement, the resulting MO metric error is about 0.0482 and the current-state self-overlap differs from identity. This is a separate pre-existing warm-start issue, not corrected or concealed by the new overlap choice. Independent solves have orthonormal orbitals to numerical precision. Before using the optimized option for quantitative sequential dynamics, the warm-start orbital transfer and inactive-core overlap treatment require attention. `--sequential` exposes the issue and checks agreement between direct overlap evaluation and the public request/result route.

## Unit tests

The CASSCF tests in `unittests/test_libra_py/test_packages/test_pyscf/test_casscf_overlaps.py` compare both representations against an independent small determinant expansion. The restricted-response tests in `test_tddft_overlaps.py` cover the analytic two-electron orbital rotation, singlet normalization, the squared closed-shell determinant, orthogonality at identical geometry, invalid amplitude dimensions, the determinant-minor formula against an explicit multielectron determinant sum, and the public displaced-geometry TDA request.

For an installed build containing this revision:

```bash
python -m pytest -q unittests/test_libra_py/test_packages/test_pyscf
```

To test this checkout while retaining installed compiled modules, run from the repository root:

```bash
python - <<'PY'
from pathlib import Path
import libra_py
import pytest
libra_py.__path__.insert(0, str(Path.cwd() / "src/libra_py"))
raise SystemExit(pytest.main(["-q", "unittests/test_libra_py/test_packages/test_pyscf"]))
PY
```

See [VALIDATION.md](VALIDATION.md) for the completed small calculations and their numerical diagnostics.
