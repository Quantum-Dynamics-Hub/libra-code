# Hamiltonians

Hamiltonian code is responsible for constructing and transforming tensor
quantities stored in `dyn.core.TensorStorage`.

`TensorStorage` is the single owner of Hamiltonian data:

- `ham_adi`, `ham_dia`
- `hvib_adi`, `hvib_dia`
- `nac_adi`, `nac_dia`
- `dc1_*`, `d1ham_*`, `d2ham_*`
- `basis_transform`, overlaps, ordering, and phase corrections

The package intentionally does not keep separate `adiabatic`, `diabatic`, or
`state` modules while those concepts are only storage fields and tensor
operations. Representation-specific kernels can be added later as focused
modules once they contain real behavior.

`HamiltonianEngine` evaluates a model function, writes returned tensors into
the active trajectory/TBF storage slices, applies representation transforms,
and builds vibronic Hamiltonians as:

```text
Hvib = H - i * NAC
```

Hamiltonian modules should not generally own electronic-amplitude
transformations. Those live in `dyn.transformations.basis_rotation` and are
called by propagation or higher-level workflows. A Hamiltonian model may still
read amplitudes from `TensorStorage` when the Hamiltonian is genuinely
amplitude-dependent, such as in exact-factorization or other nonlinear
mean-field constructions.

## nHamiltonian Translation Map

The old C++ `nHamiltonian` class used parent/child Hamiltonian nodes. The
Python dyn prototype keeps one flat `TensorStorage` object instead, with
trajectory/TBF slices selected by `traj.id` and `traj.tbf_ids`.

Translated method groups:

- `_aux`: `aux.py` provides shape checks, dict-or-attribute result adapters,
  active-slice access, and active-slice writes.
- `_basic`: `basic.py` provides allocation, reset, content-copy, and memory
  status helpers for Hamiltonian storage.
- `_compute_diabatic`: `diabatic.py::compute_diabatic` calls a model function,
  validates optional returned tensors, writes them into storage, and builds
  `hvib_dia` when needed.
- `_compute_adiabatic`: `adiabatic.py::compute_adiabatic` handles direct model
  returns, and `compute_adiabatic_from_diabatic` uses
  `dyn.utils.linalg.generalized_eigh` to solve `H_dia U = S_dia U E` and
  transform first-derivative data.
- `_basis_transform`: moved to `dyn.transformations.basis_rotation`, because
  amplitude and operator rotations are representation transforms, not
  Hamiltonian construction.
- `_compute_Ehrenfest` and `_compute_Ehrenfest_forces`: `ehrenfest.py`
  provides energies, force tensors, and force vectors in both representations.

Not translated:

- parent/child traversal, `lvl`, `split`, and `full_id` routing. Use
  `TensorStorage` trajectory/TBF slices instead.
- adiabatic-state reordering updates. The C++ `update_ordering` path is tightly
  coupled to legacy state tracking and has no current Python call site, so it
  should be reintroduced only with a tested ordering/gauge workflow.
- pointer ownership modes. TensorStorage owns arrays directly.
- C++ destructor and external-reference setters. They do not have Python
  equivalents in this storage model.
