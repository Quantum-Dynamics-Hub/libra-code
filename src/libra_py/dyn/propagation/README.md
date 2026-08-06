How things move

Split by physical role, not by algorithm.

nuclear.py
• 	classical propagation
• 	velocity Verlet, Langevin, etc.
electronic.py
• 	TDSE solvers
• 	orthonormalization
• 	NAC handling
• 	`Dynamics.cpp` electronic-integrator dispatch
coupled.py
• 	coupled TBF propagation
•	Ehrenfest-like schemes
• 	multi-TBF TDSEs
integrators.py
• 	generic time integrators
• 	step splitting
• 	adaptive time stepping

✔ easily swappable
✔ testable in isolation

Electronic integrator selectors
-------------------------------

`propagate_electronic_method` implements the amplitude-based branches selected
by `DynControlParams.electronic_integrator` (the C++ local variable `method`).

For adiabatic amplitudes:

- `-1`: no electronic propagation;
- `0`, `1`, `2`: crude-split, symmetric-split, and averaged-Hamiltonian LD;
- `3`: one-point old vibronic Hamiltonian;
- `4`: two-point old/current vibronic Hamiltonian;
- `5`: two-point vibronic Hamiltonian with `T† H_new T` correction;
- `6`: the implemented C++ old-Hvib half-step without reordering;
- `7`: reordered old/new Hvib half-steps;
- `8`: intended old/new half-step form of the experimental new-LD branch;
- `10`--`15`: rotation-based equivalents of `0`--`5`;
- `100`--`115`: corresponding C++ aliases.

For diabatic amplitudes, methods `0`--`3` and `100`--`103` provide midpoint,
split, and overlap-aware non-Hermitian propagation. Density-matrix selectors
(`rep_tdse` 2 and 3 in C++) are not amplitude integrators and remain outside
this function.

The projector convention is the one used by `Dynamics.cpp`: `T @ C` maps old
dynamically consistent coefficients into the new raw basis, while
`T† @ H_new @ T` expresses a new raw Hamiltonian in the old tracked labels.

Documented differences from C++
-------------------------------

- C++ method 8 reads `Hvib` before assigning it in that branch. Python uses
  the documented old/new half-step interpretation.
- The C++ condition for rotation method 12 repeats `method == 12`, making
  alias 112 unreachable. Python accepts 112.
- C++ method 6 builds an old/new sum but then ignores it. Python preserves the
  actual old-Hvib half-step behavior and documents it explicitly.
- The greedy C++ state tracker reverses cycles longer than two because of its
  permutation-composition order. Python retains `perm[old] = new`, which is
  required by projector construction and the propagation equations above.
