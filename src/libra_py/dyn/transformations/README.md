🧠 5. What each submodule becomes
🟦 local_diabatization.py

Your current logic:

T matrices from overlap
basis tracking per TBF/trajectory
“LD propagation corrections”

✔ still exists, but now scoped properly


🟨 basis_rotation.py

General linear transforms:

C' = U† C
H' = U† H U

Used for:

adiabatic ↔ diabatic
symmetry transforms
representation changes

This is also where TensorStorage wrappers for amplitude representation changes
belong. Hamiltonian construction can depend on amplitudes only for genuinely
amplitude-dependent models; routine amplitude basis changes stay here or in
propagation orchestration.



🟥 orthogonalization.py

Numerical stability tools:

Löwdin orthogonalization
Gram-Schmidt
SVD-based cleanup


🟪 gauge.py

Phase and continuity fixes:

phase matching
NAC smoothing
sign continuity of eigenvectors


🟩 state_tracking.py

Python translations of `dyn_projectors.cpp` state-identity algorithms:

- diagonal-overlap phase corrections;
- greedy, Hungarian/Munkres-Kuhn, and stochastic reorderings;
- energy-aware overlap cost matrices;
- force-based crossing cost matrices;
- permutation matrices and active-state remapping.

The functions retain the C++ permutation convention: `perm[i]` is the new
label of old state `i`.


🟧 projectors.py

Instantaneous and cumulative projection updates from `dyn_ham.cpp`, including
the original `state_tracking_algo` option numbers. Option `-1` imports
`orthogonalized_T` from `local_diabatization.py`; options `5` and `6` implement
the SVD-based LD variants. The module also provides storage-backed
`update_proj_adi` for the dynamics engine.
