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
