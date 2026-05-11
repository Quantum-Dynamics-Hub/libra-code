How things move

Split by physical role, not by algorithm.

nuclear.py
• 	classical propagation
• 	velocity Verlet, Langevin, etc.
electronic.py
• 	TDSE solvers
• 	orthonormalization
• 	NAC handling
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

