The untouchable foundation

These modules should almost never change.


storage.py
•	TensorStorage
•	backend-agnostic tensor allocation
•	no physics
views.py
•	TensorView
•	pure indexing logic
•	no algorithms
tbf.py
•	TrajectoryBasisFunction
•	identity, metadata, lifecycle flags
•	no propagation
trajectory.py
•	Trajectory
•	bookkeeping of TBFs
•	wavefunction-level properties
manager.py
•	TrajectoryManager
•	creation/deletion logic
•	global registry

Forbidden here:
•	force evaluation
•	time propagation
•	spawning logic


