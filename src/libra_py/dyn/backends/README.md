This is crucial if we care about performance.

numpy_backend.py
•	NumPy wrappers

torch_backend.py
•	PyTorch tensors
•	autograd compatibility

Everything above talks to abstract tensor operations.


`TorchBackend(device="cpu", dtype="float64")` implements allocation and the
contractions needed by the integrated coherent Gaussian path. Complex storage
uses the matching precision; mixed real/complex einsum operands are promoted.
Constructing the facade imports Torch lazily. `NumpyBackend(dtype="float32")`
also supports explicit storage precision (default float64). Shared helpers in
`arrays.py` preserve array dtype/device and detach only at explicit I/O boundaries.
This does not imply that every older dynamics routine is Torch-compatible.
