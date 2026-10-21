# Gaussian wavepacket integral examples

Run these scripts from the repository root with `src` on your Python path:

```bash
PYTHONPATH=src python examples/libra_py/dyn/gwp/integrals/example_numpy.py
PYTHONPATH=src python examples/libra_py/dyn/gwp/integrals/example_torch.py
```

Both scripts perform the same calculations on the same input values, with
identical print labels and ordering. `example_numpy.py` uses NumPy arrays;
`example_torch.py` uses PyTorch tensors. Each demonstrates:

1. Overlap, kinetic energy, position, and spatial derivative matrix elements.
2. Wavefunction value, spatial gradient, and probability density at one point.
3. Product decomposition (all five outputs) and the triple integral.
4. Squared overlap and its derivatives with respect to packet 1's centers and widths.
5. Basis overlap and kinetic matrices, and overlaps with electronic-state labels.

NumPy computes the parameter derivatives using centered finite differences
with step `1e-6`; PyTorch computes them with autograd. They agree to numerical
precision appropriate for finite differences. Printed values retain their
native array/tensor representation, so PyTorch also shows dtype and autograd
metadata. Pass `--device cuda` to the PyTorch script to use a supported GPU.

The NumPy example requires NumPy; the PyTorch example also requires PyTorch
for scalar autograd. Neither script requires SymPy or the C++ extension.
Both import `libra_py.dyn.gwp.integrals`.

Packet phases are per coordinate and widths must be positive. The examples use
atomic units, with hbar = 1 and diagonal inverse masses. Matrix inputs have shape
`(ndof, ntraj)`, with one packet per column.

See the [API and formula documentation](../../../../../src/libra_py/dyn/gwp/README.md)
for overload conventions and details of the C++ formula corrections.

## Reading the inputs

Start with `example_numpy.py`. It names each parameter separately and passes
all eight packet parameters explicitly, in the order
`q1, p1, gamma1, alp1, q2, p2, gamma2, alp2`.

| Parameter | Meaning | Shape in these examples |
| --- | --- | --- |
| `q1`, `q2` | Packet centers in bohr | `(2,)`: one value per spatial DOF |
| `p1`, `p2` | Momenta in atomic units | `(2,)` |
| `gamma1`, `gamma2` | Phases per DOF; their sum is the total phase at hbar = 1 | `(2,)` |
| `alp1`, `alp2` | Positive widths in bohr^-2; density variance is `1/(4*alp)` | `(2,)` |
| `inv_mass` | Inverse mass per DOF; `[1, 0.5]` corresponds to masses `[1, 2]` in electron masses | `(2,)` |
| `r` | A single evaluation point in two-dimensional coordinate space | `(2,)`, not a batch of points |
| `q_basis`, etc. | Parameters for two packets, with rows = DOFs and columns = packets | `(2, 2)` |
| `states` | Electronic-state label for each packet | `(2,)`, indexed by packet rather than DOF |

Both examples stack each named parameter separately to put packets in
columns: NumPy uses `axis=1`, and PyTorch uses `dim=1`. PyTorch stacking
also preserves the autograd graph. Each script keeps the input construction
and integral calls explicit rather than hiding them behind a backend wrapper.

## Understanding the outputs

The overlap is approximately `0.705963 + 0.116059j`. Kinetic energy with the
chosen inverse masses is approximately `0.357611 + 0.078937j`. These are
transition matrix elements between different packets, so they can be complex.
The position and derivative coupling each have two components. A self-overlap
is 1; using the same basis on both sides gives Hermitian overlap and kinetic
matrices. Assigning electronic labels `[0, 1]` removes the off-diagonal overlaps
because those electronic states are orthogonal.

`gwp_value(r, q, p, alp)` evaluates a zero-phase wavefunction. The example uses
packet 1, whose phase is zero. `gwp_value(r, q, alp)` instead returns its
probability density, and the example verifies that it equals the squared
magnitude of the wavefunction. Product decomposition returns five outputs; the
new packet must be multiplied by `scale` to reproduce the pointwise product.
The triple overlap conjugates **both** first packets.

## Parameter gradients

In the PyTorch example, only `q1` and `alp1` set `requires_grad=True`. Calling `backward()` on the real
objective `abs(overlap)**2` fills their `.grad` arrays, giving sensitivities of
the squared overlap to packet centers and widths. This does not update the
parameters or optimize the packets. These parameter gradients differ from the
spatial gradient with respect to `r` computed by `gwp_deriv`.

The NumPy example perturbs each component of `q1` and `alp1` separately,
evaluating the squared overlap at plus/minus `1e-6` while keeping all other
parameters fixed. The centered differences estimate the same derivatives
as autograd without changing the original input arrays.

The squared overlap is approximately `0.511853`; its center gradient is about
`[0.2529, -0.6879]` and its width gradient about `[-0.0528, -0.1006]`.
All tensors use float64 on the selected device; complex integral outputs use
complex128. To run on CUDA explicitly:

```bash
PYTHONPATH=src python examples/libra_py/dyn/gwp/integrals/example_torch.py --device cuda
```

Both scripts include small consistency assertions and print their results.
They do not write output files. More extensive checks are documented in the
[test README](../../../../../unittests/test_libra_py/test_dyn/test_gwp/README.md).
