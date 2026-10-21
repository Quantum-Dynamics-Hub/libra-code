# Gaussian wavepacket integrals

```python
import numpy as np
from libra_py.dyn.gwp import integrals as gwp

q = np.array([0.0, 0.5])
p = np.array([0.2, -0.1])
gamma = np.zeros(2)
alp = np.ones(2)
S = gwp.gwp_overlap(q, p, gamma, alp, q, p, gamma, alp)  # 1
T = gwp.gwp_kinetic(q, p, gamma, alp, q, p, gamma, alp)
```

Pass PyTorch tensors instead of NumPy arrays to select PyTorch. Tensor device,
floating dtype and autograd graphs are preserved; NumPy users do not need
PyTorch or the C++ extension installed. Mixed tensor dtypes are promoted. All
parameters are real, widths must be finite and positive, and tensors must share
a device. Single-packet inputs are scalars, `(ndof,)` or `(ndof, 1)` arrays;
coordinate and momentum lengths must match. Widths and phases can each have
a single broadcast entry. In the eight-argument pair interfaces, use vector
phases with vector widths: a scalar third argument and an array fourth
argument select the legacy overload. For example, use `np.zeros(ndof)` for
zero phases with per-coordinate widths. Vector-valued results have the shape of the first
coordinate input. Functions accept positional arguments.

## API and conventions

For the usual eight-argument calls, arguments are
`q1, p1, gamma1, alp1, q2, p2, gamma2, alp2`, with atomic units (hbar = 1).
Each packet is

```text
G(r) = product_i (2*alp_i/pi)^(1/4)
       * exp(sum_i [-alp_i*(r_i-q_i)^2 + i*p_i*(r_i-q_i) + i*gamma_i]).
```

The phase is per coordinate, so a scalar phase in this interface is repeated
`ndof` times. To supply a total phase, distribute it among the coordinates.

| Function | Result / additional arguments |
| --- | --- |
| `gwp_overlap` | Integral `conj(G1)*G2` |
| `gwp_dipole` | Vector of integrals `conj(G1)*r_i*G2` |
| `gwp_coupling` | Vector of integrals `conj(G1)*dG2/dr_i` |
| `gwp_kinetic` | Integral of `conj(G1)*(-1/2 sum_i invM_i*d_i^2)*G2`; optional ninth inverse-mass vector, default 1 |
| `gwp_product_decomposition` | Tuple `(scale, q, p, gamma, alp)` such that `G1*G2 = scale*G` |
| `gwp_overlap` with 12 parameters | Three packets: integral `conj(G1*G2)*G3`, matching the C++ conjugation convention |
| `gwp_overlap_matrix` | Eight `(ndof,ntraj)` arrays, result `(ntraj1,ntraj2)`; state overload is `q1,p1,gamma1,alp1,state1,q2,p2,gamma2,alp2,state2` |
| `gwp_kinetic_matrix` | Eight `(ndof,ntraj)` arrays followed by inverse masses; result `(ntraj1,ntraj2)` |
| `gwp_value(r,q,p,alp)` | Wavefunction value, zero phase, hbar = 1 |
| `gwp_value(r,q,alp)` | Normalized probability density `abs(G)**2`, not the zero-momentum wavefunction |
| `gwp_deriv` | Spatial gradient of the corresponding `gwp_value` overload |

Overlap, dipole, coupling and kinetic also accept the C++ equal-width overload
`R1,P1,gamma1,R2,P2,gamma2,alp,hbar`. Here coordinates and momenta must be
vectors, and phases, width and positive hbar are scalars. These phases are
**total** phases. `gwp_value(r,R,P,gamma,alp,hbar)` uses the same convention;
its gradient is also available in Python.

As documented by the C++ interface, the legacy kinetic operator is
`-1/2 Laplacian`, even when hbar differs from 1. Multiply by `hbar**2/mass`
for the physical kinetic energy in other units. The general inverse-mass
interface supports diagonal masses, not a full mass tensor.

The matrix routines broadcast over all trajectory pairs. Memory use scales
as `ndof*ntraj1*ntraj2`; split large basis sets into blocks when needed.

## Formula verification

In one dimension let `A = a1+a2`, `dq = q2-q1`, `dp = p2-p1`,
`b = 2*a1*a2/A`, and `P = (a2*p1+a1*p2)/A`. Completing the square in
`conj(G1)*G2` gives a complex Gaussian mean and variance

```text
mu = (a1*q1+a2*q2)/A + i*dp/(2*A)
variance = 1/(2*A)
S = (4*a1*a2/A^2)^(1/4)
    * exp(-b*dq^2/2 - dp^2/(4*A) + i*(gamma2-gamma1-P*dq))
D = b*dq + i*P
<r>_12 = mu*S
<d/dr>_12 = D*S
<-1/2 d^2/dr^2>_12 = (b-D^2)*S/2.
```

For separable multidimensional packets, multiply the overlaps and sum the
kinetic prefactors (weighted by inverse masses). These expressions require
no division by the overlap and remain well-defined when it underflows.

For the pointwise product of two normalized packets,

```text
a = a1+a2
q = (a1*q1+a2*q2)/a
p = p1+p2
gamma = gamma1+gamma2 + (q2-q1)*(a2*p1-a1*p2)/a
scale = (2*a1*a2/(pi*a))^(1/4) * exp(-a1*a2/a*(q2-q1)^2).
```

The scales multiply over dimensions. The triple integral uses this product
and then the two-packet overlap; both first packets are conjugated. This is
not an integral with just G1 conjugated, despite the ambiguous older C++
notation `<G1|G2|G3>`.

The audit corrected these C++ issues:

* Two `gwp_value` overloads accumulated into an uninitialized real exponent.
  Their derivative helpers inherited that undefined behavior.
* Both product-decomposition overloads omitted `pi^(-ndof/4)`, also affecting
  the triple overlap.
* The equal-width kinetic overload used `-alp` instead of `-ndof*alp` and
  omitted `/hbar` in its imaginary prefactor.
* The scalar product-decomposition Python binding required mutable `double&`
  outputs, which Python floats cannot supply. An eight-argument overload now
  returns `(scale,q,p,gamma,alp)`; existing C++ signatures remain available.

## Tests

Rebuild the extension before comparing against it:

```bash
cmake --build _build --target gwp libra_core -j 4
PYTHONPATH=src:_build/src python -m pytest -q unittests/test_libra_py/test_dyn/test_gwp/test_integrals.py
```

The tests compare both backends against Python-exposed C++ overloads, check
independent numerical quadrature of wavefunctions and their derivatives,
verify Hermitian/anti-Hermitian symmetries and diagonal limits, and use
PyTorch Jacobians compared against centered finite differences. C++/PyTorch tests skip when the corresponding optional
module is absent. The gradient checks use scalar autograd calls and do not
require SymPy. CUDA device/dtype checks run when CUDA is
available; CPU checks always run when PyTorch is installed.

## Examples

Runnable NumPy and PyTorch examples are in
[`examples/libra_py/dyn/gwp/integrals`](../../../../examples/libra_py/dyn/gwp/integrals/README.md).
The implementation lives in `integrals.py`; the package `__init__.py` re-exports
its public functions for convenience.

For test coverage, numerical tolerances, and dependency requirements, see the
[test documentation](../../../../unittests/test_libra_py/test_dyn/test_gwp/README.md).
Public functions also document their individual overloads in `help(function)`.
