"""Independent quadrature, algebraic properties, autograd and C++ parity.

Run: PYTHONPATH=src:_build/src python -m pytest -q unittests/test_libra_py/test_dyn/test_gwp/test_integrals.py
C++ and Torch comparisons skip only when their optional dependency is absent.
"""
import sys

import numpy as np
import pytest

from libra_py.dyn.gwp import integrals as gwp


@pytest.fixture(params=['numpy', 'torch'])
def backend(request):
    """Run each numerical test with NumPy and, when installed, float64 PyTorch."""
    if request.param == 'torch':
        torch = pytest.importorskip('torch')
        return lambda x: torch.as_tensor(x, dtype=torch.float64)
    return np.asarray


@pytest.fixture
def cpp():
    """Load the rebuilt Python C++ extension; skip comparisons if it is unavailable."""
    return pytest.importorskip('cyglibra_core' if sys.platform == 'cygwin' else 'liblibra_core')


def array(x):
    """Convert comparison results to NumPy; detaching is confined to test assertions."""
    return x.detach().cpu().numpy() if hasattr(x, 'detach') else np.asarray(x)


def close(x, y, **kwargs):
    """Compare float64/complex128 results with rtol=2e-12 and atol=2e-13."""
    np.testing.assert_allclose(array(x), array(y), rtol=2e-12, atol=2e-13, **kwargs)


def packets(n=3, columns=None, seed=421):
    """Generate reproducible (q,p,gamma,alp) test inputs.

    n is the DOF count; columns=None gives (n,) vectors, otherwise (n,columns)
    basis matrices. Widths lie in [0.3,2.0) to avoid singular or extremely
    narrow Gaussians. Fixed seeds keep numerical regressions reproducible.
    """
    rng = np.random.default_rng(seed)
    shape = (n,) if columns is None else (n, columns)
    return [rng.normal(size=shape), rng.normal(size=shape),
            rng.normal(size=shape), rng.uniform(0.3, 2.0, size=shape)]


def matrix(cpp, x):
    """Copy NumPy values to C++ MATRIX, converting 1D inputs into column vectors."""
    x = np.asarray(x)
    if x.ndim == 1:
        x = x[:, None]
    result = cpp.MATRIX(*x.shape)
    for index in np.ndindex(x.shape):
        result.set(*index, float(x[index]))
    return result


def cpp_array(x):
    """Extract C++ matrix entries for comparisons, leaving scalar results unchanged."""
    if hasattr(x, 'num_of_rows'):
        return np.array([[x.get(i, j) for j in range(x.num_of_cols)]
                         for i in range(x.num_of_rows)])
    return x


@pytest.mark.parametrize('name', ['overlap', 'dipole', 'coupling', 'kinetic'])
@pytest.mark.parametrize('kind', ['scalar', 'vector', 'column', 'legacy'])
def test_cpp_integrals(cpp, backend, name, kind):
    """Compare four operators across both backends and all pair input forms.

    Covers scalar, vector, column, and equal-width legacy calls. The legacy
    case uses three DOFs and hbar=1.7 to exercise the corrected dimension
    and hbar factors. C++ vector results are flattened only for comparison.
    """
    a = packets(seed=11)+packets(seed=12)
    if kind == 'scalar':
        a = [float(x[0]) for x in a]
        ca = a
    elif kind == 'legacy':
        a = [a[0], a[1], 0.21, a[4], a[5], -0.7, 1.2, 1.7]
        ca = [matrix(cpp, x) if np.ndim(x) else x for x in a]
    else:
        if kind == 'column':
            a = [x[:, None] for x in a]
        ca = [matrix(cpp, x) for x in a]
    actual = getattr(gwp, 'gwp_'+name)(*[backend(x) for x in a])
    expected = cpp_array(getattr(cpp, 'gwp_'+name)(*ca))
    close(array(actual).reshape(-1), np.asarray(expected).reshape(-1))


def test_cpp_mass_and_matrices(cpp, backend):
    """Compare rectangular basis matrices, electronic masks and inverse masses.

    Sets contain two and four packets in three DOFs, catching swapped axes
    or assumptions of square outputs. Distinct inverse masses exercise
    coordinate weighting, both for matrix and single-pair kinetic integrals.
    """
    a = packets(columns=2)+packets(columns=4, seed=71)
    masses = np.array([0.2, 1.7, 0.8])
    ca = [matrix(cpp, x) for x in a]
    pa = [backend(x) for x in a]
    close(gwp.gwp_overlap_matrix(*pa), cpp_array(cpp.gwp_overlap_matrix(*ca)))
    close(gwp.gwp_kinetic_matrix(*pa, backend(masses)),
          cpp_array(cpp.gwp_kinetic_matrix(*ca, matrix(cpp, masses))))
    states1, states2 = [0, 1], [1, 0, 0, 1]
    s1, s2 = cpp.intList(), cpp.intList()
    for s in states1:
        s1.append(s)
    for s in states2:
        s2.append(s)
    close(gwp.gwp_overlap_matrix(*pa[:4], states1, *pa[4:], states2),
          cpp_array(cpp.gwp_overlap_matrix(*ca[:4], s1, *ca[4:], s2)))
    av = [x[:, 0] for x in a]
    close(gwp.gwp_kinetic(*[backend(x) for x in av], backend(masses)),
          cpp.gwp_kinetic(*[matrix(cpp, x) for x in av], matrix(cpp, masses)))


@pytest.mark.parametrize('kind', ['density', 'wavefunction', 'legacy'])
def test_cpp_values(cpp, backend, kind):
    """Compare density, wavefunction and legacy values against C++.

    The three/four-argument gradients are compared too. These calls exercise
    the C++ exponent accumulators fixed by the audit; the legacy value uses
    a nonzero total phase and hbar=1.6.
    """
    q, p, _, alp = packets()
    r = np.array([0.5, -0.7, 0.1])
    args = [r, q, alp] if kind == 'density' else [r, q, p, alp]
    if kind == 'legacy':
        args = [r, q, p, 0.23, 0.7, 1.6]
    ca = [matrix(cpp, x) if np.ndim(x) else x for x in args]
    pa = [backend(x) for x in args]
    close(gwp.gwp_value(*pa), cpp.gwp_value(*ca))
    if kind != 'legacy':
        close(gwp.gwp_deriv(*pa), cpp_array(cpp.gwp_deriv(*ca))[:, 0])


def test_cpp_product_and_triple(cpp, backend):
    """Compare all decomposition outputs and the three-packet integral.

    The vector C++ call mutates MATRIX outputs; its scalar Python wrapper
    returns a tuple. Both are checked, including the corrected pi factor.
    The triple overlap uses conjugation of BOTH first packets.
    """
    args = packets()+packets(seed=512)
    outputs = [cpp.MATRIX(3, 1) for _ in range(4)]
    scale = cpp.gwp_product_decomposition(*[matrix(cpp, x) for x in args], *outputs)
    actual = gwp.gwp_product_decomposition(*[backend(x) for x in args])
    close(actual[0], scale)
    for x, y in zip(actual[1:], outputs):
        close(x, cpp_array(y)[:, 0])
    args += packets(seed=99)
    close(gwp.gwp_overlap(*[backend(x) for x in args]),
          cpp.gwp_overlap(*[matrix(cpp, x) for x in args]))
    # The scalar Python binding returns the output parameters as a tuple.
    scalars = [float(x[0]) for x in args[:8]]
    for x, y in zip(gwp.gwp_product_decomposition(*map(backend, scalars)),
                    cpp.gwp_product_decomposition(*scalars)):
        close(x, y)


def wave(x, q, p, gamma, alp):
    """Evaluate the defining 1D wavefunction independently of the integral routines."""
    return (2*alp/np.pi)**0.25*np.exp(-alp*(x-q)**2+1j*(p*(x-q)+gamma))


@pytest.mark.parametrize('seed', [4, 18, 31])
def test_quadrature(backend, seed):
    """Verify analytic formulas against direct 1D numerical integration.

    For each fixed seed, 400-point Gauss-Legendre quadrature integrates the
    defining wavefunctions and explicit spatial derivatives on a wide finite
    interval. This reference does not call C++ or the completed-square
    formulas. Also check the product identity pointwise and integrate
    conj(G1*G2)*G3 to detect normalization or conjugation errors.
    """
    # Gauss-Legendre integration on an interval extending 12 Gaussian widths
    # beyond every center; evaluate the integrands directly, independently of
    # the implementation's completed-square expressions.
    a = [float(x[0]) for x in packets(1, seed=seed)+packets(1, seed=seed+1)]
    q1, p1, g1, a1, q2, p2, g2, a2 = a
    radius = 12/np.sqrt(min(a1, a2))
    lo, hi = min(q1, q2)-radius, max(q1, q2)+radius
    nodes, weights = np.polynomial.legendre.leggauss(400)
    x = (hi+lo)/2+(hi-lo)/2*nodes
    w = (hi-lo)/2*weights
    bra = wave(x, *a[:4]).conj()
    ket = wave(x, *a[4:])
    d = -2*a2*(x-q2)+1j*p2
    for name, factor in [('overlap', 1), ('dipole', x), ('coupling', d),
                         ('kinetic', -0.5*(d*d-2*a2))]:
        close(getattr(gwp, 'gwp_'+name)(*[backend(v) for v in a]),
              np.sum(w*bra*factor*ket))
    scale, q, p, g, alp = gwp.gwp_product_decomposition(*[backend(v) for v in a])
    close(wave(x, *a[:4])*ket, array(scale)*wave(x, *map(array, (q, p, g, alp))))
    third = [0.4, -0.3, 0.7, 1.1]
    close(gwp.gwp_overlap(*[backend(v) for v in a+third]),
          np.sum(w*bra*ket.conj()*wave(x, *third)))


def test_symmetries_and_diagonal(backend):
    """Check normalization, diagonal limits and operator exchange symmetries.

    Overlap, position and kinetic operators are Hermitian; the spatial
    derivative is anti-Hermitian. A separate check converts a total legacy
    phase into per-DOF phases at hbar=1.9 and compares both overloads.
    """
    a, b = packets(), packets(seed=910)
    aa, bb = [backend(x) for x in a], [backend(x) for x in b]
    close(gwp.gwp_overlap(*aa, *aa), 1)
    close(gwp.gwp_dipole(*aa, *aa), a[0])
    close(gwp.gwp_coupling(*aa, *aa), 1j*a[1])
    close(gwp.gwp_kinetic(*aa, *aa), np.sum(a[3]+a[1]**2)/2)
    for name in ['overlap', 'dipole', 'kinetic', 'coupling']:
        f = getattr(gwp, 'gwp_'+name)
        sign = -1 if name == 'coupling' else 1
        close(f(*aa, *bb), sign*array(f(*bb, *aa)).conj())
    # General and equal-width/hbar overloads must agree in multiple dimensions.
    n, hbar, width = 3, 1.9, 0.8
    general = [a[0], a[1]/hbar, np.full(n, 0.2/(n*hbar)), np.full(n, width),
               b[0], b[1]/hbar, np.full(n, -0.5/(n*hbar)), np.full(n, width)]
    legacy = [a[0], a[1], 0.2, b[0], b[1], -0.5, width, hbar]
    for name in ['overlap', 'dipole', 'coupling', 'kinetic']:
        f = getattr(gwp, 'gwp_'+name)
        close(f(*map(backend, general)), f(*map(backend, legacy)))


def check_jacobian(function, inputs, step=1e-6, rtol=1e-3, atol=1e-5):
    """Compare every real/imaginary output derivative to centered differences.

    Inputs are real float64 tensors. Tuple outputs are supported. Calling
    autograd on each real scalar output with an implicit unit seed avoids
    the symbolic-shape/SymPy import triggered by gradcheck's explicit seeds.
    No PyTorch internals are patched, and no derivative checks are skipped.
    """
    import torch

    def components(args):
        outputs = function(*args)
        if not isinstance(outputs, tuple):
            outputs = (outputs,)
        parts = []
        for output in outputs:
            parts.append(output.real.reshape(-1))
            if output.is_complex():
                parts.append(output.imag.reshape(-1))
        return torch.cat(parts)

    values = components(inputs)
    assert values.requires_grad, "outputs are detached from autograd"
    analytic = []
    for value in values:
        gradients = torch.autograd.grad(value, inputs, retain_graph=True, allow_unused=True)
        analytic.append(np.concatenate([
            np.zeros(x.numel()) if g is None else array(g).reshape(-1)
            for x, g in zip(inputs, gradients)
        ]))

    numeric = []
    # Each perturbation is evaluated on fresh tensors outside the graph.
    with torch.no_grad():
        for arg_index, arg in enumerate(inputs):
            for element in range(arg.numel()):
                plus = [x.detach().clone() for x in inputs]
                minus = [x.detach().clone() for x in inputs]
                plus[arg_index].reshape(-1)[element] += step
                minus[arg_index].reshape(-1)[element] -= step
                numeric.append(array((components(plus)-components(minus))/(2*step)))
    np.testing.assert_allclose(np.asarray(analytic), np.asarray(numeric).T,
                               rtol=rtol, atol=atol)


def test_jacobian_checker():
    """Accept correct complex/tuple derivatives and reject a detached term."""
    torch = pytest.importorskip('torch')
    x = torch.tensor([0.3, -0.7], dtype=torch.float64, requires_grad=True)
    check_jacobian(lambda x: (x**2, 1j*x), (x,))
    # Values remain correct, but detaching x removes its nonzero derivative.
    with pytest.raises(AssertionError):
        check_jacobian(lambda x: x.detach()+0*x, (x,))


@pytest.mark.parametrize('name', ['overlap', 'dipole', 'coupling', 'kinetic',
                                  'product_decomposition'])
def test_autograd(name):
    """Compare PyTorch Jacobians with finite differences for all pair parameters.

    Float64 Jacobian comparisons check real and imaginary outputs for overlap, dipole,
    coupling and kinetic integrals, and all returned decomposition values.
    This catches detached tensors or incorrect derivatives, not only values.
    """
    torch = pytest.importorskip('torch')
    args = tuple(torch.tensor(x, dtype=torch.float64, requires_grad=True)
                 for x in packets(2)+packets(2, seed=10))
    check_jacobian(getattr(gwp, 'gwp_'+name), args)


@pytest.mark.parametrize('name', ['overlap_matrix', 'kinetic_matrix'])
def test_matrix_autograd(name):
    """Check gradients through broadcast rectangular pairwise matrices.

    Different packet counts exercise broadcasting. Kinetic checks also
    differentiate the inverse masses, which are shared by every packet pair.
    """
    torch = pytest.importorskip('torch')
    args = tuple(torch.tensor(x, dtype=torch.float64, requires_grad=True)
                 for x in packets(2, 2)+packets(2, 3, seed=10))
    if name == 'kinetic_matrix':
        args += (torch.tensor([0.7, 1.2], dtype=torch.float64, requires_grad=True),)
    check_jacobian(getattr(gwp, 'gwp_'+name), args)


@pytest.mark.parametrize('density', [False, True])
def test_value_gradient(backend, density):
    """Compare spatial derivatives with centered finite differences in r.

    Test both the complex wavefunction and real probability density. The
    step is 1e-5; rtol=1e-8 and atol=1e-10 allow finite-difference error,
    unlike the tighter direct-formula comparison tolerances.
    """
    q, p, _, alp = packets()
    r = np.array([0.5, -0.3, 0.4])
    tail = [q, alp] if density else [q, p, alp]
    numerical = []
    for i in range(3):
        dr = np.eye(3)[i]*1e-5
        numerical.append((gwp.gwp_value(r+dr, *tail)-gwp.gwp_value(r-dr, *tail))/2e-5)
    np.testing.assert_allclose(array(gwp.gwp_deriv(*map(backend, [r]+tail))), numerical,
                               rtol=1e-8, atol=1e-10)


def test_invalid_inputs(backend):
    """Check rejection of bad widths, shapes, inverse masses and zero hbar.

    Width cases include zero, negative, NaN and infinity. These checks use
    the Python implementation only; some invalid C++ calls terminate the
    process rather than raising an exception.
    """
    a = list(map(backend, packets()+packets(seed=7)))
    for width in [0., -1., np.nan, np.inf]:
        bad = a.copy()
        bad[3] = backend(width)
        with pytest.raises(ValueError, match='positive'):
            gwp.gwp_overlap(*bad)
    bad = a.copy()
    bad[4] = backend([0., 1.])
    with pytest.raises(ValueError, match='dimensions'):
        gwp.gwp_overlap(*bad)
    with pytest.raises(ValueError, match='inverse mass'):
        gwp.gwp_kinetic(*a, backend([1., 2.]))
    with pytest.raises(ValueError, match='hbar'):
        gwp.gwp_overlap(a[0], a[1], 0., a[4], a[5], 0., 1., 0.)
    with pytest.raises(ValueError, match='shape'):
        gwp.gwp_overlap_matrix(*a)


@pytest.mark.parametrize('dtype', ['float32', 'float64'])
def test_torch_dtype_device(dtype):
    """Check float32/64 map to complex64/128 without leaving the input device.

    CPU is always exercised when PyTorch is installed. CUDA is also exercised
    when available; a CPU-only run makes no claim of GPU execution coverage.
    """
    torch = pytest.importorskip('torch')
    for device in ['cpu']+(['cuda'] if torch.cuda.is_available() else []):
        args = [torch.tensor(x, dtype=getattr(torch, dtype), device=device)
                for x in packets()+packets(seed=512)]
        result = gwp.gwp_overlap(*args)
        assert result.device == args[0].device
        assert result.dtype == (torch.complex64 if dtype == 'float32' else torch.complex128)
