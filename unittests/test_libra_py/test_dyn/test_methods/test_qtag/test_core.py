# Copyright (C) 2022-2026 Matthew Dutra, Alexey V. Akimov
# SPDX-License-Identifier: GPL-3.0-or-later

"""C++ parity and independent formula checks for QTAG's array primitives.

Rebuild libra_core before C++ comparisons. NumPy tests do not require Torch;
Torch tests use scalar autograd and run without SymPy. Formula references
include mixed-coordinate polynomials, unequal widths and nonzero phases.
"""
from types import SimpleNamespace
import numpy as np
import pytest
from libra_py.dyn.methods.qtag import core
from libra_py.dyn.methods.qtag._backend import asarray, zeros, xp_of


@pytest.fixture(params=["numpy", "torch"])
def backend(request):
    if request.param == "torch":
        torch = pytest.importorskip("torch")
        return lambda x: torch.as_tensor(x, dtype=torch.float64)
    return lambda x: np.asarray(x, dtype=float)


@pytest.fixture
def cpp():
    import sys
    return pytest.importorskip("cyglibra_core" if sys.platform == "cygwin" else "liblibra_core")


def numpy(x):
    return x.detach().cpu().resolve_conj().resolve_neg().numpy() if hasattr(x, "detach") else np.asarray(x)


def close(x, y, rtol=1e-10, atol=1e-11):
    np.testing.assert_allclose(numpy(x), numpy(y), rtol=rtol, atol=atol)


def packet_data():
    # Two DOFs, four trajectories, unequal widths and phases.
    rng = np.random.default_rng(318)
    return (rng.normal(size=(2, 4)), rng.normal(size=(2, 4)),
            rng.uniform(.7, 2.1, size=(2, 4)), rng.normal(size=(2, 4)))


def cm(cpp, x, complex=False):
    x = numpy(x)
    if x.ndim == 1:
        x = x[:, None]
    obj = (cpp.CMATRIX if complex else cpp.MATRIX)(*x.shape)
    for i, j in np.ndindex(x.shape):
        obj.set(i, j, x[i, j].item())
    return obj


def unmatrix(x):
    return np.array([[x.get(i, j) for j in range(x.num_of_cols)] for i in range(x.num_of_rows)])


def intlist(cpp, values):
    obj = cpp.intList()
    for n in values:
        obj.append(n)
    return obj


def polynomial(q, params=None):
    """A 2D quartic potential with variable Hessian and mixed derivatives."""
    x, y = q
    H = zeros((q.shape[1], 2, 2), q, True)
    dH = zeros((q.shape[1], 2, 2, 2), q, True)
    hh = zeros((q.shape[1], 2, 2, 2, 2), q, True)
    w = asarray([[1., .2+.1j], [.2-.1j, .7]], q, True)
    v = .3*x*x+.5*y*y+.2*x*y+.04*x**4
    H = v[:, None, None]*w
    H[:, 1, 1] = H[:, 1, 1]+.4
    dH[:, 0] = (.6*x+.2*y+.16*x**3)[:, None, None]*w
    dH[:, 1] = (y+.2*x)[:, None, None]*w
    hh[:, 0, 0] = (.6+.48*x*x)[:, None, None]*w
    hh[:, 1, 1] = 1.*w
    hh[:, 0, 1] = .2*w
    hh[:, 1, 0] = .2*w
    return core.HamiltonianData(H, dH, hh)


def cpp_model(cpp, q, params, full_id):
    values = polynomial(unmatrix(q))
    idx = list(full_id)[-1]
    d1, d2 = cpp.CMATRIXList(), cpp.CMATRIXList()
    for x in values.dH[idx]:
        d1.append(cm(cpp, x, True))
    for x in values.d2H[idx].reshape(4, 2, 2):
        d2.append(cm(cpp, x, True))
    return SimpleNamespace(ham_dia=cm(cpp, values.H[idx], True),
                           ovlp_dia=cm(cpp, np.eye(2), True),
                           d1ham_dia=d1, d2ham_dia=d2)


def test_cpp_elementary(cpp, backend):
    """Compare psi, overlaps, kinetic, embedding and corrected momentum."""
    q, p, a, s = map(backend, packet_data())
    args = [cm(cpp, x) for x in (q, p, a, s)]
    masses = backend([.3, 1.4])
    close(core.qtag_overlap_elementary(q, p, a, s), unmatrix(cpp.qtag_overlap_elementary(*args)))
    close(core.qtag_kinetic_elementary(q, p, a, s, masses),
          unmatrix(cpp.qtag_kinetic_elementary(*args, cm(cpp, masses))))
    c = np.array([[.2+.1j, .3, -.1j, .4], [-.3, .1j, .2, .1]])
    point = backend([.3, -.2])
    close(core.qtag_psi(point, q, p, a, s, c), unmatrix(cpp.qtag_psi(cm(cpp, point), *args, cm(cpp, c, True))))
    overlap = core.qtag_overlap_elementary(q, p, a, s)
    close(core.qtag_overlap([1, 0, 1, 0], overlap, 2),
          unmatrix(cpp.qtag_overlap(intlist(cpp, [1, 0, 1, 0]), cm(cpp, overlap, True), 2)))
    close(core.qtag_momentum(q, p, a, s, c[0]),
          unmatrix(cpp.qtag_momentum(*args, cm(cpp, c[0], True))))


@pytest.mark.parametrize("method", [0, 1, 2, 3])
def test_cpp_potential_and_assembly(cpp, backend, method):
    """Variable Hessians expose both C++ LHA endpoint and cross-term errors."""
    q, p, a, s = map(backend, packet_data())
    ham = cpp.nHamiltonian(2, 2, 2)
    ham.add_new_children(2, 2, 2, 4)
    ham.init_all(2, 1)
    callback = lambda q, params, full_id: cpp_model(cpp, q, params, full_id)
    ham.compute_diabatic(callback, cm(cpp, q), {}, 1)
    data = core.evaluate_model(polynomial, q)
    # The default C++ exact-coupling overload sets ABC to zero. Full assembly
    # below exercises nonzero Gaussian coupling read from model_params.
    for n1, n2 in ((0, 0), (0, 1)):
        i, j = [0, 1], [2, 3]
        pyargs = (q[:, i], p[:, i], s[:, i], a[:, i], n1, i,
                  q[:, j], p[:, j], s[:, j], a[:, j], n2, j)
        cargs = [cm(cpp, x) for x in pyargs[:4]]+[n1, intlist(cpp, i)]
        cargs += [cm(cpp, x) for x in pyargs[6:10]]+[n2, intlist(cpp, j)]
        close(core.qtag_potential(*pyargs, data, method), unmatrix(cpp.qtag_potential(*cargs, ham, method)))
    params = dict(ex_cpl_A=.25, ex_cpl_B=.6, ex_cpl_C=-.3)
    dyn = dict(qtag_pot_approx_method=method)
    S, H = core.qtag_hamiltonian_and_overlap(q, p, a, s, [0, 0, 1, 1], backend([.3, 1.4]), polynomial, params, dyn)
    cs, ch = cpp.CMATRIX(4, 4), cpp.CMATRIX(4, 4)
    cpp.qtag_hamiltonian_and_overlap(cm(cpp, q), cm(cpp, p), cm(cpp, a), cm(cpp, s),
                                     cpp.CMATRIX(2, 4), intlist(cpp, [0, 0, 1, 1]),
                                     cm(cpp, [.3, 1.4]), ham, callback, params, dyn, cs, ch)
    close(S, unmatrix(cs))
    close(H, unmatrix(ch))
    close(H, H.conj().T)
    # Same physical basis in interleaved order must give a permutation of H/S.
    order = [2, 0, 3, 1]
    ss, hh = core.qtag_hamiltonian_and_overlap(q[:, order], p[:, order], a[:, order], s[:, order],
                                              [1, 0, 1, 0], backend([.3, 1.4]), polynomial, params, dyn)
    close(ss, numpy(S)[np.ix_(order, order)])
    close(hh, numpy(H)[np.ix_(order, order)])
    cpp.qtag_hamiltonian_and_overlap(cm(cpp, q[:, order]), cm(cpp, p[:, order]),
                                     cm(cpp, a[:, order]), cm(cpp, s[:, order]),
                                     cpp.CMATRIX(2, 4), intlist(cpp, [1, 0, 1, 0]),
                                     cm(cpp, [.3, 1.4]), ham, callback, params, dyn, cs, ch)
    close(ss, unmatrix(cs))
    close(hh, unmatrix(ch))


def test_momentum_spatial_derivative(backend):
    """Independent wavefunction finite differences verify log-gradients/slopes."""
    q, p, a, s = map(backend, packet_data())
    c = asarray([.7, .2j, -.1, .3+.2j], q, True)
    ratio, slope, _ = core.momentum_fields(q, p, a, s, c)
    eps = 2e-4
    base = core.qtag_psi(q, q, p, a, s, c[None])[0]
    for d in range(2):
        shift = zeros(q.shape, q)
        shift[d] = eps
        plus = core.qtag_psi(q+shift, q, p, a, s, c[None])[0]
        minus = core.qtag_psi(q-shift, q, p, a, s, c[None])[0]
        deriv = (plus-minus)/(2*eps)
        close(ratio[d], deriv/base, rtol=2e-6, atol=2e-7)
        close(slope[d], (plus-2*base+minus)/eps**2/base-(deriv/base)**2, rtol=2e-5, atol=2e-6)


def test_quadrature_lha_and_exact_gaussian(backend):
    """Direct 2D integration checks mixed moments and exact coupling ratios."""
    q1, q2 = np.array([-.4, .2]), np.array([.5, -.3])
    p1, p2 = np.array([.1, -.2]), np.array([-.3, .4])
    a1, a2 = np.array([1.2, .8]), np.array([.9, 1.7])
    phase1, phase2 = np.array([.2, -.1]), np.array([-.3, .2])
    data = polynomial(np.stack([q1, q2], axis=1))
    x, w = np.polynomial.legendre.leggauss(150)
    mesh = np.array(np.meshgrid(10*x, 10*x, indexing='ij')).reshape(2, -1)
    weights = (100*w[:, None]*w[None, :]).reshape(-1)
    def wave(q, p, a, s):
        dq = mesh-q[:, None]
        return np.prod((a[:, None]/np.pi)**.25*np.exp(-a[:, None]*dq**2/2+1j*(p[:, None]*dq+s[:, None])), axis=0)
    integrand = wave(q1, p1, a1, phase1).conj()*wave(q2, p2, a2, phase2)
    overlap = np.sum(weights*integrand)
    # Average the two quadratic Taylor polynomials of a quartic potential.
    potential = np.zeros(mesh.shape[1], dtype=complex)
    for k, q in enumerate((q1, q2)):
        dq = mesh-q[:, None]
        potential += .5*data.H[k, 0, 0]+.5*np.einsum('d,dm->m', data.dH[k, :, 0, 0], dq)
        potential += .25*np.einsum('de,dm,em->m', data.d2H[k, :, :, 0, 0], dq, dq)
    args = [backend(x) for x in (q1, p1, phase1, a1, q2, p2, phase2, a2)]
    val = core.LHA(data.H[0], data.H[1], data.dH[0], data.dH[1], data.d2H[0], data.d2H[1],
                    *args[:4], 0, *args[4:], 0)
    close(val*overlap, np.sum(weights*integrand*potential))
    exact = .4*np.exp(-.7*(mesh+.2)**2).sum(axis=0)
    data_backend = core.evaluate_model(polynomial, backend(np.stack([q1, q2], axis=1)))
    for func in (core.BATe, core.LHAe):
        ratio = func(0, 1, *args[:4], 0, *args[4:], 1, .4, .7, -.2, data_backend)
        close(ratio*overlap, np.sum(weights*integrand*exact))


def test_cpp_coefficient_propagation(cpp, backend):
    """Generalized Hermitian evolution conserves the nonorthogonal norm."""
    q, p, a, s = map(backend, packet_data())
    S, H = core.qtag_hamiltonian_and_overlap(q, p, a, s, [0, 0, 1, 1], backend([.3, 1.4]), polynomial)
    c = asarray([.7, .2j, -.1, .3+.2j], q, True)
    actual = core.propagate_electronic_qtag(.07, c, H, S)
    cc = cm(cpp, c, True)
    cpp.propagate_electronic_qtag(.07, cc, cm(cpp, H, True), cm(cpp, S, True))
    close(actual, unmatrix(cc)[:, 0])
    close(actual.conj() @ S @ actual, c.conj() @ S @ c)
    close(core.propagate_electronic_qtag(-.07, actual, H, S), c)


def test_torch_scalar_gradients():
    """Finite-difference check of potential/kinetic gradients without SymPy."""
    torch = pytest.importorskip('torch')
    data = packet_data()
    args = [torch.tensor(x, dtype=torch.float64, requires_grad=True) for x in data]
    def objective(q, p, a, s):
        S, H = core.qtag_hamiltonian_and_overlap(q, p, a, s, [0, 1, 0, 1],
                                                torch.tensor([.3, 1.4]), polynomial,
                                                dyn_params={'qtag_pot_approx_method': 1})
        return (abs(H)**2).sum()+abs(S).sum()
    value = objective(*args)
    value.backward()
    eps = 1e-6
    for k in range(4):
        plus, minus = [x.detach().clone() for x in args], [x.detach().clone() for x in args]
        plus[k][0, 1] += eps
        minus[k][0, 1] -= eps
        close(args[k].grad[0, 1], (objective(*plus)-objective(*minus))/(2*eps), rtol=2e-5, atol=2e-6)


def test_invalid_inputs(backend):
    q, p, a, s = map(backend, packet_data())
    with pytest.raises(ValueError, match='positive'):
        core.qtag_overlap_elementary(q, p, -a, s)
    with pytest.raises(ValueError, match='label'):
        core.qtag_hamiltonian_and_overlap(q, p, a, s, [0, 2, 0, 0], backend([1, 1]), polynomial)
    def no_hessian(q, params):
        d = polynomial(q)
        return d.H, d.dH
    with pytest.raises(ValueError, match='Hessians'):
        core.qtag_hamiltonian_and_overlap(q, p, a, s, [0]*4, backend([1, 1]), no_hessian,
                                         dyn_params={'qtag_pot_approx_method': 1})
    with pytest.raises(ValueError, match='ex_cpl'):
        core.qtag_hamiltonian_and_overlap(q, p, a, s, [0]*4, backend([1, 1]), polynomial,
                                         dyn_params={'qtag_pot_approx_method': 2})
