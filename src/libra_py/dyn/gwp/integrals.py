# Copyright (C) 2026 Alexey V. Akimov
# SPDX-License-Identifier: GPL-3.0-or-later

"""Normalized Gaussian wavepacket integrals for NumPy and PyTorch.

Functions accept the positional overloads of ``liblibra_core`` (see README.md).
Passing any torch.Tensor selects PyTorch; otherwise results are NumPy arrays or
scalars. Torch calculations preserve device, floating dtype and autograd graphs.
No C++ extension is needed by this module.

In atomic units a separable normalized packet is
``G(r) = product_i (2*alp_i/pi)**(1/4) *
exp(sum_i(-alp_i*(r_i-q_i)**2 + 1j*p_i*(r_i-q_i) + 1j*gamma_i))``.
Widths alp have units of inverse length squared. The probability density's
coordinate variance is 1/(4*alp_i): larger widths mean narrower packets.
The general interface uses phases per coordinate; the legacy equal-width
interface instead uses a total phase divided by hbar. All inputs are positional.

Single-packet integrals contract all DOFs. Matrix interfaces take one packet
per column and return all pairs. See each function for its overloads, output
shape, and whether it computes a wavefunction, density, or operator integral.
"""
import math
import numpy as np

__all__ = ["gwp_value", "gwp_deriv", "gwp_product_decomposition", "gwp_overlap",
           "gwp_dipole", "gwp_coupling", "gwp_kinetic", "gwp_overlap_matrix",
           "gwp_kinetic_matrix"]


def _arrays(args):
    """Select NumPy/PyTorch and coerce real parameters without detaching tensors."""
    # Lazy import keeps NumPy usable without the optional torch dependency.
    import sys
    torch = sys.modules.get("torch")
    tensors = [] if torch is None else [x for x in args if isinstance(x, torch.Tensor)]
    if tensors:
        ref = tensors[0]
        dtype = ref.dtype
        for x in tensors:
            if x.device != ref.device:
                raise ValueError("all tensors must be on the same device")
            dtype = torch.promote_types(dtype, x.dtype)
        if dtype.is_complex:
            raise ValueError("Gaussian parameters must be real")
        if not dtype.is_floating_point:
            dtype = torch.float64
        return torch, [torch.as_tensor(x, dtype=dtype, device=ref.device) for x in args]
    arrays = [np.asarray(x) for x in args]
    if any(np.iscomplexobj(x) for x in arrays):
        raise ValueError("Gaussian parameters must be real")
    dtype = np.result_type(*[x.dtype for x in arrays], np.float32)
    return np, [np.asarray(x, dtype=dtype) for x in arrays]


def _positive(x, name):
    """Reject nonpositive or nonfinite entries in a named parameter."""
    if not bool((x > 0).all()) or not bool((x < float("inf")).all()):
        raise ValueError(name + " must be finite and positive")


def _vector(x):
    """Validate a scalar/vector/column input and flatten it to one DOF axis."""
    if x.ndim > 2 or (x.ndim == 2 and x.shape[1] != 1):
        raise ValueError("expected a scalar, (ndof,) or (ndof, 1) vector")
    if any(n == 0 for n in x.shape):
        raise ValueError("empty vectors are not supported")
    return x.reshape(-1)


def _packets(args):
    """Validate eight general packet parameters, allowing scalar widths/phases."""
    xp, a = _arrays(args)
    a = [_vector(x) for x in a]
    n = a[0].shape[0]
    if any(x.shape[0] != n for x in (a[1], a[4], a[5])):
        raise ValueError("coordinate and momentum dimensions must match")
    if any(x.shape[0] not in (1, n) for x in a):
        raise ValueError("phase and width dimensions must match coordinates")
    _positive(a[3], "alp1")
    _positive(a[7], "alp2")
    return xp, a


def _parameters(args):
    """Convert general or legacy pair inputs into atomic-unit packet vectors."""
    if len(args) != 8:
        raise TypeError("expected eight packet parameters")
    xp, a = _arrays(args)
    # Legacy equal-width overload: R1,P1,gamma1,R2,P2,gamma2,alp,hbar.
    if a[2].ndim == 0 and a[3].ndim > 0:
        q1, p1, g1, q2, p2, g2, alp, hbar = a
        _positive(hbar, "hbar")
        if hbar.ndim or alp.ndim or g2.ndim:
            raise ValueError("legacy phases, width and hbar must be scalars")
        n = _vector(q1).shape[0]
        return _packets((q1, p1/hbar, g1/(hbar*n), alp,
                         q2, p2/hbar, g2/(hbar*n), alp))
    return _packets(a)


def _moments(xp, a):
    """Return overlap and the position/derivative prefactors from Gaussian moments.

    Axis 0 is the DOF axis; trailing axes may broadcast over packet pairs.
    Returns (S, mu, D, b), where dipole = S*mu, coupling = S*D, and the
    per-DOF kinetic factor is (b-D**2)/2. mu is the complex Gaussian mean,
    D = 2*a1*a2/(a1+a2)*(q2-q1) + i*(a2*p1+a1*p2)/(a1+a2), and
    b = 2*a1*a2/(a1+a2). S contracts the DOF axis. No division by S is
    used, so overlap underflow does not cause division by zero.
    """
    q1, p1, g1, a1, q2, p2, g2, a2 = a
    width = a1 + a2
    b = 2*a1*a2/width
    dq, dp = q2-q1, p2-p1
    momentum = (a2*p1+a1*p2)/width
    log_s = (0.25*xp.log(4*a1*a2/width**2) - 0.5*b*dq**2
             - dp**2/(4*width) + 1j*(g2-g1-momentum*dq))
    s = xp.exp(log_s.sum(axis=0))
    mean = (a1*q1+a2*q2)/width + 0.5j*dp/width
    derivative = b*dq + 1j*momentum
    return s, mean, derivative, b


def gwp_overlap(*args):
    """Integrate the overlap of two normalized Gaussian wavepackets.

    Computes ``S12 = integral conj(G1(r))*G2(r) dr`` over all coordinates.
    For identical packets S12 = 1; swapping packets conjugates S12.

    Parameters
    ----------
    *args : real scalars, array_like, or torch.Tensor
        General order: ``q1, p1, gamma1, alp1, q2, p2, gamma2, alp2``.
        ``q`` is the packet center, ``p`` its momentum, ``gamma`` its phase
        per coordinate, and ``alp`` its finite, positive Gaussian width.
        Coordinates and momenta have matching lengths: scalar for one DOF,
        ``(ndof,)``, or ``(ndof, 1)``. Widths and phases have one entry per
        DOF or a single broadcast entry. Atomic units (hbar = 1) apply.
        Use vector phases when using vector widths: a scalar third argument
        and an array fourth argument select the legacy overload instead.

        Legacy order: ``R1, P1, gamma1, R2, P2, gamma2, alp, hbar``.
        Here R/P are vectors, phases are scalar total phases (not per DOF),
        and both packets share a scalar width and positive scalar hbar.

        A third packet can be appended in the general convention, giving
        twelve arguments ``q1,p1,gamma1,alp1,...,q3,p3,gamma3,alp3``. This
        computes ``integral conj(G1*G2)*G3 dr``: BOTH first packets are
        conjugated, matching the C++ implementation.

    Returns
    -------
    complex scalar or zero-dimensional torch.Tensor
        The integrated overlap, with all spatial dimensions contracted.

    Raises
    ------
    TypeError
        The argument count does not match a supported overload.
    ValueError
        Packet shapes, widths, hbar, real-valued parameters, or tensor
        devices violate the input conventions.

    Notes
    -----
    NumPy inputs produce NumPy results. Any tensor parameter selects PyTorch
    with promoted tensor dtype and a shared device; autograd is preserved.
    See the module docstring for normalization and phase conventions.
    """
    if len(args) == 12:
        scale, q, p, gamma, alp = gwp_product_decomposition(*args[:8])
        return scale*gwp_overlap(q, p, gamma, alp, *args[8:])
    xp, a = _parameters(args)
    return _moments(xp, a)[0]


def _component_result(result, first):
    """Restore a component-wise result to the first coordinate input shape."""
    shape = getattr(first, "shape", np.shape(first))
    return result.reshape(shape)


def gwp_dipole(*args):
    """Integrate the position operator between two wavepackets.

    Parameters
    ----------
    *args : real scalars, array_like, or torch.Tensor
        General order: ``q1, p1, gamma1, alp1, q2, p2, gamma2, alp2``.
        ``q`` is the packet center, ``p`` its momentum, ``gamma`` its phase
        per coordinate, and ``alp`` its finite, positive Gaussian width.
        Coordinates and momenta have matching lengths: scalar for one DOF,
        ``(ndof,)``, or ``(ndof, 1)``. Widths and phases have one entry per
        DOF or a single broadcast entry. Atomic units (hbar = 1) apply.
        Use vector phases when using vector widths: a scalar third argument
        and an array fourth argument select the legacy overload instead.

        Legacy order: ``R1, P1, gamma1, R2, P2, gamma2, alp, hbar``.
        Here R/P are vectors, phases are scalar total phases (not per DOF),
        and both packets share a scalar width and positive scalar hbar.

    Returns
    -------
    complex ndarray or torch.Tensor
        Component i is ``integral conj(G1)*r_i*G2 dr``.
        The result has the shape of q1 (R1 in the legacy overload);
        scalar input produces a zero-dimensional result.

    Raises
    ------
    TypeError
        The argument count does not match a supported overload.
    ValueError
        Packet shapes, widths, hbar, real-valued parameters, or tensor
        devices violate the input conventions.

    Notes
    -----
    NumPy inputs produce NumPy results. Any tensor parameter selects PyTorch
    with promoted tensor dtype and a shared device; autograd is preserved.
    See the module docstring for normalization and phase conventions.

    This is a position matrix element; no charge factor is included.
    For identical packets the result is the packet center q.
    """
    xp, a = _parameters(args)
    s, mean, _, _ = _moments(xp, a)
    return _component_result(s*mean, args[0])


def gwp_coupling(*args):
    """Integrate the spatial derivative acting on the second packet.

    Parameters
    ----------
    *args : real scalars, array_like, or torch.Tensor
        General order: ``q1, p1, gamma1, alp1, q2, p2, gamma2, alp2``.
        ``q`` is the packet center, ``p`` its momentum, ``gamma`` its phase
        per coordinate, and ``alp`` its finite, positive Gaussian width.
        Coordinates and momenta have matching lengths: scalar for one DOF,
        ``(ndof,)``, or ``(ndof, 1)``. Widths and phases have one entry per
        DOF or a single broadcast entry. Atomic units (hbar = 1) apply.
        Use vector phases when using vector widths: a scalar third argument
        and an array fourth argument select the legacy overload instead.

        Legacy order: ``R1, P1, gamma1, R2, P2, gamma2, alp, hbar``.
        Here R/P are vectors, phases are scalar total phases (not per DOF),
        and both packets share a scalar width and positive scalar hbar.

    Returns
    -------
    complex ndarray or torch.Tensor
        Component i is ``integral conj(G1)*dG2/dr_i dr``.
        The result has the shape of q1 (R1 in the legacy overload);
        scalar input produces a zero-dimensional result.

    Raises
    ------
    TypeError
        The argument count does not match a supported overload.
    ValueError
        Packet shapes, widths, hbar, real-valued parameters, or tensor
        devices violate the input conventions.

    Notes
    -----
    NumPy inputs produce NumPy results. Any tensor parameter selects PyTorch
    with promoted tensor dtype and a shared device; autograd is preserved.
    See the module docstring for normalization and phase conventions.

    The derivative is with respect to the integration coordinate r, not
    a packet parameter. No momentum factor -i*hbar is included. Swapping
    packets gives the negative complex conjugate; a diagonal element is
    i*p in atomic units.
    """
    xp, a = _parameters(args)
    s, _, derivative, _ = _moments(xp, a)
    return _component_result(s*derivative, args[0])


def gwp_kinetic(*args):
    """Integrate the kinetic operator with diagonal inverse masses.

    Computes ``<G1| -1/2 sum_i invM_i * d_i^2 |G2>``. Derivatives act
    on G2; off-diagonal elements may be complex, even for real masses.

    Parameters
    ----------
    *args : real scalars, array_like, or torch.Tensor
        General order: ``q1, p1, gamma1, alp1, q2, p2, gamma2, alp2``.
        ``q`` is the packet center, ``p`` its momentum, ``gamma`` its phase
        per coordinate, and ``alp`` its finite, positive Gaussian width.
        Coordinates and momenta have matching lengths: scalar for one DOF,
        ``(ndof,)``, or ``(ndof, 1)``. Widths and phases have one entry per
        DOF or a single broadcast entry. Atomic units (hbar = 1) apply.
        Use vector phases when using vector widths: a scalar third argument
        and an array fourth argument select the legacy overload instead.

        Legacy order: ``R1, P1, gamma1, R2, P2, gamma2, alp, hbar``.
        Here R/P are vectors, phases are scalar total phases (not per DOF),
        and both packets share a scalar width and positive scalar hbar.

        An optional ninth argument ``invM`` supplies inverse masses as a
        scalar, ``(ndof,)``, or ``(ndof, 1)`` array. The default is one
        for each coordinate. This is not a full inverse-mass matrix.

    Returns
    -------
    complex scalar or zero-dimensional torch.Tensor
        Sum of the kinetic matrix elements over all DOFs. In atomic units
        this is energy for the supplied inverse masses.

    Raises
    ------
    TypeError
        The argument count does not match a supported overload.
    ValueError
        Packet shapes, widths, hbar, real-valued parameters, or tensor
        devices violate the input conventions.

    Notes
    -----
    NumPy inputs produce NumPy results. Any tensor parameter selects PyTorch
    with promoted tensor dtype and a shared device; autograd is preserved.
    See the module docstring for normalization and phase conventions.

    The legacy overload retains the C++ operator -1/2 Laplacian when
    invM is omitted, even for hbar != 1. Multiply that result by
    hbar**2/mass for physical kinetic energy with a common mass. If invM
    is supplied, multiply by hbar**2 only; masses are already included.
    """
    xp, a = _parameters(args[:8])
    if len(args) not in (8, 9):
        raise TypeError("expected eight parameters and optional inverse masses")
    mass = 1
    if len(args) == 9:
        xp, converted = _arrays((*a, args[8]))
        a, mass = converted[:8], _vector(converted[8])
        if mass.shape[0] not in (1, a[0].shape[0]):
            raise ValueError("inverse mass dimensions must match coordinates")
    s, _, derivative, b = _moments(xp, a)
    return 0.5*s*((b-derivative**2)*mass).sum(axis=0)


def gwp_product_decomposition(*args):
    """Express the pointwise product of two packets as a scaled packet.

    Neither packet is conjugated: ``G1(r)*G2(r) = scale*G(r)``.

    Parameters
    ----------
    *args : real scalars, array_like, or torch.Tensor
        Exactly ``q1,p1,gamma1,alp1,q2,p2,gamma2,alp2``, in atomic units.
        q/p have matching scalar, ``(ndof,)`` or ``(ndof, 1)`` shapes.
        gamma is a phase per DOF; alp is a finite positive width per DOF.
        Widths and phases may each contain a single broadcast entry.
        There is no legacy equal-width/hbar overload for this function.

    Returns
    -------
    scale : real scalar or zero-dimensional torch.Tensor
        Normalization prefactor, including the pi**(-ndof/4) factor.
    q, p, gamma, alp : ndarray or torch.Tensor
        Parameters of the normalized product packet, each with q1's
        shape. The width is alp1+alp2, momentum is p1+p2, and center is
        weighted by the input widths. Phases include the shift of center.
        Inputs are not mutated; outputs are returned as a five-item tuple.

    Raises
    ------
    TypeError
        Not exactly eight inputs were supplied.
    ValueError
        Packet dimensions, widths, parameter types, or devices are invalid.

    Notes
    -----
    NumPy inputs produce NumPy results. Any tensor parameter selects PyTorch
    with promoted tensor dtype and a shared device; autograd is preserved.
    See the module docstring for normalization and phase conventions.
    """
    if len(args) != 8:
        raise TypeError("expected eight packet parameters")
    xp, a = _packets(args)
    q1, p1, g1, a1, q2, p2, g2, a2 = a
    width = a1+a2
    q = (a1*q1+a2*q2)/width
    gamma = g1+g2+(q2-q1)*(a2*p1-a1*p2)/width
    log_scale = 0.25*xp.log(2*a1*a2/(math.pi*width))-a1*a2/width*(q2-q1)**2
    scale = xp.exp(log_scale.sum())
    # Broadcast constant widths/phases to every coordinate before returning.
    return (scale, *[_component_result(x+0*q, args[0]) for x in (q, p1+p2, gamma, width)])


def _value_parameters(args):
    """Normalize value/gradient overloads and report whether density is requested."""
    xp, a = _arrays(args)
    if len(a) == 6:
        r, q, p, gamma, alp, hbar = a
        _positive(hbar, "hbar")
        if any(x.ndim for x in (gamma, alp, hbar)):
            raise ValueError("gamma, alp, hbar must be scalars in the six-argument overload")
        gamma = gamma/(hbar*_vector(q).shape[0])
        p = p/hbar
    elif len(a) == 4:
        r, q, p, alp = a
        gamma = 0
    elif len(a) == 3:
        r, q, alp = a
        p, gamma = 0*q, 0
    else:
        raise TypeError("expected three, four or six arguments")
    xp, a = _packets((q, p, gamma, alp, r, p, gamma, alp))
    return xp, a, len(args) == 3


def gwp_value(*args):
    """Evaluate a normalized wavefunction or its probability density.

    Parameters
    ----------
    *args : real scalars, array_like, or torch.Tensor
        Supported positional overloads:

        * ``r, q, p, alp``: wavefunction with zero phase, hbar = 1.
        * ``r, q, alp``: probability density |G(r)|**2, independent of
          momentum/phase. This is NOT the zero-momentum wavefunction.
        * ``r, q, p, gamma, alp, hbar``: wavefunction with scalar total
          phase gamma, a common scalar width alp, and positive scalar hbar.

        r is the evaluation point and q the packet center; each is scalar,
        ``(ndof,)``, or ``(ndof, 1)`` with matching length. p has the same
        length. Widths must be finite and positive; the three/four-argument
        overloads accept either per-DOF widths or a broadcast scalar width.
        A vector r represents one multidimensional point, not a batch.

    Returns
    -------
    scalar or zero-dimensional torch.Tensor
        Complex wavefunction value for four/six arguments, or real density
        for three. The wavefunction has integral |G|**2 = 1; the density
        itself integrates to 1 over all space.

    Raises
    ------
    TypeError
        The argument count does not match a supported overload.
    ValueError
        Packet shapes, widths, hbar, real-valued parameters, or tensor
        devices violate the input conventions.

    Notes
    -----
    NumPy inputs produce NumPy results. Any tensor parameter selects PyTorch
    with promoted tensor dtype and a shared device; autograd is preserved.
    See the module docstring for normalization and phase conventions.
    """
    xp, a, density = _value_parameters(args)
    q, p, gamma, alp, r = a[:5]
    log_g = 0.25*xp.log(2*alp/math.pi)-alp*(r-q)**2+1j*(p*(r-q)+gamma)
    if density:
        return xp.exp(2*log_g.real.sum())
    return xp.exp(log_g.sum())


def gwp_deriv(*args):
    """Evaluate the spatial gradient of a wavefunction or density.

    Parameters
    ----------
    *args : real scalars, array_like, or torch.Tensor
        The same overloads as :func:`gwp_value`: ``r,q,p,alp``, ``r,q,alp``,
        or ``r,q,p,gamma,alp,hbar``. r/q/p have matching scalar, ``(ndof,)``
        or ``(ndof, 1)`` lengths. Widths are positive; the six-argument
        overload requires scalar gamma/alp/hbar and uses a total phase.

    Returns
    -------
    ndarray or torch.Tensor
        Gradient with respect to r, with r's input shape. Wavefunction
        gradients are complex; density gradients are real. These are
        spatial derivatives, not gradients with respect to q or p.

    Raises
    ------
    TypeError
        The argument count does not match a supported overload.
    ValueError
        Packet shapes, widths, hbar, real-valued parameters, or tensor
        devices violate the input conventions.

    Notes
    -----
    NumPy inputs produce NumPy results. Any tensor parameter selects PyTorch
    with promoted tensor dtype and a shared device; autograd is preserved.
    See the module docstring for normalization and phase conventions.

    For a wavefunction, component i is
    ``(-2*alp_i*(r_i-q_i) + 1j*p_i/hbar)*G(r)``. For the density it is
    ``-4*alp_i*(r_i-q_i)*|G(r)|**2``. The six-argument derivative overload
    is a Python convenience beyond the C++ derivative overloads.
    """
    xp, a, density = _value_parameters(args)
    q, p, _, alp, r = a[:5]
    factor = -4*alp*(r-q) if density else -2*alp*(r-q)+1j*p
    return _component_result(factor*gwp_value(*args), args[0])


def _matrix_parameters(args):
    """Validate packet matrices and add pair axes for (ndof,ntraj1,ntraj2) broadcasting."""
    xp, a = _arrays(args)
    if any(x.ndim != 2 for x in a):
        raise ValueError("matrix parameters must have shape (ndof, ntraj)")
    if any(x.shape != a[0].shape for x in a[:4]) or any(x.shape != a[4].shape for x in a[4:]):
        raise ValueError("each packet set must have matching shapes")
    if a[0].shape[0] != a[4].shape[0] or any(0 in x.shape for x in a):
        raise ValueError("packet sets must have the same nonzero ndof and nonempty trajectories")
    _positive(a[3], "alp1")
    _positive(a[7], "alp2")
    return xp, [x[:, :, None] for x in a[:4]] + [x[:, None, :] for x in a[4:]]


def gwp_overlap_matrix(*args):
    """Compute every pairwise overlap between two packet sets.

    Parameters
    ----------
    *args : array_like or torch.Tensor
        Packet matrices in order
        ``q1,p1,gamma1,alp1,q2,p2,gamma2,alp2``.
        All four arrays in set 1 have shape ``(ndof, ntraj1)``; those in
        set 2 have shape ``(ndof, ntraj2)``. Rows are spatial DOFs; columns
        are packets. Both sets must be nonempty and share ndof, but their
        packet counts can differ. q/p/gamma/alp are centers, momenta,
        phases per DOF, and positive finite widths. Scalar broadcasting
        of packet parameters is not supported by the matrix interfaces.
        Atomic units (hbar = 1) apply.

        To include orthonormal electronic states, use ten arguments:
        ``q1,p1,gamma1,alp1,state1,q2,p2,gamma2,alp2,state2``.
        state1/state2 contain one integer label per packet, with shapes
        ``(ntraj1,)`` and ``(ntraj2,)``. Different labels give zero overlap.

    Returns
    -------
    complex ndarray or torch.Tensor, shape (ntraj1, ntraj2)
        Entry [i,j] is <G1_i|G2_j>, optionally multiplied by the electronic
        state Kronecker delta. Using the same set twice gives a Hermitian
        matrix with unit diagonal, up to floating-point error.

    Raises
    ------
    TypeError
        The argument count does not match a supported overload.
    ValueError
        Packet shapes, widths, hbar, real-valued parameters, or tensor
        devices violate the input conventions.

    Notes
    -----
    NumPy inputs produce NumPy results. Any tensor parameter selects PyTorch
    with promoted tensor dtype and a shared device; autograd is preserved.
    See the module docstring for normalization and phase conventions.

    Temporary storage scales as ndof*ntraj1*ntraj2. For large bases, call
    this function on blocks of columns to limit memory use.
    """
    if len(args) == 10:
        packets = args[:4]+args[5:9]
        states = (args[4], args[9])
    elif len(args) == 8:
        packets, states = args, None
    else:
        raise TypeError("expected eight parameters, or ten including electronic states")
    xp, a = _matrix_parameters(packets)
    result = _moments(xp, a)[0]
    if states is not None:
        if xp is np:
            s1, s2 = [np.asarray(s) for s in states]
        else:
            s1, s2 = [xp.as_tensor(s, device=result.device) for s in states]
        if s1.shape != (result.shape[0],) or s2.shape != (result.shape[1],):
            raise ValueError("one electronic state is required per trajectory")
        result = result*(s1[:, None] == s2[None, :])
    return result


def gwp_kinetic_matrix(*args):
    """Compute pairwise kinetic matrix elements for two packet sets.

    Parameters
    ----------
    *args : array_like or torch.Tensor
        Packet matrices in order
        ``q1,p1,gamma1,alp1,q2,p2,gamma2,alp2``.
        All four arrays in set 1 have shape ``(ndof, ntraj1)``; those in
        set 2 have shape ``(ndof, ntraj2)``. Rows are spatial DOFs; columns
        are packets. Both sets must be nonempty and share ndof, but their
        packet counts can differ. q/p/gamma/alp are centers, momenta,
        phases per DOF, and positive finite widths. Scalar broadcasting
        of packet parameters is not supported by the matrix interfaces.
        Atomic units (hbar = 1) apply.

        A required ninth argument ``invM`` supplies diagonal inverse
        masses, shared by all packets: scalar, ``(ndof,)`` or ``(ndof, 1)``.
        Use ones for unit masses. No electronic-state overload is provided.

    Returns
    -------
    complex ndarray or torch.Tensor, shape (ntraj1, ntraj2)
        Entry [i,j] is <G1_i|-1/2 sum_k invM_k*d_k**2|G2_j>.
        Using the same set twice gives a Hermitian kinetic matrix.

    Raises
    ------
    TypeError
        The argument count does not match a supported overload.
    ValueError
        Packet shapes, widths, hbar, real-valued parameters, or tensor
        devices violate the input conventions.

    Notes
    -----
    NumPy inputs produce NumPy results. Any tensor parameter selects PyTorch
    with promoted tensor dtype and a shared device; autograd is preserved.
    See the module docstring for normalization and phase conventions.

    Temporary storage scales as ndof*ntraj1*ntraj2. Inverse masses are
    coordinate-dependent diagonal coefficients, not a dense mass tensor.
    """
    if len(args) != 9:
        raise TypeError("expected eight packet matrices and inverse masses")
    xp, converted = _arrays(args)
    xp, a = _matrix_parameters(converted[:8])
    mass = _vector(converted[8])
    if mass.shape[0] not in (1, a[0].shape[0]):
        raise ValueError("inverse mass dimensions must match coordinates")
    s, _, derivative, b = _moments(xp, a)
    return 0.5*s*((b-derivative**2)*mass[:, None, None]).sum(axis=0)
