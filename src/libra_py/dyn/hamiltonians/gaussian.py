"""Assemble coherent Gaussian Hamiltonians from standard model/storage data.

Local electronic Hamiltonians and derivatives belong to TensorStorage and are
filled by HamiltonianEngine. This module assembles transient packet-pair S,H;
Gaussian overlap/kinetic integrals and potential moments live in dyn.gwp.
"""
from types import SimpleNamespace
from .engine import HamiltonianEngine
from .aux import as_result_mapping
from ..backends.arrays import asarray, xp_of, adjoint, indices
from ..gwp.storage import read_packets
from ..gwp.packets import spatial_overlap, spatial_kinetic
from ..gwp.potentials import taylor_ratio, gaussian_coupling_ratio


def gaussian_hamiltonian(storage, traj, model, method=0, exact_coupling=None):
    """Evaluate local model data and return active packet overlap/Hamiltonian.

    ``model`` follows HamiltonianEngine's ``model(R,P,storage,traj)`` contract;
    existing dyn.models objects work directly. Method 0 is BAT, 1 LHA,
    2 LHA with exact off-diagonal Gaussian coupling, 3 BAT with exact coupling.
    Methods 2/3 require ``exact_coupling=(A,B,C)`` representing the SAME
    off-diagonal potential sum_d A exp[-B*(r_d-C)^2] for each state pair.
    Local model matrices must be finite/Hermitian and LHA requires analytic
    Hessians. Labels may be unsorted, counts unequal and slots noncontiguous.
    """
    if method not in (0, 1, 2, 3):
        raise ValueError('potential method must be 0, 1, 2 or 3')
    if method in (2, 3) and exact_coupling is None:
        raise ValueError('exact coupling methods require (A,B,C)')
    slots, (q, p, a, s), labels, c, mass = read_packets(storage, traj)
    active = SimpleNamespace(id=traj.id, tbf_ids=slots, metadata=traj.metadata)
    n, d, ns = len(slots), storage.ndof, storage.nstates

    def checked_model(**kwargs):
        result = dict(as_result_mapping(model(**kwargs)))
        expected = {'H_dia': (n, ns, ns), 'dH_dia': (n, d, ns, ns)}
        if method in (1, 2):
            expected['d2H_dia'] = (n, d, d, ns, ns)
        for name, shape in expected.items():
            if result.get(name) is None:
                raise ValueError(f'Gaussian method {method} requires model {name}')
            value = asarray(result[name], q, True)
            if value.shape != shape:
                raise ValueError(f'{name} must have shape {shape}')
            xp = xp_of(value)
            if not bool(xp.isfinite(value).all()) or not bool(xp.allclose(value, adjoint(value), rtol=1e-6, atol=1e-8)):
                raise ValueError(f'{name} must be finite and Hermitian')
            result[name] = value
        return result

    HamiltonianEngine(storage.backend).evaluate(active, storage, checked_model, rep='diabatic')
    h, dh = storage.ham_dia[traj.id, slots], storage.d1ham_dia[traj.id, slots]
    hh = storage.d2ham_dia[traj.id, slots] if method in (1, 2) else None
    st = indices(labels, q)
    mask = st[:, None] == st[None, :]
    spatial = spatial_overlap(q, p, a, s)
    S = spatial*mask
    H = spatial_kinetic(q, p, a, s, mass)*mask
    for n1 in sorted(set(labels)):
        i = indices([k for k, label in enumerate(labels) if label == n1], q)
        for n2 in sorted(set(labels)):
            j = indices([k for k, label in enumerate(labels) if label == n2], q)
            args = (q[:, i], p[:, i], a[:, i], q[:, j], p[:, j], a[:, j])
            if method in (2, 3) and n1 != n2:
                v = gaussian_coupling_ratio(*args, exact_coupling)
            else:
                v = taylor_ratio(*args, h[i, n1, n2], h[j, n1, n2],
                                 dh[i, :, n1, n2], dh[j, :, n1, n2],
                                 None if hh is None else hh[i, :, :, n1, n2],
                                 None if hh is None else hh[j, :, :, n1, n2])
            H[i[:, None], j[None, :]] = H[i[:, None], j[None, :]]+v*spatial[i[:, None], j[None, :]]
    return S, H
