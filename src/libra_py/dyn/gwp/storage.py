"""Access coherent, state-labeled Gaussian TBFs in shared TensorStorage.

Each active slot holds one packet on act_states_dia. Its coherent coefficient
is ampl_dia[slot, label], not a separately normalized electronic wavefunction.
Pair overlap/Hamiltonian matrices are transient, so storage lifecycle operations
only need to handle their usual single TBF axis. No parallel state is owned here.
"""
from ..backends.arrays import copy, indices, scalar, xp_of
from .packets import validate_basis, validate_states


def active_slots(storage, traj):
    """Return live slots in trajectory order, including noncontiguous slots."""
    if not traj.active:
        raise ValueError('trajectory is inactive')
    slots = [k for k in traj.tbf_ids if bool(storage.alive[traj.id, k])]
    if not slots:
        raise ValueError('trajectory has no active Gaussian TBFs')
    return slots


def read_packets(storage, traj):
    """Return slots, (q,p,alpha,phase), labels, coefficients, and inverse mass.

    Packet arrays are transposed to (ndof,npackets) for integral kernels.
    Advanced indexing makes temporary slices, preserving Torch graphs. A common
    positive mass per DOF is required for a coherent wavefunction.
    """
    slots = active_slots(storage, traj)
    if not storage.gaussian_vars_status:
        raise ValueError('call allocate_gaussian_vars before Gaussian dynamics')
    packet = tuple(getattr(storage, key)[traj.id, slots].T for key in
                   ('q', 'p', 'gaussian_alpha', 'gaussian_phase'))
    packet = validate_basis(*packet)
    labels = validate_states([int(scalar(v)) for v in storage.act_states_dia[traj.id, slots]],
                             len(slots), storage.nstates)
    amplitudes = storage.ampl_dia[traj.id, slots]
    rows = indices(range(len(slots)), packet[0])
    st = indices(labels, packet[0])
    c = amplitudes[rows, st]
    mask = indices(range(storage.nstates), packet[0])[None, :] == st[:, None]
    if bool((abs(amplitudes*(~mask)) > 0).any()):
        raise ValueError('state-labeled Gaussian TBFs require one coefficient per slot')
    mass = storage.iM[traj.id, slots]
    xp = xp_of(mass)
    if not bool(xp.isfinite(c).all()) or not bool(xp.isfinite(mass).all()):
        raise ValueError('coefficients and inverse masses must be finite')
    if not bool((mass > 0).all()) or not bool((mass == mass[0]).all()):
        raise ValueError('Gaussian TBFs in one trajectory require common positive inverse masses')
    return slots, packet, labels, c, mass[0]


def write_coefficients(storage, traj, slots, labels, coeff):
    """Write coherent coefficients and per-slot density matrices, without normalization.

    Clone buffers before assignment to keep previously recorded Torch graphs
    valid. Existing TensorViews continue to resolve the current storage buffers.
    Per-slot dm_dia is local bookkeeping; coherent observables need cross terms.
    """
    c = coeff.reshape(-1)
    st = indices(labels, c)
    states = indices(range(storage.nstates), c)
    amplitudes = c[:, None]*(states[None, :] == st[:, None])
    storage.ampl_dia = copy(storage.ampl_dia)
    storage.ampl_dia[traj.id, slots] = amplitudes
    storage.dm_dia = copy(storage.dm_dia)
    storage.dm_dia[traj.id, slots] = amplitudes[:, :, None]*amplitudes[:, None, :].conj()


def write_packets(storage, traj, slots, packet):
    """Commit temporary (ndof,npackets) arrays to their owning TBF slots."""
    for name, values in zip(('q', 'p', 'gaussian_alpha', 'gaussian_phase'), packet):
        dest = copy(getattr(storage, name))
        dest[traj.id, slots] = values.T
        setattr(storage, name, dest)
