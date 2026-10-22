"""Coherent Gaussian observables, including inter-packet overlap terms.

The independent-TBF population/energy estimators in observables are unsuitable
for a coherent expansion. Here a trajectory's norm is C†SC, its energy C†HC,
and populations include all overlaps between packets with the same state label.
"""
import numpy as np
from ..backends.arrays import asarray, copy, indices, stack, xp_of
from ..gwp.storage import read_packets
from ..gwp.packets import wavefunction


def populations(labels, coeff, overlap, target_states):
    """Return unnormalized state populations in the requested label order."""
    c = coeff.reshape(-1)
    labels = indices(labels, c)
    values = []
    for state in target_states:
        selected = c*(labels == state)
        values.append((selected.conj() @ overlap @ selected).real)
    return stack(values)


def gaussian_snapshot(storage, traj, overlap, hamiltonian, time, target_states=None):
    """Copy one coherent snapshot, using storage's (npackets,ndof) layout.

    Arrays stay on their original backend/device. Call savers.save_step with
    this mapping to detach only at the explicit disk boundary. Slots and labels
    are saved alongside values, so noncontiguous/pruned topologies are explicit.
    ``coeffs`` are scalar packet coefficients; ``Etot`` is not divided by norm.
    """
    slots, (q, p, a, s), labels, c, mass = read_packets(storage, traj)
    targets = list(range(storage.nstates)) if target_states is None else list(target_states)
    if not targets or any(int(x) != x or not 0 <= x < storage.nstates for x in targets):
        raise ValueError('target_states must be valid electronic labels')
    values = dict(time=asarray(time, q), step=indices(storage.timestep, q),
                  q=q.T, p=p.T, gaussian_alpha=a.T, gaussian_phase=s.T,
                  coeffs=c, states=indices(labels, q), slot_ids=indices(slots, q),
                  target_states=indices(targets, q),
                  norm=(c.conj() @ overlap @ c).real,
                  Etot=(c.conj() @ hamiltonian @ c).real,
                  pops=populations(labels, c, overlap, targets))
    return {key: copy(value) for key, value in values.items()}


def reconstruct_density(history, xmin, xmax, npoints, snaps=None, which_states=None):
    """Reconstruct stored Gaussian snapshots on a Cartesian spatial grid.

    ``history`` stacks gaussian_snapshot mappings along a leading time axis.
    ``snaps`` indexes those records (index zero is the initial state); each
    record carries its own slot/state labels. Output density and psi have shape
    (nsnap,nstates,*npoints), coordinates have (ndof,npoints_total). Phases are
    retained. An arbitrary ordered subset of electronic states can be selected.
    """
    q = history['q']
    d = q.shape[-1]
    if len(xmin) != d or len(xmax) != d or len(npoints) != d:
        raise ValueError('grid needs one bound/count per DOF')
    if any(int(n) != n or n < 2 for n in npoints) or any(b <= a for a, b in zip(xmin, xmax)):
        raise ValueError('grid counts must be >=2 and upper bounds exceed lower bounds')
    axes = [asarray(np.linspace(a, b, int(n)), q) for a, b, n in zip(xmin, xmax, npoints)]
    points = stack(xp_of(q).meshgrid(*axes, indexing='ij')).reshape(d, -1)
    snaps = [0, len(q)-1] if snaps is None else list(snaps)
    states = [int(x) for x in history['target_states'][0]] if which_states is None else list(which_states)
    values = []
    for k in snaps:
        if not 0 <= k < len(q):
            raise ValueError('snapshot index is outside history')
        labels, c = history['states'][k], history['coeffs'][k]
        coefficients = stack([c*(labels == state) for state in states])
        values.append(wavefunction(points, q[k].T, history['p'][k].T,
                                   history['gaussian_alpha'][k].T,
                                   history['gaussian_phase'][k].T, coefficients)
                      .reshape(len(states), *npoints))
    psi = stack(values)
    return dict(axes=axes, points=points, psi=psi, density=abs(psi)**2, snaps=snaps,
                states=states, times=history['time'][indices(snaps, q)])
