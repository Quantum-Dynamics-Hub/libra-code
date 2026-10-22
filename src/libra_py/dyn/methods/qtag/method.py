"""QTAG-specific splitting and coherent observables for DynamicsEngine.

No run loop, clock, saver or independent numerical state is owned here.
Temporary pair matrices are reused only within an uninterrupted engine.run().
"""
from ...backends.arrays import scalar
from ...gwp.storage import read_packets, write_packets, write_coefficients
from ...hamiltonians.gaussian import gaussian_hamiltonian
from ...observables.gaussian import populations, gaussian_snapshot
from ...propagation.gaussian import propagate_basis
from ...propagation.nonorthogonal import propagate_coefficients


class QTAGMethod:
    """Implementation of the common prepare/advance/finish/snapshots hooks.

    advance performs electronic half-step and basis motion for every trajectory.
    The engine then advances its clock once before finish evaluates models at
    the new time and performs the final electronic half-step. Pair matrices are
    indexed by trajectory id; prepare discards them after external edits.
    """
    def __init__(self, storage, model, params):
        self.storage, self.model, self.params = storage, model, params
        self.method = params.get('qtag_pot_approx_method', 0)
        self.rcond = params.get('overlap_rcond', 1e-12)
        mp = getattr(model, 'params', {})
        self.exact_coupling = params.get('exact_coupling')
        keys = ('ex_cpl_A', 'ex_cpl_B', 'ex_cpl_C')
        if self.exact_coupling is None and all(k in mp for k in keys):
            self.exact_coupling = tuple(mp[k] for k in keys)
        self.matrices = {}

    def _assemble(self, traj):
        return gaussian_hamiltonian(self.storage, traj, self.model, self.method, self.exact_coupling)

    def prepare(self, trajectories):
        """Rebuild numerical caches from current storage and topology."""
        self.matrices = {t.id: self._assemble(t) for t in trajectories}

    def snapshots(self, trajectories, time):
        """Return coherent snapshots; independent-path estimators do not apply."""
        records = {}
        for traj in trajectories:
            result = gaussian_snapshot(self.storage, traj, *self.matrices[traj.id], time,
                                       self.params.get('target_states'))
            traj.norm = scalar(result['norm'])
            records[traj.id] = result
        return records

    def advance(self, trajectories, dt):
        """First electronic half-step and moving-basis projection."""
        work = []
        for traj in trajectories:
            slots, packet, labels, c, mass = read_packets(self.storage, traj)
            S, H = self.matrices[traj.id]
            c = propagate_coefficients(dt/2, c, H, S, self.rcond)
            states = self.params.get('states', sorted(set(labels)))
            pars = dict(self.params, dt=dt, iM=mass, states=states)
            pops = populations(labels, c, S, states)
            moved = propagate_basis(*packet, labels, c, pars, pops)
            work.append((traj, slots, labels, moved))
        for traj, slots, labels, moved in work:
            write_packets(self.storage, traj, slots, moved[:4])
            write_coefficients(self.storage, traj, slots, labels, moved[5])

    def finish(self, trajectories, dt):
        """Evaluate at the new engine time and complete electronic propagation."""
        for traj in trajectories:
            slots, packet, labels, c, mass = read_packets(self.storage, traj)
            S, H = self._assemble(traj)
            c = propagate_coefficients(dt/2, c, H, S, self.rcond)
            write_coefficients(self.storage, traj, slots, labels, c)
            self.matrices[traj.id] = (S, H)
