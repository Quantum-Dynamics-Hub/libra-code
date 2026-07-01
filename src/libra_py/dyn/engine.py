# electronic.py

from dataclasses import dataclass
import numpy as np



# ============================================================
# Electronic state container
# ============================================================

@dataclass
class ElectronicState:
    coeff: np.ndarray          # C matrix
    rho: np.ndarray           # density matrix
    active_states: np.ndarray
    energies: np.ndarray
    nac: np.ndarray


# ============================================================
# TSHEngine
# ============================================================

class TSHEngine:
    """
    Full trajectory surface hopping engine.
    """

    def __init__(self, ham_engine, params, rng=None):

        self.ham = ham_engine
        self.params = params
        self.rng = rng or np.random.default_rng()

        # internal state
        self.state = None

    # --------------------------------------------------------
    # Initialization
    # --------------------------------------------------------

    def initialize(self, q, p, electronic_state: ElectronicState):
        self.q = q
        self.p = p
        self.state = electronic_state

        # initial Hamiltonian evaluation
        self.ham.evaluate(self.q, None, mode=0)
        self.state = self._sync_from_hamiltonian()

    # --------------------------------------------------------
    # Main propagation step
    # --------------------------------------------------------

    def step(self, dt):

        # ====================================================
        # 1. Nuclear propagation (half step)
        # ====================================================
        self._propagate_nuclei_half(dt)

        # ====================================================
        # 2. Hamiltonian update (q changed)
        # ====================================================
        ham_state = self.ham.evaluate(self.q, None, mode=0)
        self._update_from_ham(ham_state)

        # ====================================================
        # 3. Electronic propagation (coefficients)
        # ====================================================
        self._propagate_electronic(dt)

        # ====================================================
        # 4. Momentum update (forces depend on electrons)
        # ====================================================
        self._compute_forces()
        self._propagate_momenta_half(dt)

        # ====================================================
        # 5. Hamiltonian update (vibronic correction)
        # ====================================================
        ham_state = self.ham.evaluate(self.q, self.p, mode=1)
        self._update_from_ham(ham_state)

        # ====================================================
        # 6. Hopping decision
        # ====================================================
        self._surface_hopping(dt)

        # ====================================================
        # 7. Final consistency update
        # ====================================================
        ham_state = self.ham.evaluate(self.q, self.p, mode=1)
        self._update_from_ham(ham_state)

    # --------------------------------------------------------
    # Nuclear propagation
    # --------------------------------------------------------

    def _propagate_nuclei_half(self, dt):
        self.q += 0.5 * dt * self.p

    def _propagate_momenta_half(self, dt):
        self.p += 0.5 * dt * self.state_forces()

    # --------------------------------------------------------
    # Electronic propagation
    # --------------------------------------------------------

    def _propagate_electronic(self, dt):
        C = self.state.coeff

        H = self.state.energies
        NAC = self.state.nac

        # placeholder TDSE propagation
        dC = -1j * H @ C * dt + NAC @ C * dt
        self.state.coeff += dC

    # --------------------------------------------------------
    # Forces
    # --------------------------------------------------------

    def _compute_forces(self):
        # placeholder
        self.forces = -np.real(self.state.nac)

    def state_forces(self):
        return self.forces

    # --------------------------------------------------------
    # Hopping
    # --------------------------------------------------------

    def _surface_hopping(self, dt):
        """
        Minimal FSSH-like logic placeholder.
        """

        pops = np.abs(self.state.coeff) ** 2
        active = self.state.active_states

        g = self._compute_hopping_probabilities()

        new_states = active.copy()

        for traj in range(len(active)):
            r = self.rng.random()
            cum = 0.0

            for j in range(len(pops)):
                if j == active[traj]:
                    continue
                cum += g[traj, j]

                if r < cum:
                    new_states[traj] = j
                    break

        self.state.active_states = new_states

    # --------------------------------------------------------
    # hopping probabilities (placeholder)
    # --------------------------------------------------------

    def _compute_hopping_probabilities(self):
        nst = len(self.state.coeff)
        return np.abs(np.random.randn(nst, nst))

    # --------------------------------------------------------
    # sync helpers
    # --------------------------------------------------------

    def _sync_from_hamiltonian(self):
        return self.ham.backend.get_state()

    def _update_from_ham(self, ham_state):
        self.state.energies = ham_state.energies
        self.state.nac = ham_state.nac
