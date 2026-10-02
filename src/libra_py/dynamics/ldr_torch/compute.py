# *********************************************************************************
# * Copyright (C) 2025 Daeho Han and Alexey V. Akimov
# *
# * This file is distributed under the terms of the GNU General Public License
# * as published by the Free Software Foundation, either version 3 of
# * the License, or (at your option) any later version.
# * See the file LICENSE in the root directory of this distribution
# * or <http://www.gnu.org/licenses/>.
# ***********************************************************************************
"""
.. module:: compute
   :platform: Unix, Windows
   :synopsis: This module implements functions for doing local diabatic representation (LDR) dynamics with PyTorch
       List of functions:
           * sech # temporary here
           * Martens_model # temporary here
           * gaussian_wavepacket
       List of classes:
           * ldr_solver

.. moduleauthor:: Daeho Han and Alexey V. Akimov

"""

__author__ = "Daeho Han, Alexey V. Akimov"
__copyright__ = "Copyright 2025 Alexey V. Akimov"
__credits__ = ["Daeho Han", "Alexey V. Akimov"]
__license__ = "GNU-3"
__version__ = "1.0"
__maintainer__ = "Alexey V. Akimov"
__email__ = "alexvakimov@gmail.com"
__url__ = "https://github.com/Quantum-Dynamics-Hub/libra-code"



import torch


class ldr_solver:
    """Dynamics in a fixed, nonorthogonal Gaussian/electronic basis.

    ``elec_ampl[j, n]`` is the overlap of electronic basis state ``(j, n)``
    with the desired initial electronic state, not an expansion coefficient.
    Supply all states, shape ``(nstates, ngrids)``, for a fixed electronic
    reference state. A legacy vector of shape ``(ngrids,)`` supplies overlaps
    only in ``istate``; omitted overlaps are zero. The default assumes
    coordinate-independent orthonormal electronic states.

    ``alpha`` contains positive Gaussian exponents: a scalar or a vector
    of length ``ndof`` gives uniform widths; shape ``(ngrids, ndof)`` gives
    a separate width for each center and dimension. A one-element vector
    broadcasts over dimensions. The nuclear basis functions are normalized.

    Energies have shape ``(nstates, ngrids)`` and electronic overlaps have
    shape ``(nstates * ngrids, nstates * ngrids)``, in state-major order.
    Inputs are converted to float64/complex128 on ``device``.
    """

    def __init__(self, params):
        self.prefix = params.get("prefix", "ldr-solution")
        self.device = torch.device(params.get("device", "cuda" if torch.cuda.is_available() else "cpu"))
        self.hbar = 1.0
        self.hamiltonian_scheme = "symmetrized"
        self.q0 = torch.as_tensor(params.get("q0", [0.0]), dtype=torch.float64, device=self.device)
        self.p0 = torch.as_tensor(params.get("p0", [0.0]), dtype=torch.float64, device=self.device)
        self.k = torch.as_tensor(params.get("k", [0.001]), dtype=torch.float64, device=self.device)
        self.mass = torch.as_tensor(params.get("mass", [2000.0]), dtype=torch.float64, device=self.device)
        self.alpha = torch.as_tensor(params.get("alpha", [18.0]), dtype=torch.float64, device=self.device)
        self.qgrid = torch.as_tensor(params.get("qgrid", [[-10 + i * 0.1] for i in range(int((10 - (-10)) / 0.1) + 1)] ), dtype=torch.float64, device=self.device) #(N, D)
        self.ngrids = len(self.qgrid) # N
        self.ndof = self.qgrid.shape[1] 
        if not torch.all(torch.isfinite(self.alpha) & (self.alpha > 0)):
            raise ValueError("Gaussian widths must be positive and finite")
        if self.alpha.ndim == 0 or (
                self.alpha.ndim == 1 and self.alpha.numel() in (1, self.ndof)):
            self._alpha_grid = self.alpha.expand(self.ngrids, self.ndof)
        elif self.alpha.shape == (self.ngrids, self.ndof):
            self._alpha_grid = self.alpha
        else:
            raise ValueError("alpha must be scalar, have shape (ndof,), "
                             "or have shape (ngrids, ndof)")
        self.nstates = params.get("nstates", 2)
        self.istate = params.get("istate", 0)
        if not 0 <= self.istate < self.nstates:
            raise ValueError("istate must index one of nstates electronic states")
        amplitudes = torch.as_tensor(
            params.get("elec_ampl", torch.ones(self.ngrids)),
            dtype=torch.cdouble, device=self.device,
        )
        if amplitudes.shape == (self.ngrids,):
            self.elec_ampl = torch.zeros(
                self.nstates, self.ngrids, dtype=torch.cdouble, device=self.device
            )
            self.elec_ampl[self.istate] = amplitudes
        elif amplitudes.shape == (self.nstates, self.ngrids):
            self.elec_ampl = amplitudes
        else:
            raise ValueError("elec_ampl must have shape (ngrids,) or (nstates, ngrids)")

        self.save_every_n_steps = params.get("save_every_n_steps", 1)
        self.properties_to_save = params.get("properties_to_save", ["time", "population_right"])
        self.dt = params.get("dt", 0.01)
        self.nsteps = params.get("nsteps", 500)
        self.ndim = self.nstates * self.ngrids

        energies = torch.as_tensor(
            params.get("E", torch.zeros(self.nstates, self.ngrids)),
            dtype=torch.cdouble, device=self.device,
        )
        if energies.shape != (self.nstates, self.ngrids):
            raise ValueError("E must have shape (nstates, ngrids)")
        if energies.is_complex():
            if torch.any(energies.imag != 0):
                raise ValueError("E must contain real electronic energies")
            energies = energies.real
        self.E = energies.to(dtype=torch.float64)

        if "s_elec" in params:
            self.s_elec = torch.as_tensor(
                params["s_elec"], dtype=torch.cdouble, device=self.device
            )
        else:
            # The same electronic state overlaps with itself at EVERY center.
            self.s_elec = torch.kron(
                torch.eye(self.nstates, dtype=torch.cdouble, device=self.device),
                torch.ones(self.ngrids, self.ngrids, dtype=torch.cdouble, device=self.device),
            )
        if self.s_elec.shape != (self.ndim, self.ndim):
            raise ValueError("s_elec must have shape (nstates * ngrids, nstates * ngrids)")
        self.s_elec = self.s_elec.contiguous()

        # Computed with LDR methods
        self.C0 = torch.zeros(self.ndim, dtype=torch.cdouble, device=self.device)
        self.C_curr = torch.zeros(self.ndim, dtype=torch.cdouble, device=self.device)

        self.s_nucl = torch.eye(self.ngrids, dtype=torch.cdouble, device=self.device)
        self.t_nucl = torch.zeros(self.ngrids, self.ngrids, dtype=torch.cdouble, device=self.device)

        self.S, self.H = torch.zeros(self.ndim, self.ndim, dtype=torch.cdouble, device=self.device), torch.zeros(self.ndim, self.ndim, dtype=torch.cdouble, device=self.device)
        self.U = torch.zeros(self.ndim, self.ndim, dtype=torch.cdouble, device=self.device)
        self.S_half = torch.zeros(self.ndim, self.ndim, dtype=torch.cdouble, device=self.device)
        
        self.time = []
        self.kinetic_energy = []
        self.potential_energy = []
        self.total_energy = []
        self.average_pos = []
        self.population_right = []
        self.denmat = []
        self.norm = []
        self.C_save = []

    def chi_overlap(self):
        """
        Compute nuclear overlap matrix s_nucl[i, j] for the mesh qmesh 
        from the Gaussian basis, g(x; q) = \exp(-\alpha * (x-q)**2).
        """
        delta = self.qgrid[:, None, :] - self.qgrid[None, :, :]    # (N, N, D)
        if self.alpha.ndim < 2:
            exponent = -0.5 * torch.sum(self.alpha * delta**2, dim=2)
            self.s_nucl = torch.exp(exponent)
        else:
            alpha_i = self._alpha_grid[:, None, :]
            alpha_j = self._alpha_grid[None, :, :]
            alpha_sum = alpha_i + alpha_j
            beta = alpha_i * alpha_j / alpha_sum
            prefactor = torch.sqrt(2 * torch.sqrt(alpha_i * alpha_j) / alpha_sum)
            self.s_nucl = prefactor.prod(dim=2) * torch.exp(
                -torch.sum(beta * delta**2, dim=2))

    def chi_kinetic(self):
        """
        Compute nuclear kinetic energy matrix t_nucl[i,j] = <g(x; qgrid[i]) | T | g(x; qgrid[j])>,
        with T = \sum_{\nu} -0.5* m_ν^{-1} \partial^{2}/\partial x_{\nu}^2.
        """
        delta = self.qgrid[:, None, :] - self.qgrid[None, :, :]               # (N, N, D)
        if self.alpha.ndim < 2:
            tau = self.alpha / (2.0 * self.mass) * (1.0 - self.alpha * delta**2)
        else:
            alpha_i = self._alpha_grid[:, None, :]
            alpha_j = self._alpha_grid[None, :, :]
            beta = alpha_i * alpha_j / (alpha_i + alpha_j)
            tau = beta / self.mass * (1.0 - 2.0 * beta * delta**2)
        tau_sum = torch.sum(tau, dim=2)                                       # (N, N)
    
        self.t_nucl = self.s_nucl * tau_sum                                   # (N, N)

    def build_compound_overlap(self):
        """
        Build the compound nuclear-electronic overlap matrix self.S (ndim, ndim)
        """
        N, s, ndim = self.ngrids, self.nstates, self.ndim
    
        # Reshape s_elec[a, b] -> (i, n, j, m) with:
        #   a = i * N + n
        #   b = j * N + m
        s_elec_4d = self.s_elec.view(s, N, s, N) # (i, n, j, m)
    
        s_nucl_4d = self.s_nucl[None, :, None, :] # (1, n, 1, m)
    
        S_4d = s_elec_4d * s_nucl_4d
    
        # Reshape back to (ndim, ndim) with compound indices
        self.S = S_4d.reshape(ndim, ndim)

    def build_compound_hamiltonian(self):
        """
        Build the compound nuclear-electronic Hamiltonian self.H (ndim, ndim) using different schemes.
        """
        N, s, ndim = self.ngrids, self.nstates, self.ndim
        scheme = self.hamiltonian_scheme
        s_elec_4d = self.s_elec.view(s, N, s, N)      # (s, N, s, N)
        T_4d = self.t_nucl[None, :, None, :]          # (1, N, 1, N)
        S_4d = self.s_nucl[None, :, None, :]          # (1, N, 1, N)
    
        if scheme == 'as_is': # For showing the original non-Hermitian form, not intended to use
            E_j_4d = self.E[None, None, :, :]   # (1, 1, s, N)
            bracket_4d = T_4d + E_j_4d * S_4d
        elif scheme == 'symmetrized':
            E_i_4d = self.E[:, :, None, None]   # (s, N, 1, 1)
            E_j_4d = self.E[None, None, :, :]   # (1, 1, s, N)
            E_avg_4d = 0.5 * (E_i_4d + E_j_4d)  # (s, N, s, N)
            bracket_4d = T_4d + E_avg_4d * S_4d
        elif scheme == 'diagonal':
            # Build Kronecker deltas for electronic and nuclear indices
            delta_ij = torch.eye(s, device=self.device)[:, None, :, None]  # (s, 1, s, 1)
            delta_nm = torch.eye(N, device=self.device)[None, :, None, :]  # (1, N, 1, N)
            delta_4d = delta_ij * delta_nm
            
            E_j_4d = self.E[None, None, :, :] # (1, 1, s, N)
            bracket_4d = T_4d + E_j_4d * S_4d * delta_4d

        else:
            raise ValueError(f"Unknown Hamiltonian scheme: {scheme}")
    
        H_4d = s_elec_4d * bracket_4d
        self.H = H_4d.reshape(ndim, ndim)

    def compute_propagator(self):
        """
        Compute U = S^-1/2 exp(-i S^-1/2 H S^-1/2 dt / hbar) S^1/2.

        The overlap must be positive definite; linearly dependent basis
        functions must be removed before propagation.
  
        """
        S = self.S
        H = self.H
        dt = self.dt
    
        evals_S, evecs_S = torch.linalg.eigh(S)
        if not torch.all(torch.isfinite(evals_S) & (evals_S > 0)):
            raise ValueError("The compound overlap S must be positive definite")
        evecs_S = evecs_S.to(dtype=torch.cdouble)
        # eigh returns eigenvectors as COLUMNS of Q: S Q = Q diag(evals_S).
        # Thus S^p = Q diag(evals_S^p) Q^dagger, not Q^dagger diag(...) Q.
        # Broadcasting the eigenvalue vector below scales Q's columns.
        self.S_half = (evecs_S * evals_S.sqrt()) @ evecs_S.conj().T
        S_invhalf = (evecs_S * evals_S.rsqrt()) @ evecs_S.conj().T
    
        H_ortho = S_invhalf @ H @ S_invhalf
    
        evals_H, evecs_H = torch.linalg.eigh(H_ortho)
    
        exp_diag = torch.diag(torch.exp(-1j * evals_H * dt / self.hbar))
        U_ortho = evecs_H @ exp_diag @ evecs_H.conj().T
    
        self.U = S_invhalf @ U_ortho @ self.S_half


    def initialize_C(self):
        """
        Project a normalized Gaussian onto the compound basis and normalize.

        The target nuclear wavefunction is proportional to
        exp[-sum(alpha0 * (q-q0)^2) + i p0.(q-q0)/hbar], where
        alpha0 = sqrt(k * mass)/2. Electronic overlaps are supplied through
        ``elec_ampl`` (see the class docstring). Build S before calling this
        method. Solve S C0 = b with b_a = <psi_a|Psi0>; the overlaps themselves
        are NOT expansion coefficients in a nonorthogonal basis.
        """
        alpha0 = 0.5 * torch.sqrt(self.k * self.mass)
        if torch.any(self.alpha <= 0) or torch.any(alpha0 <= 0):
            raise ValueError("Gaussian widths must be positive")
        delta = self.qgrid - self.q0
        alpha = self._alpha_grid
        width = alpha + alpha0
        momentum = self.p0 / self.hbar
        # Analytic overlap of normalized Gaussians, in coordinates relative
        # to q0 to avoid cancellation between large absolute positions.
        prefactor = torch.prod((4 * alpha * alpha0 / width**2)**0.25, dim=1)
        exponent = torch.sum(
            -alpha * alpha0 / width * delta**2
            -momentum**2 / (4 * width)
            +1j * alpha / width * momentum * delta,
            dim=1,
        )
        nuclear_overlap = prefactor * torch.exp(exponent)
        b = (self.elec_ampl * nuclear_overlap[None, :]).reshape(self.ndim)
        coefficients = torch.linalg.solve(self.S, b)
        norm_squared = torch.vdot(coefficients, self.S @ coefficients).real
        if not torch.isfinite(norm_squared) or norm_squared <= 0:
            raise ValueError("The initial state must have a finite, nonzero projected norm")
        self.C0 = coefficients / norm_squared.sqrt()
    
    def propagate(self):
        """
        Propagate coefficient.
        """
        # Initialize first step with normalized initial wavefunction
        self.C_curr = self.C0.clone()

        print(F"step = 0")
        self.save_results(0)
        
        for step in range(1, self.nsteps):
            C_vec = self.C_curr.clone()
            self.C_curr = self.U @ C_vec

            if step % self.save_every_n_steps == 0:
                print(F"step = {step}")
                self.save_results(step)

    def save_results(self, step):
        if "time" in self.properties_to_save:
            self.time.append(step*self.dt)
        if "norm" in self.properties_to_save:
            overlap = torch.matmul(self.S, self.C_curr)
            self.norm.append(torch.sqrt(torch.vdot(self.C_curr, overlap)))
        if "population_right" in self.properties_to_save:
            self.population_right.append(self.compute_populations())
        if "denmat" in self.properties_to_save:
            self.denmat.append(self.compute_denmat())
        if "kinetic_energy" in self.properties_to_save:
            self.kinetic_energy.append(self.compute_kinetic_energy())
        if "potential_energy" in self.properties_to_save:
            self.potential_energy.append(self.compute_potential_energy())
        if "total_energy" in self.properties_to_save:
            self.total_energy.append(self.compute_total_energy())
        if "average_pos" in self.properties_to_save:
            self.average_pos.append(self.compute_average_pos())
        if "C_save" in self.properties_to_save:
            self.C_save.append(self.C_curr)
        
    def compute_populations(self):
        """
        Compute electronic state population for a single step.
        """
        N, s = self.ngrids, self.nstates
        C_vec = self.C_curr
        
        # Compute SC once: shape (ndim,)
        SC = self.S @ C_vec
    
        C_blocks = C_vec.view(s, N)
        SC_blocks = SC.view(s, N)
    
        # Compute P[i] = sum_j <C_j|S_{ji}|C_i> = Re[ sum_N (C_j*) * SC_j ]
        P = torch.sum(C_blocks.conj() * SC_blocks, dim=1).real
    
        return P

    def compute_denmat(self):
        """
        Compute electronic density matrix for a single step using the orthogonalization.
        """
        N, s = self.ngrids, self.nstates
        C_vec = self.C_curr
    
        # Orthogonalize coefficients: C_ortho = S^{1/2} C
        C_ortho = self.S_half @ C_vec
      
        C_blocks = C_ortho.view(s, N)
    
        rho = C_blocks @ C_blocks.conj().T # (s, s)
    
        return rho

    def compute_kinetic_energy(self):
        """
        Compute nuclear kinetic energy as C^+ T C / C^+ S C for a single step.
        """
        N, s, ndim = self.ngrids, self.nstates, self.ndim
    
        # Rebuild compound kinetic matrix: T_4d * s_elec_4d
        s_elec_4d = self.s_elec.view(s, N, s, N)
        T_4d = self.t_nucl[None, :, None, :]
        T_compound = (s_elec_4d * T_4d).reshape(ndim, ndim)
    
        C_vec = self.C_curr
    
        numer = torch.vdot(C_vec, T_compound @ C_vec).real
        denom = torch.vdot(C_vec, self.S @ C_vec).real
    
        return numer / denom
    
    
    def compute_potential_energy(self):
        """
        Compute potential energy as C^+ V C / C^+ S C for a single step.
        """
        N, s, ndim = self.ngrids, self.nstates, self.ndim
    
        s_elec_4d = self.s_elec.view(s, N, s, N)
        S_4d = self.s_nucl[None, :, None, :]
        E_j_4d = self.E[None, None, :, :]  # (1,1,j,m)
    
        V_compound = (s_elec_4d * (E_j_4d * S_4d)).reshape(ndim, ndim)
    
        C_vec = self.C_curr
    
        numer = torch.vdot(C_vec, V_compound @ C_vec).real
        denom = torch.vdot(C_vec, self.S @ C_vec).real
    
        return numer / denom
    
    
    def compute_total_energy(self):
        """
        Compute total energy as C^+ H C / C^+ S C for a single step.
        """
        C_vec = self.C_curr
    
        numer = torch.vdot(C_vec, self.H @ C_vec).real
        denom = torch.vdot(C_vec, self.S @ C_vec).real
    
        return numer / denom

    def compute_average_pos(self):
        """
        Compute average position as <q_i> = \sum_i C^+ Q C / C^+ S C for a single step.
        """
        N, s, ndim = self.ngrids, self.nstates, self.ndim
        
        C_vec = self.C_curr

        denom = torch.vdot(C_vec, self.S @ C_vec).real
        s_elec_4d = self.s_elec.view(s, N, s, N)
        
        avg_q = []
        for idof in range(self.ndof):
            alpha_i = self._alpha_grid[:, None, idof]
            alpha_j = self._alpha_grid[None, :, idof]
            q_med = (alpha_i * self.qgrid[:, None, idof]
                     + alpha_j * self.qgrid[None, :, idof]) / (alpha_i + alpha_j)
            q_nucl = self.s_nucl * q_med 
            Q_4d = q_nucl[None, :, None, :]
            Q_4d_compound = s_elec_4d * Q_4d
            Q_compound = Q_4d_compound.reshape(ndim, ndim)

            numer = torch.vdot(C_vec, Q_compound @ C_vec).real
            avg_q.append(numer / denom)

        return avg_q

    def save(self):
        torch.save( {"q0":self.q0,
                     "p0":self.p0,
                     "k":self.k,
                     "mass":self.mass,
                     "alpha":self.alpha,
                     "qgrid":self.qgrid,
                     "nstates":self.nstates,
                     "istate":self.istate,
                     "s_nucl":self.s_nucl,
                     "t_nucl":self.t_nucl,
                     "E":self.E,
                     "s_elec":self.s_elec,
                     "S":self.S,
                     "H":self.H,
                     "U":self.U,
                     "C_save":self.C_save,
                     "save_every_n_steps":self.save_every_n_steps,
                     "hamiltonian_scheme": self.hamiltonian_scheme,
                     "dt":self.dt, "nsteps":self.nsteps,
                     "time":self.time,
                     "kinetic_energy":self.kinetic_energy,
                     "potential_energy":self.potential_energy,
                     "total_energy":self.total_energy,
                     "average_pos":self.average_pos,
                     "population_right":self.population_right,
                     "denmat":self.denmat,
                     "norm":self.norm
                    }, F"{self.prefix}.pt" )

    def buildSH(self):
        self.chi_overlap()
        self.chi_kinetic()
        self.build_compound_overlap()
        self.build_compound_hamiltonian()
    
    def solve(self):
        print("Building overlap and Hamiltonian matrices")
        self.buildSH()
        print("Computing the time propagator")
        self.compute_propagator()
        print("Initializing Coefficients")
        self.initialize_C()
        print("Propagating Coefficients")
        self.propagate()
        self.save()
