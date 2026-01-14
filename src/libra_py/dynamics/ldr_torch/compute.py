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
    def __init__(self, params):
        self.prefix = params.get("prefix", "ldr-solution")
        self.device = params.get("device", torch.device("cuda" if torch.cuda.is_available() else "cpu"))
        self.hbar = 1.0
        self.hamiltonian_scheme = "symmetrized"
        self.q0 = torch.tensor(params.get("q0", [0.0]), dtype=torch.float64, device=self.device)
        self.p0 = torch.tensor(params.get("p0", [0.0]), dtype=torch.float64, device=self.device)
        self.k = torch.tensor(params.get("k", [0.001]), dtype=torch.float64, device=self.device)
        self.mass = torch.tensor(params.get("mass", [2000.0]), dtype=torch.float64, device=self.device)
        self.alpha = torch.tensor(params.get("alpha", [18.0]), dtype=torch.float64, device=self.device)
        self.qgrid = torch.tensor(params.get("qgrid", [[-10 + i * 0.1] for i in range(int((10 - (-10)) / 0.1) + 1)] ), dtype=torch.float64, device=self.device) #(N, D)
        self.ngrids = len(self.qgrid) # N
        self.ndof = self.qgrid.shape[1] 
        self.nstates = params.get("nstates", 2)
        self.istate = params.get("istate", 0)
        self.elec_ampl = params.get("elec_ampl", torch.tensor([1.0+0.j]*self.ngrids, dtype=torch.cdouble))

        self.save_every_n_steps = params.get("save_every_n_steps", 1)
        self.properties_to_save = params.get("properties_to_save", ["time", "population_right"])
        self.dt = params.get("dt", 0.01)
        self.nsteps = params.get("nsteps", 500)
        self.ndim = self.nstates * self.ngrids

        self.E = params.get("E", torch.zeros(self.nstates, self.ngrids, device=self.device) )

        s_elec_default = torch.zeros(self.ndim, self.ndim, dtype=torch.cdouble, device=self.device)
        for i in range(self.nstates):
            start, end = i * self.ngrids, (i + 1) * self.ngrids
            s_elec_default[start:end, start:end] = torch.eye(self.ngrids, device=self.device)            
        self.s_elec = params.get("s_elec", s_elec_default )   

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
        exponent = -0.5 * torch.sum(self.alpha * delta**2, dim=2)  # (N, N)
        self.s_nucl = torch.exp(exponent)

    def chi_kinetic(self):
        """
        Compute nuclear kinetic energy matrix t_nucl[i,j] = <g(x; qgrid[i]) | T | g(x; qgrid[j])>,
        with T = \sum_{\nu} -0.5* m_ν^{-1} \partial^{2}/\partial x_{\nu}^2.
        """
        delta = self.qgrid[:, None, :] - self.qgrid[None, :, :]               # (N, N, D)
        tau = self.alpha / (2.0 * self.mass) * (1.0 - self.alpha * delta**2)  # (N, N, D)
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
    
        s_nucl_4d = self.s_nucl.unsqueeze(0).unsqueeze(2) # (1, n, 1, m)
    
        S_4d = s_elec_4d * s_nucl_4d
    
        # Reshape back to (ndim, ndim) with compound indices
        self.S = S_4d.permute(0, 1, 2, 3).reshape(ndim, ndim)

    def build_compound_hamiltonian(self):
        """
        Build the compound nuclear-electronic Hamiltonian self.H (ndim, ndim) using different schemes.
        """
        N, s, ndim = self.ngrids, self.nstates, self.ndim
        scheme = self.hamiltonian_scheme
        s_elec_4d = self.s_elec.view(s, N, s, N)      # (s, N, s, N)
        T_4d = self.t_nucl.unsqueeze(0).unsqueeze(2)  # (1, N, 1, N)
        S_4d = self.s_nucl.unsqueeze(0).unsqueeze(2)  # (1, N, 1, N)
    
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
            delta_ij = torch.eye(s, device=self.device).unsqueeze(1).unsqueeze(3)  # (s, 1, s, 1)
            delta_nm = torch.eye(N, device=self.device).unsqueeze(0).unsqueeze(2)  # (1, N, 1, N)
            delta_4d = delta_ij * delta_nm
            
            E_j_4d = self.E[None, None, :, :] # (1, 1, s, N)
            bracket_4d = T_4d + E_j_4d * S_4d * delta_4d

        else:
            raise ValueError(f"Unknown Hamiltonian scheme: {scheme}")
    
        H_4d = s_elec_4d * bracket_4d
        self.H = H_4d.reshape(ndim, ndim)

    def compute_propagator(self):
        """
        Compute the exponential propagator matrix U = exp(-i H dt) in the non-orthogonal basis
        using the Lowdin orthonormalization.
  
        """
        S = self.S
        H = self.H
        dt = self.dt
    
        evals_S, evecs_S = torch.linalg.eigh(S)
    
        self.S_half = (evecs_S @ torch.diag(evals_S.sqrt().to(dtype=torch.cdouble)) @ evecs_S.T).to(dtype=torch.cdouble)
        S_invhalf = (evecs_S @ torch.diag((1.0 / evals_S).sqrt().to(dtype=torch.cdouble)) @ evecs_S.T).to(dtype=torch.cdouble)
    
        H_ortho = S_invhalf @ H @ S_invhalf
    
        evals_H, evecs_H = torch.linalg.eigh(H_ortho)
    
        exp_diag = torch.diag(torch.exp(-1j * evals_H * dt))
        U_ortho = evecs_H @ exp_diag @ evecs_H.conj().T
    
        self.U = S_invhalf @ U_ortho @ self.S_half


    def initialize_C(self):
        """
        Initialize coefficient vector self.C0 at t=0, assuming:
        - electronic state self.istate
        - nuclear wavefunction is a Gaussian centered at self.q0 and self.p0, i.e.
        chi0(q;q0,p0) = exp( - alpha0 * (q - self.q0)**2 + i * self.p0 * (q - self.q0) )
        alpha0 = 0.5/s_q **2; s_q = (1/(self.k * self.mass)) **0.25
        
        For the Gaussian overlap calculation, consult with the formula in Begušić, T.; Vaníček, J. J. Chem. Phys. 2020, 153 (18), 184110.

        Sets:
            self.C0 : complex-valued coefficient vector of shape (ndim)
        """
        N, ist = self.ngrids, self.istate

        q0     = self.q0.to(torch.cdouble)
        p0     = self.p0.to(torch.cdouble)
        qgrid  = self.qgrid.to(torch.cdouble)
        alpha  = self.alpha.to(torch.cdouble)

        s_q = (1.0 / (self.k*self.mass) ) ** 0.25
        alpha0 = 1 / ( 2 * s_q **2 )
        alpha0 = alpha0.to(torch.cdouble)

        # Width matrix
        Ag, A = torch.diag(2.j * self.alpha), torch.diag(2.j * alpha0)
        delta_A = A - Ag.conj()
        delta_A_inv = torch.torch.linalg.inv(delta_A)

        # Compute Gaussian nuclear wavefunction at each grid point
        for n in range(N):
            index = ist * N + n
    
            xi0, xig = p0 - torch.matmul(A, q0), -torch.matmul(Ag, qgrid[n])
            delta_xi = xi0 - xig.conj()
            delta_eta = -0.5 * torch.dot(xi0 + p0, q0) + 0.5 * torch.dot(xig, qgrid[n]).conj()
            exponent = -1.j * 0.5 * torch.dot(delta_xi, torch.matmul(delta_A_inv, delta_xi)) + 1.j * delta_eta
    
            self.C0[index] = self.elec_ampl[n] * torch.exp(exponent)
    
        # Normalize
        overlap = torch.matmul(self.S, self.C0)
        norm = torch.sqrt(torch.vdot(self.C0, overlap))
    
        self.C0 /= norm
    
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
            q_med = 0.5 * (self.qgrid[:, None, idof] + self.qgrid[None,:,idof])
            q_nucl = self.s_nucl * q_med 
            Q_4d = q_nucl[None, :, None, :]
            Q_4d_compound = s_elec_4d * Q_4d
            Q_compound = Q_4d_compound.permute(0, 1, 2, 3).reshape(ndim, ndim)

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

