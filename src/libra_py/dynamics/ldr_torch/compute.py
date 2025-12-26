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
        self.Hamiltonian_scheme = "symmetrized"
        self.q0 = torch.tensor(params.get("q0", [0.0]), dtype=torch.float64, device=self.device)
        self.p0 = torch.tensor(params.get("p0", [0.0]), dtype=torch.float64, device=self.device)
        self.k = torch.tensor(params.get("k", [0.001]), dtype=torch.float64, device=self.device)
        self.mass = torch.tensor(params.get("mass", [2000.0]), dtype=torch.float64, device=self.device)
        self.alpha = torch.tensor(params.get("alpha", [18.0]), dtype=torch.float64, device=self.device)
        self.qgrid = torch.tensor(params.get("qgrid", [[-10 + i * 0.1] for i in range(int((10 - (-10)) / 0.1) + 1)] ), dtype=torch.float64, device=self.device) #(N, D)
        self.ngrids = len(self.qgrid) # N
        self.nstates = params.get("nstates", 2)
        self.istate = params.get("istate", 0)
        self.elec_ampl = params.get("elec_ampl", torch.tensor([1.0+0.j]*self.ngrids, dtype=torch.cdouble))

        self.save_every_n_steps = params.get("save_every_n_steps", 1)
        self.properties_to_save = params.get("properties_to_save", ["time", "population_right"])
        self.dt = params.get("dt", 0.01)
        self.nsteps = params.get("nsteps", 500)
        self.ndim = self.nstates * self.ngrids

        self.E = params.get("E", torch.zeros(self.nstates, self.ngrids, device=self.device) )

        Selec_default = torch.zeros(self.ndim, self.ndim, dtype=torch.cdouble, device=self.device)
        for i in range(self.nstates):
            start, end = i * self.ngrids, (i + 1) * self.ngrids
            Selec_default[start:end, start:end] = torch.eye(self.ngrids, device=self.device)            
        self.Selec = params.get("Selec", Selec_default )   

        # Computed with LDR methods
        self.C0 = torch.zeros(self.ndim, dtype=torch.cdouble, device=self.device)
        self.Ccurr = torch.zeros(self.ndim, dtype=torch.cdouble, device=self.device)

        self.Snucl = torch.eye(self.ngrids, dtype=torch.cdouble, device=self.device)
        self.Tnucl = torch.zeros(self.ngrids, self.ngrids, dtype=torch.cdouble, device=self.device)

        self.S, self.H = torch.zeros(self.ndim, self.ndim, dtype=torch.cdouble, device=self.device), torch.zeros(self.ndim, self.ndim, dtype=torch.cdouble, device=self.device)
        self.U = torch.zeros(self.ndim, self.ndim, dtype=torch.cdouble, device=self.device)
        self.Shalf = torch.zeros(self.ndim, self.ndim, dtype=torch.cdouble, device=self.device)
        
        self.time = []
        self.kinetic_energy = []
        self.potential_energy = []
        self.total_energy = []
        self.population_right = []
        self.denmat = []
        self.norm = []
        self.C_save = []

    def chi_overlap(self):
        """
        Compute nuclear overlap matrix Snucl[i, j] for the mesh qmesh 
        from the Gaussian basis, g(x; q) = \exp(-\alpha * (x-q)**2).
        """
        delta = self.qgrid[:, None, :] - self.qgrid[None, :, :]    # (N, N, D)
        exponent = -0.5 * torch.sum(self.alpha * delta**2, dim=2)  # (N, N)
        self.Snucl = torch.exp(exponent)

    def chi_kinetic(self):
        r"""
        Compute nuclear kinetic energy matrix Tnucl[i,j] = <g(x; qgrid[i]) | T | g(x; qgrid[j])>,
        with T = \sum_{\nu} -0.5* m_ν^{-1} \partial^{2}/\partial x_{\nu}^2.
        """
        delta = self.qgrid[:, None, :] - self.qgrid[None, :, :]               # (N, N, D)
        tau = self.alpha / (2.0 * self.mass) * (1.0 - self.alpha * delta**2)  # (N, N, D)
        tau_sum = torch.sum(tau, dim=2)                                       # (N, N)
    
        self.Tnucl = self.Snucl * tau_sum                                     # (N, N)
    
    def build_compound_overlap(self):
        """
        Build the compound nuclear-electronic overlap matrix self.S (ndim, ndim)
        """
        N, s, ndim = self.ngrids, self.nstates, self.ndim
    
        # Reshape Selec[a, b] -> (i, n, j, m) with:
        #   a = i * N + n
        #   b = j * N + m
        Selec4D = self.Selec.view(s, N, s, N) # (i, n, j, m)
    
        Snucl4D = self.Snucl.unsqueeze(0).unsqueeze(2) # (1, n, 1, m)
    
        S4D = Selec4D * Snucl4D
    
        # Reshape back to (ndim, ndim) with compound indices
        self.S = S4D.permute(0, 1, 2, 3).reshape(ndim, ndim)

    def build_compound_hamiltonian(self):
        """
        Build the compound nuclear-electronic Hamiltonian self.H (ndim, ndim) using different schemes.
        """
        N, s, ndim = self.ngrids, self.nstates, self.ndim
        scheme = self.Hamiltonian_scheme
        Selec4D = self.Selec.view(s, N, s, N)       # (s, N, s, N)
        T4D = self.Tnucl.unsqueeze(0).unsqueeze(2)  # (1, N, 1, N)
        S4D = self.Snucl.unsqueeze(0).unsqueeze(2)  # (1, N, 1, N)
    
        if scheme == 'as_is': # For showing the original non-Hermitian form, not intended to use
            Ej4D = self.E[None, None, :, :]           # (1, 1, s, N)
            bracket4D = T4D + Ej4D * S4D
        elif scheme == 'symmetrized':
            Ei4D = self.E[:, :, None, None]   # (s, N, 1, 1)
            Ej4D = self.E[None, None, :, :]   # (1, 1, s, N)
            Eavg4D = 0.5 * (Ei4D + Ej4D)      # (s, N, s, N)
            bracket4D = T4D + Eavg4D * S4D
        elif scheme == 'diagonal':
            # Build Kronecker deltas for electronic and nuclear indices
            delta_ij = torch.eye(s, device=self.device).unsqueeze(1).unsqueeze(3)  # (s, 1, s, 1)
            delta_nm = torch.eye(N, device=self.device).unsqueeze(0).unsqueeze(2)  # (1, N, 1, N)
            delta4D = delta_ij * delta_nm
            
            Ej4D = self.E[None, None, :, :] # (1, 1, s, N)
            bracket4D = T4D + Ej4D * S4D * delta4D

        else:
            raise ValueError(f"Unknown Hamiltonian scheme: {scheme}")
    
        H4D = Selec4D * bracket4D
        self.H = H4D.reshape(ndim, ndim)

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
        self.Ccurr = self.C0.clone()

        print(F"step = 0")
        self.save_results(0)
        
        for step in range(1, self.nsteps):
            Cvec = self.Ccurr.clone()
            self.Ccurr = self.U @ Cvec

            if step % self.save_every_n_steps == 0:
                print(F"step = {step}")
                self.save_results(step)

    def save_results(self, step):
        if "time" in self.properties_to_save:
            self.time.append(step*self.dt)
        if "norm" in self.properties_to_save:
            overlap = torch.matmul(self.S, self.Ccurr)
            self.norm.append(torch.sqrt(torch.vdot(self.Ccurr, overlap)))
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
        if "C_save" in self.properties_to_save:
            self.C_save.append(self.Ccurr)
        
    def compute_populations(self):
        """
        Compute electronic state population for a single step.
        """
        N, s = self.ngrids, self.nstates
        Cvec = self.Ccurr
        
        # Compute SC once: shape (ndim,)
        SC = self.S @ Cvec
    
        C_blocks = Cvec.view(s, N)
        SC_blocks = SC.view(s, N)
    
        # Compute P[i] = sum_j <C_j|S_{ji}|C_i> = Re[ sum_N (C_j*) * SC_j ]
        P = torch.sum(C_blocks.conj() * SC_blocks, dim=1).real
    
        return P
    
    def compute_denmat_raw(self):
        """
        Compute electronic density matrix for a single step without the orthogonalization.
        """
        N, s = self.ngrids, self.nstates
        Cvec = self.Ccurr
    
        # Compute SC once: shape (ndim,)
        SC = self.S @ Cvec
    
        C_blocks = Cvec.view(s, N)
        SC_blocks = SC.view(s, N)
    
        # \rho_ij = \sum_n C_i,n^* (SC)_j,n
        rho = SC_blocks.conj() @ C_blocks.T 
    
        return rho

    def compute_denmat(self):
        """
        Compute electronic density matrix for a single step without the orthogonalization.
        """
        N, s = self.ngrids, self.nstates
        Cvec = self.Ccurr
    
        # Orthogonalize coefficients: C_ortho = S^{1/2} C
        C_ortho = self.S_half @ Cvec
      
        C_blocks = C_ortho.view(s, N)
    
        rho = C_blocks @ C_blocks.conj().T # (s, s)
    
        return rho

    def compute_kinetic_energy(self):
        """
        Compute nuclear kinetic energy as C^+ T C / C^+ S C for a single step.
        """
        N, s, ndim = self.ngrids, self.nstates, self.ndim
    
        # Rebuild compound kinetic matrix: T4D * Selec4D
        Selec4D = self.Selec.view(s, N, s, N)
        T4D = self.Tnucl.unsqueeze(0).unsqueeze(2)  # (1, n, 1, m)
        T4D_compound = Selec4D * T4D
        T_compound = T4D_compound.permute(0, 1, 2, 3).reshape(ndim, ndim)
    
        Cvec = self.Ccurr
    
        numer = torch.vdot(Cvec, T_compound @ Cvec).real
        denom = torch.vdot(Cvec, self.S @ Cvec).real
    
        return numer / denom
    
    
    def compute_potential_energy(self):
        """
        Compute potential energy as C^+ V C / C^+ S C for a single step.
        """
        N, s, ndim = self.ngrids, self.nstates, self.ndim
    
        Selec4D = self.Selec.view(s, N, s, N)
        S4D = self.Snucl.unsqueeze(0).unsqueeze(2)  # (1, n, 1, m)
        Ej4D = self.E[None, None, :, :]  # (1,1,j,m)
    
        V4D_compound = Selec4D * (Ej4D * S4D)
        V_compound = V4D_compound.permute(0, 1, 2, 3).reshape(ndim, ndim)
    
        Cvec = self.Ccurr
    
        numer = torch.vdot(Cvec, V_compound @ Cvec).real
        denom = torch.vdot(Cvec, self.S @ Cvec).real
    
        return numer / denom
    
    
    def compute_total_energy(self):
        """
        Compute total energy as C^+ H C / C^+ S C for a single step.
        """
        Cvec = self.Ccurr
    
        numer = torch.vdot(Cvec, self.H @ Cvec).real
        denom = torch.vdot(Cvec, self.S @ Cvec).real
    
        return numer / denom
        
    def save(self):
        torch.save( {"q0":self.q0,
                     "p0":self.p0,
                     "k":self.k,
                     "mass":self.mass,
                     "alpha":self.alpha,
                     "qgrid":self.qgrid,
                     "nstates":self.nstates,
                     "istate":self.istate,
                     "Snucl":self.Snucl,
                     "Tnucl":self.Tnucl,
                     "E":self.E,
                     "Selec":self.Selec,
                     "S":self.S,
                     "H":self.H,
                     "U":self.U,
                     "C_save":self.C_save,
                     "save_every_n_steps":self.save_every_n_steps,
                     "Hamiltonian_scheme": self.Hamiltonian_scheme,
                     "dt":self.dt, "nsteps":self.nsteps,
                     "time":self.time,
                     "kinetic_energy":self.kinetic_energy,
                     "potential_energy":self.potential_energy,
                     "total_energy":self.total_energy,
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

