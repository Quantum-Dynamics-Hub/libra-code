# *********************************************************************************
# * Copyright (C) 2025 Alexey V. Akimov
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
   :synopsis: This module implements functions for doing exact on-the-grid dynamics with PyTorch
       List of functions:
           * sech # temporary here
           * Martens_model # temporary here
           * gaussian_wavepacket
       List of classes:
           * exact_tdse_solver

.. moduleauthor:: Alexey V. Akimov

"""

__author__ = "Alexey V. Akimov"
__copyright__ = "Copyright 2025 Alexey V. Akimov"
__credits__ = ["Alexey V. Akimov"]
__license__ = "GNU-3"
__version__ = "1.0"
__maintainer__ = "Alexey V. Akimov"
__email__ = "alexvakimov@gmail.com"
__url__ = "https://github.com/Quantum-Dynamics-Hub/libra-code"



import torch
import torch.fft
import numpy as np


def sech(x):
  return 1 / torch.cosh(x)

def Martens_model(q, params):
    """
    q - Tensor(ndof)

    Martens_model1 is just this one but with Vc = 0.0
    """
    #params = {"Va": 0.00625, "Vb": 0.0106}
    Va = params.get("Va", 0.00625)
    Vb = params.get("Vb", 0.0106)
    Vc = params.get("Vc", 0.0)
    return Va * (sech(2.0*q[0]))**2 + 0.5 * Vb * (q[1] + Vc * (q[0]**2 - 1.0 ) )**2


# Define Tully's simple avoided crossing diabatic potential matrix
def tully_potential_matrix(Q, params):
    """
    Q: Tensor with shape [ndof, Ngrid]
    Returns diabatic potential matrix [2, 2, Ngrid]
    """
    x = Q[0]  # Assume 1D nuclear coordinate
    
    A = params.get("A", 0.01)
    B = params.get("B", 1.6)
    C = params.get("C", 0.005)
    D = params.get("D", 1.0)
   

    V11 = torch.where(
    x >= 0,
    A * (1 - torch.exp(-B * x)),
    -A * (1 - torch.exp(B * x))
    )
 
    #V11 = A * (1 - torch.exp(-B * x))  # Diabatic state 1 potential
    V22 = -V11                         # Diabatic state 2 potential (mirror)
    V12 = C * torch.exp(-D * x**2)    # Coupling between diabatic states
    
    shape = x.shape
    Vmat = torch.zeros((*shape, 2, 2), dtype=torch.cfloat)
    Vmat[..., 0, 0] = V11
    Vmat[..., 1, 1] = V22
    Vmat[..., 0, 1] = V12
    Vmat[..., 1, 0] = torch.conj(V12)

    return Vmat


def gaussian_wavepacket(q, params):
    """
    q = tensor( [ndof, N_1, N_2, ... N_ndof] )
    
    """
    hbar = 1.0
    ndof = q.shape[0]
    sz = len(q.shape)
    
    mass = torch.tensor(params.get("mass", [2000.0, 2000.0]) )
    omega = torch.tensor(params.get("omega", [0.004, 0.004]) )
    sigma = 1.0 / torch.sqrt( 2.0 * mass * omega )
    q0 = torch.tensor(params.get("q0", [-1.0, 0.0]) )
    p0 = torch.tensor(params.get("p0", [3.0, 0.0]) )

    # Reshape q, p, sigma to be compatible in shape with q
    sigma = sigma.view(ndof, *[1]*(sz - 1) )
    q0 = q0.view(ndof, *[1]*(sz - 1) )
    p0 = p0.view(ndof, *[1]*(sz - 1) )

    # Do the calculations:        
    phase = 1j * torch.sum(p0 * q, dim=0, keepdim=False) / hbar
    envelope = torch.exp(-torch.sum( 0.25*(q - q0)**2 / sigma**2 , dim=0, keepdim=False) )  # because it is wavefunction
    norm = torch.prod(1.0 / ( (sigma * (2 * torch.pi)**0.5 )**0.5 ) , dim=0, keepdim=False)
    
    return (norm * envelope * torch.exp(phase))  # this also reduces the first dimension of q



class exact_tdse_solver:
    def __init__(self, params):
        self.prefix = params.get("prefix", "exact-solution")
        self.grid_size = torch.tensor(params.get("grid_size", [4, 4]))
        self.ndim = len(self.grid_size)
        self.q_min = torch.tensor(params.get("q_min", [-10.0] * self.ndim))
        self.q_max = torch.tensor(params.get("q_max", [10.0] * self.ndim))
        self.save_every_n_steps = params.get("save_every_n_steps", 1)
        self.dt = params.get("dt", 0.01)
        self.nsteps = params.get("nsteps", 500)
        self.mass = torch.tensor(params.get("mass", [1.0] * self.ndim))
        self.potential_fn = params.get("potential_fn", None)
        self.potential_fn_params = params.get("potential_fn_params", None)
        self.psi0_fn = params.get("psi0_fn", None)        
        self.psi0_fn_params = params.get("psi0_fn_params", None)
        self.psi = None
        self.psi_k = None
        self.prob_density = None        
        self.device = params.get("device", torch.device("cuda" if torch.cuda.is_available() else "cpu"))
        self.hbar = 1.0
        self.psi_all = None
        self.time = []
        self.kinetic_energy = []
        self.potential_energy = []
        self.total_energy = []
        self.population_right = []
        self.norm = []
    

    def initialize_grids(self):
        print("Initializing grids")
        print("grid_size = ", self.grid_size)
        # Real-space grid
        q_axes = [torch.linspace(self.q_min[i], self.q_max[i], self.grid_size[i]) for i in range(self.ndim)]
        q_grids = torch.meshgrid(*q_axes, indexing="ij")
        self.dq = torch.tensor([q_axes[i][1] - q_axes[i][0] for i in range(self.ndim)])
        self.Q = torch.stack(q_grids)
        print("Q grid = ", self.Q.shape)
        print("dq = ", self.dq)

        # Momentum-space grid        
        self.dk = 2 * torch.pi / (self.grid_size * self.dq)
        k_axes = [torch.fft.fftshift(torch.arange(-self.grid_size[i] // 2, self.grid_size[i] // 2)) * self.dk[i] for i in range(self.ndim)]
        k_grids = torch.meshgrid(*k_axes, indexing="ij")
        self.K = torch.stack(k_grids)
        print("K grid = ", self.K.shape)
        print("dk = ", self.dk)

        # Volume elements
        self.dV = self.dq.prod()
        self.dVk = self.dk.prod()  #(self.dq / self.grid_size).prod()
        print("dV = ", self.dV)
        print("dVk = ", self.dVk)

        size = [ int((self.nsteps - self.nsteps % self.save_every_n_steps ) / self.save_every_n_steps)+1 ]
        for i in range(self.ndim):
            size.append(self.grid_size[i])
        self.psi_all = torch.zeros(size, dtype=torch.cfloat)

    def initialize_operators(self):
        view_shape = [self.ndim] + [1] * self.ndim
        self.T = 0.5 * torch.sum((self.hbar * self.K) ** 2 / self.mass.view(view_shape), dim=0)
        self.V = self.potential_fn(self.Q, self.potential_fn_params) if self.potential_fn else torch.zeros_like(self.T)
        self.psi = self.psi0_fn(self.Q, self.psi0_fn_params) if self.psi0_fn else torch.zeros_like(self.V)
        self.psi_k = torch.fft.fftn(self.psi) #, norm='forward')
        self.psi_all[0] = self.psi
        self.expV_half = torch.exp(-0.5j * self.V * self.dt / self.hbar)
        self.expT = torch.exp(-1j * self.T * self.dt / self.hbar)

    def propagate(self):        
        for step in range(self.nsteps):

            if step % self.save_every_n_steps == 0:
                istep = int(step / self.save_every_n_steps)
                self.psi_all[istep] = self.psi           
                self.prob_density = torch.abs(self.psi) ** 2            
                self.psi_k = torch.fft.fftn(self.psi) # norm='forward')
                KE = torch.sum(torch.abs(self.psi_k) ** 2 * self.T) * (self.dq/self.grid_size).prod()
                PE = torch.sum(self.prob_density * self.V) * self.dV  
                nrm = torch.sum(torch.abs(self.psi) ** 2 ) * self.dV            
                x_coords = self.Q[0]
                right_mask = self.Q[0] > 0
                pop_right = torch.sum(self.prob_density[right_mask]) * self.dV

                self.norm.append( nrm )
                self.time.append(step * self.dt)
                self.kinetic_energy.append( KE.real.item() )
                self.potential_energy.append( PE.real.item() )
                self.total_energy.append( KE + PE )
                self.population_right.append(pop_right.item())
                        
                print(f"Step {step}: KE = {KE:.4f}, PE = {PE:.4f}, Total = {KE + PE:.4f}, Norm = {nrm:.4f}")

            #============== Propagate ===================
            self.psi *= self.expV_half
            self.psi_k = torch.fft.fftn(self.psi) # norm='forward')
            self.psi_k *= self.expT
            self.psi = torch.fft.ifftn(self.psi_k) # norm='forward')
            self.psi *= self.expV_half

    def save(self):
        torch.save( {"grid_size":self.grid_size, 
                     "ndim":self.ndim, 
                     "q_min":self.q_min, "q_max":self.q_max, 
                     "save_every_n_steps": self.save_every_n_steps,
                     "dt":self.dt, "nsteps":self.nsteps,
                     "mass":self.mass,
                     "psi":self.psi, "psi_k":self.psi_k,
                     "prob_density":self.prob_density,
                     "psi_all":self.psi_all,
                     "time":self.time,
                     "Q":self.Q, "K":self.K, "dq":self.dq, "dk":self.dk,
                     "dV":self.dV, "dVk":self.dVk,
                     "kinetic_energy":self.kinetic_energy,
                     "potential_energy":self.potential_energy,
                     "total_energy":self.total_energy,
                     "population_right":self.population_right,
                     "norm":self.norm,
                     "V":self.V, "T":self.T
                    }, F"{self.prefix}.pt" )
    
    def solve(self):
        self.initialize_grids()
        self.initialize_operators()
        self.propagate()
        self.save()
        



class exact_tdse_solver_multistate:
    def __init__(self, params):
        """
        Initializes the TDSE solver.

        Parameters:
            params (dict): Dictionary of simulation parameters:
                - prefix (str): Filename prefix for output
                - grid_size (list[int]): Number of points per spatial dimension
                - q_min (list[float]), q_max (list[float]): Spatial bounds
                - save_every_n_steps (int): Interval for recording data
                - dt (float): Time step
                - nsteps (int): Number of time steps
                - mass (list[float]): Masses per dimension
                - potential_fn_params (dict): Parameters for potential energy surface
                - psi0_fn_params (dict): Parameters for initial wavepacket
                - device (torch.device): Computation device
                - Nstates (int): Number of electronic states
                - representation (str): "diabatic" or "adiabatic"
                - initial_state_index (int): Which state to initialize
                - method (str): Propagation scheme ('miller-colton', 'split-operator', 'crank-nicolson')
        """

        self.prefix = params.get("prefix", "exact-solution")
        self.grid_size = torch.tensor(params.get("grid_size", [64, 64]))
        self.ndim = len(self.grid_size)
        self.q_min = torch.tensor(params.get("q_min", [-10.0]*self.ndim))
        self.q_max = torch.tensor(params.get("q_max", [10.0]*self.ndim))
        self.save_every_n_steps = params.get("save_every_n_steps", 1)
        self.dt = params.get("dt", 0.01)
        self.nsteps = params.get("nsteps", 500)
        self.mass = torch.tensor(params.get("mass", [2000.0]*self.ndim))
        self.potential_fn = params.get("potential_fn", None)
        self.potential_fn_params = params.get("potential_fn_params", {})
        self.psi0_fn = params.get("psi0_fn", None)        
        self.psi0_fn_params = params.get("psi0_fn_params", {})
        self.device = params.get("device", torch.device("cuda" if torch.cuda.is_available() else "cpu"))
        self.hbar = 1.0
        self.Nstates = params.get("Nstates", 2)
        self.representation = params.get("representation", "diabatic")
        self.initial_state_index = params.get("initial_state_index", 0)
        self.method = params.get("method", "miller-colton").lower()

        self.time = []
        self.kinetic_energy = []
        self.potential_energy = []
        self.total_energy = []
        self.population_right = []
        self.norm = []

    def initialize_grids(self):
        """
        Constructs real- and momentum-space grids, initializes kinetic energy operator
        and wavefunction on grid in chosen representation.
        """
        print("Initializing grids")
        print("grid_size = ", self.grid_size)
        self.ngrid = int(self.grid_size.prod().item())
        print("self.ngrid = ", self.ngrid)

        # Real-space grid
        q_axes = [torch.linspace(self.q_min[i], self.q_max[i], self.grid_size[i]) for i in range(self.ndim)]
        q_grids = torch.meshgrid(*q_axes, indexing="ij")
        self.dq = torch.tensor([q_axes[i][1] - q_axes[i][0] for i in range(self.ndim)])
        self.Q = torch.stack(q_grids)
        print("Q = ", self.Q.shape)
        print("dq = ", self.dq)

        # Momentum-space grid        
        self.dk = 2 * torch.pi / (self.grid_size * self.dq)
        k_axes = [torch.fft.fftshift(torch.arange(-self.grid_size[i] // 2, self.grid_size[i] // 2)) * self.dk[i] for i in range(self.ndim)]
        k_grids = torch.meshgrid(*k_axes, indexing="ij")
        self.K = torch.stack(k_grids)
        print("K = ", self.K.shape)
        print("dk = ", self.dk)

        # Volume elements
        self.dV = self.dq.prod()
        self.dVk = self.dk.prod()  #(self.dq / self.grid_size).prod()
        print("dV = ", self.dV)
        print("dVk = ", self.dVk)
    

        # Allocate storage:
        self.nsnaps = self.nsteps // self.save_every_n_steps + 1  # how many snapshots to save

        # Diabatic properties
        self.psi_r_dia = torch.zeros((*self.grid_size, self.Nstates), dtype=torch.cfloat)
        self.psi_k_dia = torch.zeros_like(self.psi_r_dia)
        self.psi_r_dia_all = torch.zeros(( self.nsnaps , *self.grid_size, self.Nstates), dtype=torch.cfloat)
        self.psi_k_dia_all = torch.zeros((self.nsnaps, *self.grid_size, self.Nstates), dtype=torch.cfloat)
        self.rho_dia_all = torch.zeros((self.nsnaps, self.Nstates, self.Nstates), dtype=torch.cfloat)

        # Adiabatic properties
        self.psi_r_adi = torch.zeros((*self.grid_size, self.Nstates), dtype=torch.cfloat)
        self.psi_k_adi = torch.zeros_like(self.psi_r_adi)
        self.psi_r_adi_all = torch.zeros(( self.nsnaps , *self.grid_size, self.Nstates), dtype=torch.cfloat)
        self.psi_k_adi_all = torch.zeros((self.nsnaps, *self.grid_size, self.Nstates), dtype=torch.cfloat)
        self.rho_adi_all = torch.zeros((self.nsnaps, self.Nstates, self.Nstates), dtype=torch.cfloat)


    def update_adi_r(self):
        """
        Convert diabatic to adiabatic r-space wavefunction: C_adi = U * C_dia
        """
        self.psi_r_adi = torch.einsum("...ij, ...j->...i", self.eigvecs, self.psi_r_dia)

    def update_dia_r(self):
        """
        Convert adiabatic to diabatic r-space wavefunction: C_dia = U.H * C_adi
        """
        self.psi_r_dia = torch.einsum("...ij, ...j->...i", self.eigvecs.conj().transpose(-2,-1), self.psi_r_adi)


    def transform_r2k(self, rep):
        """
        Compute the k-space diabatic and adiabatic wavefunctions from the r-space counterparts
        rep: 0 - diabatic, 1 - adiabatic
        """
        # Define which axes are spatial (everything except last one = Nstates)
        spatial_dims = tuple(range(self.ndim))

        # Apply FFT over spatial dimensions for all states at once
        if rep==0:
            self.psi_k_dia = torch.fft.fftn(self.psi_r_dia, dim=spatial_dims ) # norm='forward')
        elif rep==1:
            self.psi_k_adi = torch.fft.fftn(self.psi_r_adi, dim=spatial_dims ) # norm='forward')

    def transform_k2r(self, rep):
        """
        Compute the k-space diabatic and adiabatic wavefunctions from the r-space counterparts
        rep: 0 - diabatic, 1 - adiabatic
        """
        # Define which axes are spatial (everything except last one = Nstates)
        spatial_dims = tuple(range(self.ndim))

        # Apply inverse FFT over same spatial dimensions
        if rep==0:
            self.psi_r_dia = torch.fft.ifftn(self.psi_k_dia, dim=spatial_dims ) # norm='forward')
        elif rep==1:
            self.psi_r_adi = torch.fft.ifftn(self.psi_k_adi, dim=spatial_dims ) # norm='forward')

        
    def initialize_operators(self):        
        # Kinetic energy operator - same for all states
        self.T = 0.5 * torch.sum((self.hbar * self.K) ** 2 / self.mass.view(self.ndim, *[1]*self.ndim), dim=0)  # T: [*grid_size]
        print("T.shape = ", self.T.shape)

        # Element-wise exponentials are fine here
        self.expT = torch.exp(-1j * self.T * self.dt / self.hbar) # expT: [*grid_size]
        print("expT.shape = ", self.expT.shape)
        
        # Potential energy operator - matrix for multiple states
        self.V = self.potential_fn(self.Q, self.potential_fn_params) # V: [*grid_size, Nstates, Nstates]
        print("V.shape = ", self.V.shape)

        # dia  <-> adi transformation for all points
        # V U = E U => H = U.H * E * U
        # E = <psi_adi | H | psi_adi >
        # V = <psi_dia | H | psi_dia >
        # | psi_adi > U = | psi_dia >
        # | Psi > = | psi_dia> C_dia = | psi_adi > C_adi = | psi_adi > U C_dia
        # so C_adi = U C_dia 

        self.eigvals, self.eigvecs = torch.linalg.eigh(self.V)  # V: [*grid_size, Nstates, Nstates]
        self.exp_diag = torch.exp(-0.5j * self.dt * self.eigvals)
        self.expV_half = self.eigvecs.conj().transpose(-2, -1) @ torch.diag_embed(self.exp_diag) @ self.eigvecs
        print("expV_half.shape = ", self.expV_half.shape)

        # Initialize the wavefunction in r-space
        wfn = self.psi0_fn(self.Q, self.psi0_fn_params)  

        if self.representation == "diabatic":
            self.psi_r_dia[..., self.initial_state_index] = wfn
            self.update_adi_r()  # dia -> adi
        elif self.representation == "adiabatic":
            self.psi_r_adi[..., self.initial_state_index] = wfn
            self.update_dia_r()  # adi -> dia
        #print("psi.shape = ", self.psi_r_dia.shape)
 
        # Update the k-space wavefunctions:
        self.transform_r2k(0) # psi_r_dia -> psi_k_dia
        self.transform_r2k(1) # psi_r_adi -> psi_k_adi


    def propagate(self):
        """
        Time-propagates the wavefunction using the selected propagation method.
        Supported methods: 'miller-colton', 'split-operator', 'crank-nicolson'
        """
        for step in range(self.nsteps):
        
            #====================== Saving and computing properties ==================
        
            if step % self.save_every_n_steps == 0:
                istep = int(step / self.save_every_n_steps)
               
                # Diabatic r-space wavefunctions
                self.psi_r_dia_all[istep] = self.psi_r_dia

                # Compute other representations:
                self.update_adi_r();   self.psi_r_adi_all[istep] = self.psi_r_adi; # adiabatic in r-space
                self.transform_r2k(0); self.psi_k_dia_all[istep] = self.psi_k_dia; # diabatic in k-space
                self.transform_r2k(1); self.psi_k_adi_all[istep] = self.psi_k_adi; # adiabatic in k-space
                             

                # Diabatic density matrix
                self.rho_dia_all[istep] = torch.einsum("...i, ...j->ij", self.psi_r_dia, self.psi_r_dia.conj() ) * self.dV

                # Adiabatic density matrix
                self.rho_adi_all[istep] = torch.einsum("...i, ...j->ij", self.psi_r_adi, self.psi_r_adi.conj() )  * self.dV

                # Kinetic energy
                KE = torch.einsum("...i,...,...i->", self.psi_k_dia.conj(), self.T, self.psi_k_dia ) * (self.dq/self.grid_size).prod()

                # Full potential energy: PE = ∫ ψ* V ψ dx
                PE = torch.einsum("...i,...ij,...j->", self.psi_r_dia.conj(), self.V, self.psi_r_dia) * self.dV


                nrm = torch.sum(torch.abs(self.psi_r_dia) ** 2 ) * self.dV            
                x_coords = self.Q[0]
                right_mask = self.Q[0] > 0
                #pop_right = torch.sum(self.prob_density[right_mask]) * self.dV
                
                self.norm.append( nrm )
                self.time.append(step * self.dt)
                self.kinetic_energy.append( KE.real.item() )
                self.potential_energy.append( PE.real.item() )
                self.total_energy.append( KE + PE )
                #self.population_right.append(pop_right.item())
                        
                print(f"Step {step}: Norm = {nrm:.4f}")
        
        
            #=================== Doing computations ==========================        
            if self.method == "crank-nicolson":
                # Use Crank-Nicolson scheme for time propagation
                # self.crank_nicolson_step()
                pass # not implemented

            else:
                # Apply half step of potential evolution in adiabatic representation
                # self.potential_half_step()

                if self.method == "miller-colton":
                    pass # not  implemented
                    # Miller-Colton propagation in momentum space
                    #for i in range(self.Nstates):
                    #    self.psi_k[i] = torch.fft.fftn(self.psi[i])
                    #    self.psi_k[i] *= torch.exp(-1j * self.T * self.dt / (2 * self.hbar))
                    #    self.psi[i] = torch.fft.ifftn(self.psi_k[i])

                elif self.method == "split-operator":
                    # Half-step propagation in real space
                    self.psi_r_dia = torch.einsum("...ij,...j->...i", self.expV_half, self.psi_r_dia) #self.expV_half @ self.psi 

                    # Full-step in reciprocal space
                    self.transform_r2k(0) 

                    # Multiply by expT — make sure it broadcasts correctly
                    self.psi_k_dia *= self.expT.unsqueeze(-1)  # if expT.shape == [*,], expand to [*, 1]

                    # Apply inverse FFT over same spatial dimensions
                    self.transform_k2r(0)

                    # Another half of the propagation in real space
                    self.psi_r_dia = torch.einsum("...ij,...j->...i", self.expV_half, self.psi_r_dia) #self.expV_half @ self.psi


    def save(self):
        """Saves the grid, wavefunction, and observables to disk."""
        torch.save( {"grid_size":self.grid_size, 
                     "ndim":self.ndim, 
                     "q_min":self.q_min, "q_max":self.q_max, 
                     "save_every_n_steps": self.save_every_n_steps,
                     "dt":self.dt, "nsteps":self.nsteps,
                     "mass":self.mass,
                     "psi_r_adi":self.psi_r_adi,
                     "psi_r_dia":self.psi_r_dia,
                     "psi_k_adi":self.psi_k_adi,
                     "psi_k_dia":self.psi_k_dia,
                     "rho_dia_all":self.rho_dia_all,
                     "rho_adi_all":self.rho_adi_all,
                     "psi_r_dia_all":self.psi_r_dia_all,
                     "psi_r_adi_all":self.psi_r_adi_all,
                     "psi_k_dia_all":self.psi_k_dia_all,
                     "psi_k_adi_all":self.psi_k_adi_all,
                     "time":self.time,
                     "Q":self.Q, "K":self.K, "dq":self.dq, "dk":self.dk,
                     "dV":self.dV, "dVk":self.dVk,
                     "kinetic_energy":self.kinetic_energy,
                     "potential_energy":self.potential_energy,
                     "total_energy":self.total_energy,
                     "population_right":self.population_right,
                     "norm":self.norm,
                     "V":self.V, "T":self.T,
                    }, F"{self.prefix}.pt" )
                    
    def solve(self):
        """Runs the full simulation and saves results."""
        self.initialize_grids()
        self.initialize_operators()
        self.propagate()
        self.save()


