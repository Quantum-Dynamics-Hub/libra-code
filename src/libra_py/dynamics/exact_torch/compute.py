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
        

