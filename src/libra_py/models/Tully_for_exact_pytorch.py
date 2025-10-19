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
.. module:: Tully_for_exact_pytorch
   :platform: Unix, Windows
   :synopsis: This module implements 3 Tully models written in PyTorch way for the use 
              with the exact DVR calculations. This format is compatible with Libra's DVR
              implementation
.. moduleauthor:: Alexey V. Akimov

"""

import torch

# Define Tully's simple avoided crossing diabatic potential matrix
def tully_model1_simple_avoided_crossing(Q, params):
    """
    Q:tensor( [ndof, N_1, N_2, ... N_ndof] )
    Returns diabatic potential matrix [N_1, N_2, ... N_ndof, 2, 2]
    """
    x = Q[0]  # Assume 1D nuclear coordinate

    A = params.get("A", 0.01)
    B = params.get("B", 1.6)
    C = params.get("C", 0.005)
    D = params.get("D", 1.0)


    # Diabatic state 1 potential
    V11 = torch.where(
    x >= 0,
    A * (1 - torch.exp(-B * x)),
    -A * (1 - torch.exp(B * x))
    )
    V22 = -V11                         # Diabatic state 2 potential (mirror)
    V12 = C * torch.exp(-D * x**2)    # Coupling between diabatic states

    shape = x.shape
    Vmat = torch.zeros((*shape, 2, 2), dtype=torch.cfloat)
    Vmat[..., 0, 0] = V11
    Vmat[..., 1, 1] = V22
    Vmat[..., 0, 1] = V12
    Vmat[..., 1, 0] = torch.conj(V12)

    return Vmat

def tully_model2_dual_avoided_crossing(Q, params):
    """
    Tully Model 2: Dual Avoided Crossing (diabatic representation)

    Q:tensor( [ndof, N_1, N_2, ... N_ndof] )
    Returns: diabatic potential matrix with shape [Ngrid, 2, 2]
             (same trailing [2,2] convention as in your model-1 function)
    """
    x = Q[0]  # Assume 1D nuclear coordinate

    # Standard parameter set used in Tully's paper
    A  = params.get("A", 0.10)
    B  = params.get("B", 0.28)
    C  = params.get("C", 0.015)
    D  = params.get("D", 0.06)
    E0 = params.get("E0", 0.05)  # Energy offset on state 2

    # Diabatic potentials
    V11 = torch.zeros_like(x)                          # V_11(x) = 0
    V22 = -A * torch.exp(-B * x**2) + E0              # V_22(x)
    V12 = C * torch.exp(-D * x**2)                    # Coupling V_12(x) = V_21(x)

    shape = x.shape
    Vmat = torch.zeros((*shape, 2, 2), dtype=torch.cfloat, device=x.device)
    Vmat[..., 0, 0] = V11
    Vmat[..., 1, 1] = V22
    Vmat[..., 0, 1] = V12
    Vmat[..., 1, 0] = torch.conj(V12)

    return Vmat

def tully_model3_extended_coupling(Q, params):
    """
    Tully Model 3: Extended Coupling with Reflection (diabatic representation)

    Q: Tensor with shape [ndof, Ngrid]
    Returns: diabatic potential matrix with shape [Ngrid, 2, 2]
    """
    x = Q[0]  # 1D nuclear coordinate

    # Canonical parameters from Tully (1990)
    A = params.get("A", 6.0e-4)
    B = params.get("B", 0.10)
    C = params.get("C", 0.90)

    # Diabatic diagonal elements: constant ±A
    V11 = torch.full_like(x, A)
    V22 = -V11

    # Off-diagonal coupling: piecewise in x
    V12 = torch.where(
        x < 0,
        B * torch.exp(C * x),
        B * (2.0 - torch.exp(-C * x))
    )

    # Assemble Hamiltonian
    shape = x.shape
    Vmat = torch.zeros((*shape, 2, 2), dtype=torch.cfloat, device=x.device)
    Vmat[..., 0, 0] = V11
    Vmat[..., 1, 1] = V22
    Vmat[..., 0, 1] = V12
    Vmat[..., 1, 0] = torch.conj(V12)

    return Vmat


