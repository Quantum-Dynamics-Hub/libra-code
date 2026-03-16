# *********************************************************************************
# * Copyright (C) 2026 Vladimir Mandelshtam and Alexey V. Akimov
# *
# * This file is distributed under the terms of the GNU General Public License
# * as published by the Free Software Foundation, either version 3 of
# * the License, or (at your option) any later version.
# * See the file LICENSE in the root directory of this distribution
# * or <http://www.gnu.org/licenses/>.
# *
# *********************************************************************************/

"""
.. module:: grids
   :platform: Unix, Windows
   :synopsis: this module implements construction of different kinds of grids for quantum dynamics calculations

.. moduleauthor:: Vladimir Mandelshtam, Alexey V. Akimov, ChatGPT

"""


import torch


def pair_rpl_vec(r_all, r_trial, sigma_all, sigma_trial, rrange1):
    """
    Compute vectorized pairwise repulsive energy between all grid points
    and a trial point.

    This function evaluates the repulsive interaction

        E = (sigma_i / r_ij)^(d + 9) + (sigma_trial / r_ij)^(d + 9)

    for all grid points simultaneously using PyTorch tensor operations.

    Parameters
    ----------
    r_all : torch.Tensor, shape (d1, NG)
        Coordinates of all grid points.

    r_trial : torch.Tensor, shape (d1,)
        Coordinates of the trial grid point.

    sigma_all : torch.Tensor, shape (NG,)
        Sigma values associated with each grid point.

    sigma_trial : float or torch.Tensor
        Sigma value corresponding to the trial point.

    rrange1 : torch.Tensor, shape (d1,)
        Scaling factors applied to coordinate differences when computing
        distances.

    Returns
    -------
    torch.Tensor, shape (NG,)
        Repulsive interaction energy between each grid point and the trial point.

    Notes
    -----
    Distance is computed as

        r = sqrt( sum_i ( (r_i - r_trial_i) * rrange1_i )^2 )

    A small numerical offset (1e-12) is added to avoid singularities
    when r → 0.
    """

    d1 = r_all.shape[0]

    diff = (r_all - r_trial[:, None]) * rrange1[:, None]
    dist = torch.sqrt(torch.sum(diff**2, dim=0) + 1e-12)

    energy = (sigma_all / dist) ** (d1 + 9) + (sigma_trial / dist) ** (d1 + 9)

    return energy


def construct_grid_qrg_fast(N_MC, rg, V_, P, rrange1, device="cpu", step=0.01,):
    """
    Optimize a quantum reference grid (QRG) using Monte Carlo sampling.

    This routine performs a stochastic optimization of grid points by
    minimizing the total pairwise repulsive energy between them.
    The algorithm is equivalent to the Fortran routine `Construct_grid_QRG`.

    The optimization proceeds by repeatedly proposing random displacements
    of individual grid points and accepting only energy-lowering moves
    (zero-temperature Monte Carlo).

    Parameters
    ----------
    N_MC : int
        Number of Monte Carlo sweeps per grid point.
        Total number of trial moves = NG * N_MC.

    rg : torch.Tensor, shape (d1, NG)
        Coordinates of grid points.

    V_ : torch.Tensor, shape (NG,)
        Potential energy at each grid point.

    P : callable
        Function evaluating the probability density at a coordinate.
        Should return

            P_value, V_value = P(r)

        where
        - P_value : probability density
        - V_value : potential energy

    rrange1 : torch.Tensor, shape (d1,)
        Coordinate scaling factors used in distance calculations.

    device : str, optional
        Device for computation ('cpu' or 'cuda').

    step : double 
        Step size for MC steps

    Returns
    -------
    rg : torch.Tensor, shape (d1, NG)
        Optimized grid coordinates.

    sigma : torch.Tensor, shape (NG,)
        Sigma values computed from the probability density.

    V_ : torch.Tensor, shape (NG,)
        Potential energy values at grid points.

    Notes
    -----
    Sigma values are computed as

        sigma_i = P(r_i)^(-1/d)

    where d is the dimensionality of the coordinate space.

    The repulsive interaction between two points is

        E_ij = (sigma_i / r_ij)^(d + 9) + (sigma_j / r_ij)^(d + 9)

    Only moves that reduce total energy are accepted:

        accept if ΔE <= 0

    The step size of trial moves is automatically adjusted
    to maintain an acceptance ratio between 40% and 60%.

    Performance
    -----------
    The energy change calculation is vectorized using PyTorch
    broadcasting, eliminating Python loops and enabling GPU
    acceleration for large grids.

    This reduces the runtime significantly for large NG.
    """

    rg = rg.to(device)
    V_ = V_.to(device)

    d1, NG = rg.shape

    sigma = torch.zeros(NG, device=device)

    # Initialize sigma
    for i in range(NG):
        P_val, V_val = P(rg[:, i])
        sigma[i] = P_val ** (-1.0 / d1)
        V_[i] = V_val

    print("Start optimizing QRG")

    accept = 0
    #step = 0.01

    for n in range(NG * N_MC):

        # choose random point
        k = torch.randint(0, NG, (1,)).item()

        delr = torch.rand(d1, device=device)
        r_trial = rg[:, k] + step * (2 * delr - 1)

        P_trial, V_trial = P(r_trial)

        if P_trial > 1e-20:

            sigma_trial = P_trial ** (-1.0 / d1)

            # vectorized interaction energy
            E_trial = pair_rpl_vec(rg, r_trial, sigma, sigma_trial, rrange1)
            E_old = pair_rpl_vec(rg, rg[:, k], sigma, sigma[k], rrange1)

            # remove self interaction
            E_trial[k] = 0.0
            E_old[k] = 0.0

            delE = torch.sum(E_trial - E_old)

            if delE <= 0:

                rg[:, k] = r_trial
                sigma[k] = sigma_trial
                V_[k] = V_trial

                accept += 1

        # adaptive step
        if n % (NG * 10) == 0 and n > 0:

            acc_ratio = accept / (NG * 10)

            if acc_ratio < 0.4:
                step *= 0.9
            elif acc_ratio > 0.6:
                step *= 1.1

            accept = 0

            print(f"Step {n} | step_size = {step:.4f} | energy = {torch.sum(E_trial).item():.4f}")

    return rg, sigma, V_





