# *********************************************************************************
# * Copyright (C) 2018-2025 Alexey V. Akimov
# *
# * This file is distributed under the terms of the GNU General Public License
# * as published by the Free Software Foundation, either version 3 of
# * the License, or (at your option) any later version.
# * See the file LICENSE in the root directory of this distribution
# * or <http://www.gnu.org/licenses/>.
# *
# *********************************************************************************/

"""
.. module:: normal_modes
   :platform: Unix, Windows
   :synopsis: ========= Theory ==============
       Convention on the normal modes transformations:
       U^T * H * U = W^2       (1)
       x = U * q               (2a)
       so q = U^T * x          (2b)
       x = M^(1/2) * (R-R_ave) (3)

       Here:
       H - mass-scaled Hessian (dynamical matrix): H_ij = (d^2E/dx_i * dx_j )/sqrt(m_i * m_j)
       U, q - normal modes (q are columns of U)
       x - mass-scaled and shifted coordinates
       R - Cartesian coordinates

       Note: to relate to the definitions here (as of 11/7/2018), use: U = D^T
       https://www.pci.uni-heidelberg.de//tc/usr/mctdh/doc/vcham/latex/nmode_coo.pdf

       Also useful:
       https://nanohub.org/courses/FATM/01a/asset/256

.. moduleauthor:: Alexey V. Akimov


"""

import os
import math
import sys
if sys.platform == "cygwin":
    from cyglibra_core import *
elif sys.platform == "linux" or sys.platform == "linux2":
    from liblibra_core import *

from . import units
from . import md_align


import numpy as np



def covariance_matrix(X, M, flag=1):
    """
    Compute the mass-weighted covariance matrix of trajectory data.

    The covariance matrix is defined as

        K_ij = 1/2 ⟨ sqrt(m_i m_j) x_i x_j ⟩

    or, if centered (flag = 1),

        K_ij = 1/2 ⟨ sqrt(m_i m_j)
                     (x_i - ⟨x_i⟩)(x_j - ⟨x_j⟩) ⟩

    where ⟨·⟩ denotes averaging over trajectory steps.

    Parameters
    ----------
    X : np.ndarray, shape (ndof, nsteps)
        Trajectory data (e.g., coordinates, velocities, or accelerations).

    M : array_like, shape (ndof,)
        Masses associated with each degree of freedom.

    flag : {0, 1}, optional
        Controls centering:
        - 0 : use raw data (no mean subtraction)
        - 1 : use fluctuations around the mean (default)

    Returns
    -------
    K : np.ndarray, shape (ndof, ndof)
        Mass-weighted covariance matrix averaged over the trajectory.

    Notes
    -----
    - Mass weighting is performed using sqrt(m_i).
    - The prefactor 1/2 is included to match the legacy implementation.
    - Averaging is done over the trajectory dimension (axis=1).
    - This form is commonly used in harmonic analysis,
      PCA, and normal-mode–like decompositions of MD data.

    Examples
    --------
    Compute the mass-weighted covariance of Cartesian coordinates:

    >>> R, E = read_trajectory("traj.xyz")
    >>> M = names_to_masses(E, PT)
    >>> K = covariance_matrix(R, M, flag=1)

    Diagonalize the covariance matrix to obtain principal modes:

    >>> eigvals, eigvecs = np.linalg.eigh(K)

    Use uncentered data (no subtraction of the mean):

    >>> K_raw = covariance_matrix(R, M, flag=0)
    """

    X = np.asarray(X, dtype=float)
    M = np.asarray(M, dtype=float)

    if X.ndim != 2:
        raise ValueError("X must be a 2D array of shape (ndof, nsteps)")
    if M.ndim != 1 or M.shape[0] != X.shape[0]:
        raise ValueError("M must be a 1D array of length ndof")

    # Center data if requested
    if flag == 1:
        X = X - X.mean(axis=1, keepdims=True)
    elif flag != 0:
        raise ValueError("flag must be 0 (no centering) or 1 (centering)")

    # Mass weighting
    sM = np.sqrt(M)[:, None]   # shape (ndof, 1)
    Xmw = sM * X               # mass-weighted data

    # Covariance with prefactor 1/2
    K = 0.5 * (Xmw @ Xmw.T) / Xmw.shape[1]

    return K



def normal_modes_from_eigenvectors(eta, M):
    """
    Convert mass-weighted eigenvectors to Cartesian normal modes.

    The transformation is

        Δx = M^{-1/2} η

    where:
        - η are mass-weighted eigenvectors
        - M is the Cartesian mass vector

    Parameters
    ----------
    eta : np.ndarray, shape (ndof, nmodes) or (ndof,)
        Mass-weighted eigenvectors, typically obtained from
        diagonalizing a mass-weighted covariance or Hessian matrix.

    M : array_like, shape (ndof,)
        Cartesian mass vector (each atomic mass repeated for x, y, z).

    Returns
    -------
    dx : np.ndarray, shape (ndof, nmodes) or (ndof,)
        Cartesian normal modes.

    Notes
    -----
    - No normalization is enforced; normalization conventions depend
      on whether the eigenvectors came from a covariance or Hessian.
    - This operation is equivalent to left-multiplying by a diagonal
      matrix M^{-1/2}.
    - For Hessian-based normal modes, this recovers physical
      displacement directions in Cartesian space.
    - For Hessian-based normal modes, frequencies are obtained from eigenvalues of
      M^{−1/2} H M^{−1/2}
    - For covariance-based modes, eigenvalues correspond to variances, not frequencies
    - Translational and rotational modes should be projected out before diagonalization

    Examples
    --------
    From a covariance matrix:

    >>> R, E = read_trajectory("traj.xyz")
    >>> M = names_to_masses(E, PT)
    >>> K = covariance_matrix(R, M)
    >>> w, eta = np.linalg.eigh(K)
    >>> dx = normal_modes_from_eigenvectors(eta, M)

    Single mode extraction:

    >>> dx0 = normal_modes_from_eigenvectors(eta[:, 0], M)

    Reshape modes to (nat, 3):

    >>> nat = len(E)
    >>> modes_xyz = dx.reshape(3*nat, -1).T.reshape(-1, nat, 3)
    """

    eta = np.asarray(eta, dtype=float)
    M = np.asarray(M, dtype=float)

    if M.ndim != 1:
        raise ValueError("M must be a 1D array of length ndof")

    ndof = M.shape[0]

    if eta.ndim == 1:
        if eta.shape[0] != ndof:
            raise ValueError("eta and M must have the same length")
        return eta / np.sqrt(M)

    elif eta.ndim == 2:
        if eta.shape[0] != ndof:
            raise ValueError("eta must have shape (ndof, nmodes)")
        return eta / np.sqrt(M)[:, None]

    else:
        raise ValueError("eta must be a 1D or 2D array")


def compute_thermal_normal_modes(
    R,
    M,
    method="QHA",
    V=None,
    A=None,
    F=None,
    temperature=None,
    remove_rigid=False,
):
    """
    Unified covariance-based normal-mode analysis at finite temperature.

    This function implements several related but conceptually distinct
    normal-mode definitions based on covariance matrices sampled from MD.

    ----------------------------------------------------------------------
    METHODS IMPLEMENTED
    ----------------------------------------------------------------------

    PCA (Essential Dynamics)
    -----------------------
    Matrix:
        C^r_ij = < δr_i δr_j >

    Interpretation:
        Eigenvectors → collective motions with largest positional variance
        Eigenvalues  → variances

    Use when:
        - Interested in dominant conformational changes
        - No physical frequencies required

    ----------------------------------------------------------------------

    QHA (Quasi-Harmonic Analysis)
    -----------------------------
    Matrix:
        K^r_ij = 1/2 < sqrt(m_i m_j) δr_i δr_j >

    Frequencies:
        ω_i^2 = k_B T / λ_i

    Interpretation:
        Effective free-energy curvature at finite temperature

    Use when:
        - Near-equilibrium fluctuations
        - Moderate anharmonicity
        - Want thermally averaged vibrational spectrum

    Reference:
        Brooks, B. R.; Janezic, D.; Karplus, M. Harmonic Analysis of Large Systems. I. Methodology. 
        J Comput Chem 1995, 16 (12), 1522–1542. https://doi.org/10.1002/jcc.540161209.

    ----------------------------------------------------------------------

    Strachan-V (Velocity covariance)
    --------------------------------
    Matrix:
        K^v_ij = 1/2 < sqrt(m_i m_j) v_i v_j >

    Frequencies:
        (2πν_i)^2 = λ_i^v / λ_i^r

    Interpretation:
        Eigenvalues represent kinetic energy per mode

    Use when:
        - Well-converged velocities
        - Avoid explicit temperature dependence

    Reference:
        Strachan, A. Normal Modes and Frequencies from Covariances in Molecular Dynamics or Monte Carlo 
        Simulations. The Journal of Chemical Physics 2004, 120 (1), 1–4. https://doi.org/10.1063/1.1635364.

    ----------------------------------------------------------------------

    Strachan-A (Acceleration covariance)
    ------------------------------------
    Matrix:
        K^a_ij = 1/2 < sqrt(m_i m_j) a_i a_j >

    Frequencies:
        (2πν_i)^4 = λ_i^a / λ_i^r

    Interpretation:
        Force-derived vibrational modes

    Use when:
        - Reliable forces
        - Strong anharmonicity

    Reference:
        Strachan, A. Normal Modes and Frequencies from Covariances in Molecular Dynamics or Monte Carlo                                                                        Simulations. The Journal of Chemical Physics 2004, 120 (1), 1–4. https://doi.org/10.1063/1.1635364.

    ----------------------------------------------------------------------

    Pereverzev Thermal Hessian
    --------------------------
    Matrix:
        H_ij(T) = β < F_i F_j > ,   β = 1 / (k_B T)

    Mass-weighted form:
        H̃_ij = β < sqrt(m_i) F_i sqrt(m_j) F_j >

    Eigenproblem:
        H̃ η = ω^2 η

    Interpretation:
        Finite-temperature generalization of the Hessian

    Use when:
        - Strong anharmonicity
        - Liquids / fluxional clusters
        - Instantaneous Hessians are noisy or unavailable

    Reference:
        - Pereverzev, A.; Sewell, T. D. Obtaining the Hessian from the Force Covariance Matrix: Application to 
          Crystalline Explosives PETN and RDX. The Journal of Chemical Physics 2015, 142 (13), 134110. https://doi.org/10.1063/1.4916614.

        - Martinez, M.; Gaigeot, M.-P.; Borgis, D.; Vuilleumier, R. Extracting Effective Normal Modes from Equilibrium Dynamics 
          at Finite Temperature. The Journal of Chemical Physics 2006, 125 (14), 144106. https://doi.org/10.1063/1.2346678.

    ----------------------------------------------------------------------

    PARAMETERS
    ----------------------------------------------------------------------
    R : ndarray, shape (3N, nsteps)
        Cartesian coordinates (Bohr).

    M : ndarray, shape (3N,)
        Mass vector (atomic units).

    method : str
        One of:
        "PCA", "QHA", "Strachan-V", "Strachan-A", "Pereverzev"

    V : ndarray, optional
        Velocities, shape (3N, nsteps).

    A : ndarray, optional
        Accelerations, shape (3N, nsteps).

    F : ndarray, optional
        Forces, shape (3N, nsteps), required for Pereverzev.

    temperature : float, optional
        Temperature in Kelvin (required for QHA and Pereverzev).

    remove_rigid : bool, default False
        Project out translations and rotations using md_align.project_out_rigid.

    ----------------------------------------------------------------------
    RETURNS
    ----------------------------------------------------------------------
    results : dict
        Keys:
        - eigvals       : eigenvalues (a.u.)
        - eigvecs       : mass-weighted eigenvectors
        - modes_cart    : Cartesian normal modes
        - frequencies   : frequencies (cm^-1), if defined
    """

    ndof, nsteps = R.shape
    sqrtM = np.sqrt(M)

    # --------------------------------------------------
    # Optional rigid-body projection
    # --------------------------------------------------
    if remove_rigid:
        R_proj = np.zeros_like(R)
        for t in range(nsteps):
            R_proj[:, t] = md_align.project_out_rigid(R[:, t], R[:, t], M)
        R = R_proj

        if V is not None:
            Vp = np.zeros_like(V)
            for t in range(nsteps):
                Vp[:, t] = md_align.project_out_rigid(V[:, t], R[:, t], M)
            V = Vp

        if A is not None:
            Ap = np.zeros_like(A)
            for t in range(nsteps):
                Ap[:, t] = md_align.project_out_rigid(A[:, t], R[:, t], M)
            A = Ap

        if F is not None:
            Fp = np.zeros_like(F)
            for t in range(nsteps):
                Fp[:, t] = md_align.project_out_rigid(F[:, t], R[:, t], M)
            F = Fp


    # --------------------------------------------------
    # Covariance helper
    # --------------------------------------------------
    def covariance(X, center):
        if center:
            X = X - X.mean(axis=1, keepdims=True)
        return (X @ X.T) / X.shape[1]

    # --------------------------------------------------
    # Positional covariance (always needed)
    # --------------------------------------------------
    Crr = covariance(R, center=True)
    Kx = 0.5 * (sqrtM[:, None] * sqrtM[None, :]) * Crr

    eigvals_x, eigvecs_x = np.linalg.eigh(Kx)
    idx = np.argsort(eigvals_x)
    eigvals_x = eigvals_x[idx]
    eigvecs_x = eigvecs_x[:, idx]

    results = {
        "eigvals": eigvals_x,
        "eigvecs": eigvecs_x,
        "modes_cart": eigvecs_x / sqrtM[:, None],
    }

    # --------------------------------------------------
    # Method-specific frequencies
    # --------------------------------------------------
    if method == "PCA":
        return results

    if method == "QHA":
        if temperature is None:
            raise ValueError("Temperature required for QHA")
        omega = np.sqrt(units.kB * temperature / eigvals_x)
        results["frequencies"] = omega / units.inv_cm2Ha
        return results

    if method == "Strachan-V":
        Kv = 0.5 * (sqrtM[:, None] * sqrtM[None, :]) * covariance(V, center=False)
        eigvals_v, _ = np.linalg.eigh(Kv)
        eigvals_v = eigvals_v[idx]
        omega = np.sqrt(eigvals_v / eigvals_x)
        results["frequencies"] = omega / units.inv_cm2Ha
        return results

    if method == "Strachan-A":
        Ka = 0.5 * (sqrtM[:, None] * sqrtM[None, :]) * covariance(A, center=False)
        eigvals_a, _ = np.linalg.eigh(Ka)
        eigvals_a = eigvals_a[idx]
        omega = (eigvals_a / eigvals_x) ** 0.25
        results["frequencies"] = omega / units.inv_cm2Ha
        return results

    if method == "Pereverzev":
        if F is None or temperature is None:
            raise ValueError("Forces and temperature required for Pereverzev")
        beta = 1.0 / (units.kB * temperature)
        inv_sqrtM = 1.0 / sqrtM

        Cff = covariance(F, center=False)
        H_eff = beta * (inv_sqrtM[:, None] * inv_sqrtM[None, :]) * Cff
        #A = F/M[:, None]
        #H_eff = beta * (sqrtM[:, None] * sqrtM[None, :]) * covariance(A, center=False)
        eigvals, eigvecs = np.linalg.eigh(H_eff)
        idx = np.argsort(eigvals)
        results["eigvals"] = eigvals[idx]
        results["eigvecs"] = eigvecs[:, idx]
        results["modes_cart"] = eigvecs[:, idx] / sqrtM[:, None]
        results["frequencies"] = np.sqrt(eigvals[idx]) / units.inv_cm2Ha
        return results

    raise ValueError(f"Unknown method: {method}")



def write_nmd(
    filename,
    R_eq,
    modes,
    atom_names,
    masses=None,
    eigenvalues=None,
    title="Normal modes",
):
    """
    Write a ProDy/VMD-compatible NMD (Normal Mode Data) file.

    This implementation follows the official ProDy NMD specification:
    http://www.bahargroup.org/prody/manual/reference/dynamics/nmdfile.html

    Parameters
    ----------
    filename : str
        Output .nmd file.

    R_eq : ndarray, shape (3N,)
        Equilibrium Cartesian coordinates (e.g., Å).

    modes : ndarray, shape (3N, nmodes)
        Normal mode eigenvectors. If masses are provided, these are
        assumed to be mass-weighted eigenvectors η.

    atom_names : list[str], length N
        Atom or element names (written as `atomnames`).

    masses : ndarray, shape (3N,), optional
        Cartesian mass vector. If provided, modes are converted via:
            Δx = M^{-1/2} η

    eigenvalues : ndarray, shape (nmodes,), optional
        Eigenvalues associated with the modes.
        Used to define a scaling factor per mode.

    title : str
        Descriptive title written using the `name` field.

    Examples
    --------
    Compute the mass-weighted covariance of Cartesian coordinates:

    >>> write_nmd(
    ...     filename="modes.nmd",
    ...     R_eq=R[:, 0],          # equilibrium geometry
    ...     modes=eigvecs_x,       # mass-weighted eigenvectors η
    ...     atom_names=E,          # e.g. ["Ti", "O", "O", ...]
    ...     masses=M,              # length 3N
    ...     eigenvalues=eigvals_x,
    ...     title="Strachan covariance normal modes",
    ... )

    """

    R_eq = np.asarray(R_eq, dtype=float)
    modes = np.asarray(modes, dtype=float)

    N = len(atom_names)
    nmodes = modes.shape[1]

    if R_eq.shape != (3 * N,):
        raise ValueError("R_eq must have shape (3N,)")

    if modes.shape[0] != 3 * N:
        raise ValueError("modes must have shape (3N, nmodes)")

    # -------------------------------------------------
    # Convert mass-weighted modes to Cartesian Δx
    # -------------------------------------------------
    if masses is not None:
        masses = np.asarray(masses, dtype=float)
        if masses.shape != (3 * N,):
            raise ValueError("masses must have shape (3N,)")
        modes_cart = modes / np.sqrt(masses[:, None])
    else:
        modes_cart = modes.copy()

    # -------------------------------------------------
    # Normalize modes (recommended for NMWiz)
    # -------------------------------------------------
    for k in range(nmodes):
        norm = np.linalg.norm(modes_cart[:, k])
        if norm > 0:
            modes_cart[:, k] /= norm

    # -------------------------------------------------
    # Scaling factors (optional, per ProDy spec)
    # -------------------------------------------------
    if eigenvalues is not None:
        eigenvalues = np.asarray(eigenvalues, dtype=float)
        if eigenvalues.shape != (nmodes,):
            raise ValueError("eigenvalues must have shape (nmodes,)")
        scales = [
            np.sqrt(1.0 / ev) if ev > 0.0 else 1.0
            for ev in eigenvalues
        ]
    else:
        scales = [1.0] * nmodes

    #--------------------------------------------------
    # Map names to ids
    #--------------------------------------------------
    name_to_id = {}
    resids = [name_to_id.setdefault(n, len(name_to_id)) for n in atom_names]


    # -------------------------------------------------
    # Write NMD file
    # -------------------------------------------------
    with open(filename, "w") as f:

        # Title
        f.write(f"name {title}\n")

        # Atom names (single line!)
        f.write("atomnames " + " ".join(atom_names) + "\n")

        # Residue names (single line!)
        f.write("resnames " + " ".join(atom_names) + "\n")

        # Residue ids (single line!)
        f.write("resids " + " ".join(str(x) for x in resids) + "\n")

        # Coordinates (single line!)
        coord_str = " ".join(f"{x:.8f}" for x in R_eq)
        f.write(f"coordinates {coord_str}\n")

        # Modes (each on one line)
        for k in range(nmodes):
            mode_str = " ".join(f"{v:.8e}" for v in modes_cart[:, k])
            f.write(f"mode {k+1} {scales[k]:.8e} {mode_str}\n")

    print(f"Wrote NMD file: {filename}")




def visualize_modes(E, R, eta, M, frequencies, params):
    """
    Visualize selected normal modes by generating XYZ trajectories.

    The visualization is performed by constructing a synthetic trajectory
    corresponding to sinusoidal motion along selected normal modes.

    Args:
        E (list of str, length N):
            Atom names.

        R (ndarray, shape (3N,)):
            Reference Cartesian geometry in atomic units.

        eta (ndarray, shape (3N, 3N)):
            Mass-weighted normal-mode eigenvectors.

        M (ndarray, shape (3N,)):
            Masses of Cartesian degrees of freedom (in a.u.).

        frequencies (ndarray, shape (3N,)):
            Mode frequencies in cm^{-1}.

        params (dict):
            Visualization parameters:
                - "scale" (float): mode amplification factor
                - "print_modes" (list[int]): indices of modes to visualize
                - "prefix" (str): output filename prefix
                - "nperiods" (int): number of oscillation periods
                - "nsteps" (int): number of frames per trajectory

    Returns:
        None
        Writes XYZ files with animated normal-mode motion (in Angstrom).

    Notes:
        Cartesian displacements are obtained from mass-weighted modes via:
            Δx = M^{-1/2} η

        Frequencies are converted internally as:
            ω (a.u.) = ν (cm^{-1}) × units.inv_cm2Ha
    """

    R = np.asarray(R)
    eta = np.asarray(eta)
    M = np.asarray(M)
    frequencies = np.asarray(frequencies)

    ndof = R.size
    nat = ndof // 3

    sqrtM = np.sqrt(M)
    scale = params["scale"]

    # Convert frequencies from cm^{-1} to atomic units
    omega_au = frequencies * units.inv_cm2Ha

    for mode in params["print_modes"]:

        filename = f"{params['prefix']}_mode{mode}.xyz"

        omega = omega_au[mode]
        if omega > 0.0:
            dt = 2.0 * math.pi * params["nperiods"] / (omega * params["nsteps"])
        else:
            dt = 0.0

        # Cartesian displacement for this mode
        dx = eta[:, mode] / sqrtM

        with open(filename, "w") as f:
            for t in range(params["nsteps"]):

                phase = math.sin(omega * t * dt)
                Rt = R + scale * phase * dx

                f.write(f"{nat}\n")
                f.write(f"mode {mode}, step {t}\n")

                Rt_xyz = Rt.reshape(nat, 3) / units.Angst

                for at in range(nat):
                    f.write(
                        f"{E[at]:2s} "
                        f"{Rt_xyz[at,0]:10.5f} "
                        f"{Rt_xyz[at,1]:10.5f} "
                        f"{Rt_xyz[at,2]:10.5f}\n"
                    )


def compute_dynmat(R, D, M, E, params):
    """

    Computes and visualizes (as the trajectories) normal modes
    using the dynamic matrix:   D_ij = [1/sqrt(m_i * m_j)]  d^2E/dR_i dR_j

    Here, H_ij = d^2E/dR_i dR_j is the Hessian

    Args:
        R ( MATRIX(ndof x nsteps-1) ): coordinates of all DOFs for all mid-timesteps
        D ( MATRIX(ndof x nsteps-1) ): the dynamic matrix
        M ( MATRIX(ndof x 1) ): masses of all DOFs
        E ( list of ndof/3 strings): atom names (elements) of all atoms
        params ( dictionary ): parameters controlling the computations, including the
            visualization (see the visualize_modes(E, R, U, w, params) description).
            Contains keyword-value pairs:

            * **params["verbosity"]** ( int ): level to control verbosity
            * **params["visualize"]** ( int ): flag to control whether we want to produce additional files (with normal modes)
               - 0 - not to
               - 1 - do it

    Returns:
        tuple: (w_a, w_inv_cm, U_a), where:

            * w_a ( MATRIX(ndof,1) ): frequencies - Hessian eigenvalues
            * w_inv_cm ( MATRIX(ndof,1) ): frequencies from - Hessian eigenvalues, in cm^-1 units
            * U_a ( MATRIX(ndof,ndof) ): Hessian eigenvectors

    Note:
        All quantities are in atomic units

    """

    verbosity = params["verbosity"]

    if verbosity > 0:
        print("========= Normal modes calculations using the provided dynamical matrix =============================")

    ndof = M.num_of_rows
    nat = ndof / 3

    w_a = MATRIX(ndof, ndof)
    U_a = MATRIX(ndof, ndof)

    solve_eigen_nosort(D, w_a, U_a, 0)

    if verbosity > 1:
        print("Dynamic matrix:")
        D.show_matrix()
        print("Its eigenvalues:")
        w_a.show_matrix()
        print("Its eigenvectors:")
        U_a.show_matrix()

    w = MATRIX(ndof, 1)
    for dof in range(0, ndof):
        if w_a.get(dof, dof) > 0.0:
            w.set(dof, 0, math.sqrt(w_a.get(dof, dof)))
    w_inv_cm = w / units.inv_cm2Ha
    if verbosity > 0:
        print("Frequencies (cm^-1)")
        w_inv_cm.show_matrix()

    if params["visualize"] > 0:
        if verbosity > 0:
            print("Visualizing modes based on dynamic matrix\n")
        visualize_modes(E, R, U_a, M, w, params)

    if verbosity > 0:
        print("========= Done with the Normal modes calculations =============================")

    return w, w_inv_cm, U_a


def get_xyz(E, R, M, U, mode):
    """

    This function returns a string in the xyz format with X, Y, Z and UX, UY, UZ
    where X,Y,Z are the coordinates, UX, UY, UZ - vectors coming from those coordinates - e.g. normal modes

    Args:
        E ( list of ndof/3 strings ): atom names (elements) of all atoms
        R ( MATRIX(ndof x nsteps-1) ): coordinates of all DOFs for all mid-timesteps
        M ( MATRIX(ndof x 1) ): masses of all DOFs
        U ( MATRIX(ndof x ndof) ): a matrix containing normal mode vectors
        mode ( int ): index of the normal mode that we want to visualize

    Returns:
        string: A string representing an xyz file

    """

    natoms = len(E)
    res = "%3i\nComment\n" % (natoms)

    for i in range(0, natoms):
        x, y, z = R.get(3 * i, 0), R.get(3 * i + 1, 0), R.get(3 * i + 2, 0)
        ux = U.get(3 * i + 0, mode) / math.sqrt(M.get(3 * i + 0))
        uy = U.get(3 * i + 1, mode) / math.sqrt(M.get(3 * i + 1))
        uz = U.get(3 * i + 2, mode) / math.sqrt(M.get(3 * i + 2))

        res = res + "%s  %5.3f  %5.3f  %5.3f  %5.3f  %5.3f  %5.3f\n" % (E[i], x, y, z, ux, uy, uz)

    return res


def get_xyz2(E, R, U, mode):
    """

    This function returns a string in the xyz format with X, Y, Z and UX, UY, UZ
    where X,Y,Z are the coordinates, UX, UY, UZ - vectors coming from those coordinates - e.g. normal modes

    Args:
        E ( list of ndof/3 string ): atom names (elements) of all atoms
        R ( MATRIX(ndof x nsteps-1) ): coordinates of all DOFs for all mid-timesteps
        U ( MATRIX(ndof x ndof) ): a matrix containing normal mode vectors
        mode ( int ): index of the normal mode that we want to visualize

    Returns:
        string: A string representing an xyz file

    """

    natoms = len(E)
    res = "%3i\nComment\n" % (natoms)

    for i in range(0, natoms):
        x, y, z = R.get(3 * i, 0), R.get(3 * i + 1, 0), R.get(3 * i + 2, 0)
        ux, uy, uz = U.get(3 * i + 0, mode), U.get(3 * i + 1, mode), U.get(3 * i + 2, mode)

        res = res + "%s  %5.3f  %5.3f  %5.3f  %5.3f  %5.3f  %5.3f\n" % (E[i], x, y, z, ux, uy, uz)

    return res
