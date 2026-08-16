# *********************************************************************************
# * Copyright (C) 2026 Somesh Chandra and Alexey V. Akimov
# *
# * This file is distributed under the terms of the GNU General Public License
# * as published by the Free Software Foundation, either version 3 of
# * the License, or (at your option) any later version.
# * See the file LICENSE in the root directory of this distribution
# * or <http://www.gnu.org/licenses/>.
# *
# *********************************************************************************/

# Samples initial conditions from the Wigner distribution.
# It reads vibrational mode data from computational chemistry calculations in mode_{}.xyz format.
# It requires mass of each atoms as {"N": 14.0, "C": 12.0, "H": 1.0}
# It requires Intial geometry in list format as, coords_eq_ang = [1.397389772,  -0.0008423226,  -0.0000282719,-1.397375322,   0.0008217057]
# Produces statistically correct initial geometries and momenta that account for both zero-point quantum energy and thermal effects.

import re
import numpy as np
from libra_py.units import kB
from libra_py import units

def read_modes_i_range(mode_start, mode_end, base_name="mode_{}.xyz"):
    """Read Cartesian normal modes from a numbered series of XYZ files.

    The first line of each file is the atom count.  The first numeric token on
    the second line is interpreted as a signed spectroscopic wavenumber in
    cm^-1.  Each following row contains an ignored atom label and three
    Cartesian displacement components.  No displacement-unit conversion is
    needed because the arbitrary scale of every mode cancels when it is
    normalized by :func:`build_mass_weighted_eigenvectors`.

    Parameters
    ----------
    mode_start, mode_end : int
        Inclusive range of mode indices.
    base_name : str, optional
        Filename pattern accepting the mode index via ``str.format``.

    Returns
    -------
    natoms : int
        Number of atoms common to all mode files.
    frequencies : numpy.ndarray, shape (nmodes,)
        Signed spectroscopic wavenumbers in cm^-1.
    D : numpy.ndarray, shape (3*natoms, nmodes)
        Cartesian mode vectors.  Coordinates are flattened in atom-major
        order: ``(x1, y1, z1, x2, y2, z2, ...)``.
    """
    displacements = []
    frequencies = []
    natoms = None

    for mode_idx in range(mode_start, mode_end + 1):
        filename = base_name.format(mode_idx)
        with open(filename, "r") as f:
            lines = f.readlines()

        if len(lines) < 3:
            raise ValueError(f"{filename} is too short to be a valid mode file")

        n = int(lines[0].strip())
        if natoms is None:
            natoms = n
        elif natoms != n:
            raise ValueError(
                f"Inconsistent natoms in {filename}: got {n}, expected {natoms}"
            )

        match = re.search(r'([-+]?\d+(?:\.\d+)?)', lines[1])
        if match is None:
            raise ValueError(f"Could not parse frequency from second line of {filename}")
        frequencies.append(float(match.group(1)))

        disp = np.zeros((natoms, 3), dtype=float)
        for i in range(natoms):
            parts = lines[i + 2].split()
            if len(parts) < 4:
                raise ValueError(
                    f"Bad format in {filename}, line {i+3}: {lines[i+2].strip()}"
                )
            disp[i, :] = [float(x) for x in parts[1:4]]

        displacements.append(disp.flatten())

    ndof = 3 * natoms
    D = np.zeros((ndof, len(displacements)), dtype=float)
    for k, vec in enumerate(displacements):
        D[:, k] = vec

    return natoms, np.array(frequencies, dtype=float), D


def build_mass_weighted_eigenvectors(D, labels, mass_map, amu_to_au=1822.888):
    """Build canonical normal-mode-to-Cartesian transformations.

    Let ``M`` be the diagonal Cartesian mass matrix and ``d_k`` a Cartesian
    mode vector.  Its normalized mass-weighted eigenvector is

    ``l_k = M^(1/2) d_k / sqrt(d_k.T M d_k)``.

    With the ``l_k`` collected as columns of ``L``, mass-weighted normal
    coordinates ``(Q, P)`` are transformed according to

    ``delta_q = M^(-1/2) L Q`` and ``p = M^(1/2) L P``.

    These coordinate and momentum transformations are canonically conjugate
    when the input modes are mutually orthogonal in the mass metric, i.e.
    ``L.T L = I``.  This routine normalizes individual columns but does not
    orthogonalize different modes.

    Parameters
    ----------
    D : numpy.ndarray, shape (3*N, nmodes)
        Cartesian displacement mode vectors, one per column.
    labels : sequence of str, length N
        Atomic labels in the same order as the Cartesian rows of ``D``.
    mass_map : mapping
        Atomic masses in unified atomic mass units (Da), keyed by label.
    amu_to_au : float, optional
        Conversion from Da to electron masses.

    Returns
    -------
    D_cart : numpy.ndarray
        ``M^(-1/2) L``, mapping normal coordinates to Cartesian displacement.
    D_p : numpy.ndarray
        ``M^(1/2) L``, mapping normal momenta to Cartesian momentum.
    mass_cart : numpy.ndarray, shape (3*N,)
        Cartesian masses in electron masses.

    Raises
    ------
    ValueError
        If a mode has near-zero mass-weighted norm.
    """
    masses = np.array([mass_map[a] * amu_to_au for a in labels], dtype=float)
    mass_cart = np.repeat(masses, 3)

    M_sqrt = np.sqrt(mass_cart)
    M_invsqrt = 1.0 / M_sqrt

    # Assumes D contains Cartesian normal-mode displacement vectors
    D_mw = D * M_sqrt[:, np.newaxis]

    norms = np.linalg.norm(D_mw, axis=0)
    if np.any(norms < 1e-15):
        raise ValueError("One or more mode vectors have near-zero norm.")

    D_mw_norm = D_mw / norms[np.newaxis, :]

    D_cart = M_invsqrt[:, np.newaxis] * D_mw_norm
    D_p = M_sqrt[:, np.newaxis] * D_mw_norm

    return D_cart, D_p, mass_cart


def generate_wigner_ics(q_eq, D_cart, D_p, omega, temperature, ntraj=1, seed=None, kB_value=None):
    """Sample thermal harmonic-oscillator Wigner initial conditions.

    Atomic units with ``hbar = 1`` are assumed.  For mode ``k``,

    ``H_k = (P_k^2 + omega_k^2 Q_k^2)/2``

    and its normalized thermal Wigner density is

    ``W_k = tanh(beta*omega_k/2)/pi``
    ``      * exp[-2*tanh(beta*omega_k/2)*H_k/omega_k]``,

    where ``beta = 1/(kB*T)``.  Therefore ``Q_k`` and ``P_k`` are independent,
    zero-mean Gaussian variables with

    ``Var(Q_k) = coth(beta*omega_k/2)/(2*omega_k)``,
    ``Var(P_k) = omega_k*coth(beta*omega_k/2)/2``.

    At ``T = 0`` the ``coth`` factor is one, which gives the ground-state
    Wigner distribution and mean energy ``omega_k/2``.  Exact zero modes are
    frozen because a free coordinate has no normalizable harmonic Wigner
    density; negative frequencies, representing imaginary modes, are rejected.

    Parameters
    ----------
    q_eq : array_like, shape (3*N,)
        Equilibrium Cartesian coordinates in bohr.
    D_cart, D_p : numpy.ndarray, shape (3*N, nmodes)
        Canonical coordinate and momentum transformations.
    omega : array_like, shape (nmodes,)
        Angular frequencies in atomic units (numerically hartree for
        ``hbar = 1``).
    temperature : float
        Temperature in kelvin; zero requests ground-state sampling.
    ntraj : int, optional
        Number of initial conditions.
    seed : int or None, optional
        Seed for NumPy's random-number generator.
    kB_value : float or None, optional
        Boltzmann constant in hartree/kelvin.  The default is
        :data:`libra_py.units.kB`.

    Returns
    -------
    list of dict
        Dictionaries containing ``traj``, ``q = q_eq + D_cart @ Q``, and
        ``p = D_p @ P``.

    Raises
    ------
    ValueError
        If a frequency is negative or no Boltzmann constant is available.
    """
    if kB_value is None:
        if kB is None:
            raise ValueError("kB_value must be provided if libra_py.units.kB is unavailable")
        kB_value = kB

    rng = np.random.default_rng(seed)

    n_modes = D_cart.shape[1]
    sigma_q = np.zeros(n_modes, dtype=float)
    sigma_p = np.zeros(n_modes, dtype=float)

    beta = 1.0 / (kB_value * temperature) if temperature > 0.0 else None

    for k in range(n_modes):
        w = omega[k]

        if w < -1e-12:
            raise ValueError(f"Negative/imaginary frequency detected at mode {k}: omega={w}")
        elif abs(w) <= 1e-12:
            continue

        coth = 1.0 / np.tanh(beta * w / 2.0) if temperature > 0.0 else 1.0
        sigma_q[k] = np.sqrt(coth / (2.0 * w))
        sigma_p[k] = np.sqrt(w * coth / 2.0)

    ics = []
    for traj in range(ntraj):
        dq_nm = rng.normal(0.0, sigma_q, size=n_modes)
        dp_nm = rng.normal(0.0, sigma_p, size=n_modes)

        ics.append({
            "traj": traj,
            "q": q_eq + D_cart @ dq_nm,
            "p": D_p @ dp_nm
        })

    return ics


def build_modes_from_hessian(hessian, masses, amu_to_au=1822.888,
                             zero_threshold=1.0e-12,
                             imaginary_threshold=1.0e-12):
    """Construct normal modes directly from a Cartesian Hessian.

    The Hessian is mass-weighted as ``M^(-1/2) H M^(-1/2)`` and diagonalized
    with :func:`numpy.linalg.eigh`.  All inputs and outputs use atomic units,
    except that ``masses`` are supplied in Da by default and converted using
    ``amu_to_au``.

    Parameters
    ----------
    hessian : array_like, shape (ndof, ndof)
        Cartesian Hessian in hartree/bohr**2.
    masses : array_like, shape (natoms,) or (ndof,)
        Atomic masses (one value per atom) or Cartesian masses (one value per
        degree of freedom), in Da.  For atomic masses, ``ndof`` must equal
        ``3*natoms``.
    amu_to_au : float, optional
        Conversion from Da to electron masses.  Set this to ``1.0`` when
        ``masses`` are already in atomic units.
    zero_threshold : float, optional
        Modes with ``abs(omega**2) <= zero_threshold`` are retained with zero
        frequency and consequently frozen by :func:`generate_wigner_ics`.
    imaginary_threshold : float, optional
        Negative eigenvalues below ``-imaginary_threshold`` are rejected.

    Returns
    -------
    omega : numpy.ndarray, shape (ndof,)
        Angular frequencies in atomic units, sorted in ascending eigenvalue
        order.  Translational/rotational zero modes are included as zeros.
    D_cart, D_p : numpy.ndarray, shape (ndof, ndof)
        Canonically conjugate transformations from mass-weighted normal-mode
        coordinates and momenta to Cartesian coordinates and momenta.
    mass_cart : numpy.ndarray, shape (ndof,)
        Cartesian masses in electron masses.

    Raises
    ------
    ValueError
        If dimensions or masses are invalid, the Hessian is not symmetric, or
        a genuine imaginary mode is present.
    """
    hessian = np.asarray(hessian, dtype=float)
    if hessian.ndim != 2 or hessian.shape[0] != hessian.shape[1]:
        raise ValueError("hessian must be a square two-dimensional array")
    if not np.all(np.isfinite(hessian)):
        raise ValueError("hessian must contain only finite values")
    if not np.allclose(hessian, hessian.T, rtol=1.0e-10, atol=1.0e-12):
        raise ValueError("hessian must be symmetric")

    ndof = hessian.shape[0]
    masses = np.asarray(masses, dtype=float).reshape(-1)
    if masses.size * 3 == ndof:
        mass_cart = np.repeat(masses, 3)
    elif masses.size == ndof:
        mass_cart = masses.copy()
    else:
        raise ValueError(
            f"masses must contain either {ndof // 3} atomic masses or "
            f"{ndof} Cartesian masses"
        )
    mass_cart *= amu_to_au
    if not np.all(np.isfinite(mass_cart)) or np.any(mass_cart <= 0.0):
        raise ValueError("all masses must be finite and positive")
    if zero_threshold < 0.0 or imaginary_threshold < 0.0:
        raise ValueError("mode thresholds must be non-negative")

    inv_sqrt_mass = 1.0 / np.sqrt(mass_cart)
    dynmat = (inv_sqrt_mass[:, None] * hessian) * inv_sqrt_mass[None, :]
    eigenvalues, eigenvectors = np.linalg.eigh(dynmat)

    imaginary = eigenvalues < -imaginary_threshold
    if np.any(imaginary):
        mode = int(np.flatnonzero(imaginary)[0])
        raise ValueError(
            "Negative Hessian eigenvalue (imaginary mode) detected at mode "
            f"{mode}: omega^2={eigenvalues[mode]}"
        )

    eigenvalues[np.abs(eigenvalues) <= zero_threshold] = 0.0
    # Tiny negative eigenvalues lying within the imaginary tolerance are
    # numerical noise and are treated as zero.
    eigenvalues[eigenvalues < 0.0] = 0.0
    omega = np.sqrt(eigenvalues)

    sqrt_mass = np.sqrt(mass_cart)
    D_cart = inv_sqrt_mass[:, None] * eigenvectors
    D_p = sqrt_mass[:, None] * eigenvectors
    return omega, D_cart, D_p, mass_cart


def generate_wigner_from_hessian(q_eq, hessian, masses, temperature,
                                  ntraj=1, seed=None, amu_to_au=1822.888,
                                  kB_value=None, zero_threshold=1.0e-12,
                                  imaginary_threshold=1.0e-12):
    """Generate Wigner initial conditions from a Cartesian Hessian.

    This convenience function combines :func:`build_modes_from_hessian` and
    :func:`generate_wigner_ics`.  The Hessian and equilibrium coordinates must
    be in atomic units; masses are in Da unless ``amu_to_au=1.0`` is used.

    Returns
    -------
    dict
        ``omega_au``, the normal-mode transformations, Cartesian masses, and
        the sampled initial conditions under the ``ics`` key.
    """
    omega, D_cart, D_p, mass_cart = build_modes_from_hessian(
        hessian, masses, amu_to_au=amu_to_au,
        zero_threshold=zero_threshold,
        imaginary_threshold=imaginary_threshold
    )
    q_eq = np.asarray(q_eq, dtype=float)
    if q_eq.ndim != 1 or q_eq.size != mass_cart.size:
        raise ValueError("q_eq must be a one-dimensional vector matching the Hessian")
    ics = generate_wigner_ics(
        q_eq, D_cart, D_p, omega, temperature, ntraj=ntraj, seed=seed,
        kB_value=kB_value
    )
    return {
        "omega_au": omega,
        "D_cart": D_cart,
        "D_p": D_p,
        "mass_cart_au": mass_cart,
        "ics": ics
    }



def prepare_wigner_from_modes(
    labels,
    q_eq,
    mode_start,
    mode_end,
    mode_file_pattern="mode_{}.xyz",
    mass_map=None,
    temperature=0.0,
    ntraj=1,
    seed=None,
    freq_to_au=None,
    amu_to_au=1822.888,
    kB_value=None
):
    """Read normal modes and generate Wigner initial conditions.

    The input wavenumber ``nu_bar`` is converted with

    ``omega_au = nu_bar_cm^-1 * inv_cm2Ha``.

    This has no missing ``2*pi``: a spectroscopic wavenumber represents the
    energy ``h*c*nu_bar = hbar*omega``, and energy equals angular frequency in
    atomic units.  Equilibrium coordinates must already be in bohr.  Masses
    are supplied in Da and converted to electron masses.

    Parameters
    ----------
    labels : sequence of str
        Atomic labels in mode-file order.
    q_eq : array_like, shape (3*N,)
        Flattened equilibrium geometry in bohr.
    mode_start, mode_end : int
        Inclusive mode-file index range.
    mode_file_pattern : str, optional
        Filename pattern accepted by ``str.format``.
    mass_map : mapping
        Atomic masses in Da, keyed by label.  Required.
    temperature : float, optional
        Sampling temperature in kelvin.
    ntraj : int, optional
        Number of initial conditions.
    seed : int or None, optional
        Random seed.
    freq_to_au : float or None, optional
        Conversion from cm^-1 to hartree; defaults to ``units.inv_cm2Ha``.
    amu_to_au : float, optional
        Conversion from Da to electron masses.
    kB_value : float or None, optional
        Boltzmann constant in hartree/kelvin.

    Returns
    -------
    dict
        Parsed frequencies, transformations, Cartesian masses, dimensions,
        and sampled initial conditions.

    Notes
    -----
    The caller must select vibrational modes consistently and normally exclude
    translations and rotations.  Negative/imaginary modes are not sampled.
    """
    if mass_map is None:
        raise ValueError("mass_map must be provided")

    if freq_to_au is None:
        if units is None:
            raise ValueError("freq_to_au must be provided if libra_py.units is unavailable")
        freq_to_au = units.inv_cm2Ha

    nat = len(labels)
    ndof = 3 * nat

    q_eq = np.array(q_eq, dtype=float)
    if q_eq.shape[0] != ndof:
        raise ValueError(f"q_eq has length {q_eq.shape[0]}, expected {ndof}")

    natoms_modes, freqs_parsed, D = read_modes_i_range(
        mode_start, mode_end, mode_file_pattern
    )

    if natoms_modes != nat:
        raise ValueError(
            f"Mode files contain {natoms_modes} atoms, but labels define {nat} atoms."
        )

    omega = freqs_parsed * freq_to_au

    D_cart, D_p, mass_cart = build_mass_weighted_eigenvectors(
        D, labels, mass_map, amu_to_au=amu_to_au
    )

    ics = generate_wigner_ics(
        q_eq=q_eq,
        D_cart=D_cart,
        D_p=D_p,
        omega=omega,
        temperature=temperature,
        ntraj=ntraj,
        seed=seed,
        kB_value=kB_value
    )

    return {
        "natoms": nat,
        "ndof": ndof,
        "freqs_cm": freqs_parsed,
        "omega_au": omega,
        "D_cart": D_cart,
        "D_p": D_p,
        "mass_cart_au": mass_cart,
        "ics": ics
    }
