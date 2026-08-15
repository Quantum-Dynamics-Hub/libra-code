#wigner.py is a Python module that generates quantum mechanical initial conditions for non-NBRA calculations
#using the Wigner phase-space distribution. 
#It reads vibrational mode data from computational chemistry calculations in mode_{}.xyz format.
# It requires mass of each atoms as {"N": 14.0, "C": 12.0, "H": 1.0}
# It requires Intial geometry in list format as, coords_eq_ang = [1.397389772,  -0.0008423226,  -0.0000282719,-1.397375322,   0.0008217057]
# Produces statistically correct initial geometries and momenta that account for both zero-point quantum energy and thermal effects.


import re
import numpy as np
from libra_py.units import kB
from libra_py import units

def read_modes_i_range(mode_start, mode_end, base_name="mode_{}.xyz"):
    """ 
    Reads vibrational mode data (displacements and frequencies) from a sequence of .xyz files.
    Input:
    mode_start, mode_end: The range of mode indices to read.
    base_name: The file naming pattern (defaults to "mode_{}.xyz").
    Process:
    It loops through the specified range of mode indices.
    For each file, it reads the number of atoms (natoms) from the first line.
    It extracts the vibrational frequency (in cm⁻¹) from the second line using a regular expression (re.search). 
    If not found, it defaults to 0.0.
    It reads the displacement vectors for each atom from the subsequent lines (ignoring the atom label in the first column) 
    and flattens them into a 1D array of length 3 * natoms.
    Output:
    natoms: The number of atoms.
    frequencies: A NumPy array of the vibrational frequencies.
    D: A displacement matrix where each column corresponds to a flattened normal mode vector.
    
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
    """
    Converts the Cartesian displacement vectors into mass-weighted, normalized eigenvectors,
    and then back into scaled Cartesian coordinates and momenta. 
    This ensures the modes are properly scaled according to the masses of the atoms involved.
    Input: 
    D: The displacement matrix.
    labels: A list of atomic symbols (e.g., ['C', 'H', 'N']).
    mass_map: A dictionary mapping atomic symbols to their masses in AMU.
    amu_to_au: Conversion factor from AMU to atomic units (electron masses).
    Process:
    Assigns atomic masses (in Atomic Mass Units, AMU) and converts them to atomic units (electron masses)
    using the factor amu2au = 1822.888.
    Creates an array of masses for each Cartesian coordinate (mass_cart) by repeating each atom's mass 3 times (for x, y, and z).
    Multiplies the Cartesian displacements (D) by the square root of the masses to get mass-weighted displacements (D_mw).
    Normalizes each mass-weighted column vector (D_mw_norm).
    Converts the normalized mass-weighted vectors back into a Cartesian representation for coordinates (D_cart) by dividing by 
    the square root of the masses, and for momenta (D_p) by multiplying by the square root of the masses.
    Outputs:
    D_cart: The appropriately scaled Cartesian transformation matrix for coordinates.
    D_p: The appropriately scaled Cartesian transformation matrix for momenta.
    mass_cart: The array of atomic masses in atomic units.
    
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
    """
    
    Generates a set of initial geometries and momenta (trajectories) sampled from a Wigner phase-space distribution, 
    which accounts for zero-point quantum energy and thermal effects for a harmonic oscillator.

    Inputs:
    q_eq: The equilibrium (ground state) Cartesian coordinates of the molecule.
    D_cart: The transformation matrix for coordinates from build_mass_weighted_eigenvectors.
    D_p: The transformation matrix for momenta from build_mass_weighted_eigenvectors.
    omega: The vibrational frequencies.
    temperature: The temperature of the system.
    ntraj: The number of trajectories (initial conditions) to generate.
    seed: A random seed for reproducibility.
    kB_value: Boltzmann constant.
    Process:
    Calculates the thermodynamic beta (1 / (kB * T)).
    For each vibrational mode, it calculates the standard deviations for position (sigma_q) and momentum (sigma_p) 
    according to the Wigner distribution for a quantum harmonic oscillator at a finite temperature. 
    (The np.tanh term accounts for thermal population of excited vibrational states).
    For each trajectory, it draws random normal mode displacements (dq_nm) and momenta (dp_nm) from normal distributions
    defined by sigma_q and sigma_p.
    It transforms these normal mode displacements and momenta back into Cartesian coordinates (q_eq + D_cart @ dq_nm) 
    and Cartesian momenta (D_p @ dp_nm).
    Outputs:
    ics: A list of dictionaries, where each dictionary contains the trajectory number ("traj"), 
    Cartesian coordinates ("q"), and Cartesian momenta ("p").
    
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
