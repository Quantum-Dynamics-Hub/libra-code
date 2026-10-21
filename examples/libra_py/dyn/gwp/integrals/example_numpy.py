"""Construct two 2D packets explicitly and evaluate their Gaussian integrals.

Run from the repository root:
    PYTHONPATH=src python examples/libra_py/dyn/gwp/integrals/example_numpy.py

All quantities use atomic units (hbar = 1). Each vector has one entry per
spatial degree of freedom (DOF). No C++ extension is required.
"""
import numpy as np

from libra_py.dyn.gwp import integrals


def main():
    """Evaluate integrals, parameter gradients, and two-packet basis matrices."""
    # Packet 1: q is position in bohr; p is momentum in atomic units.
    q1 = np.array([0.0, 0.5], dtype=float)
    p1 = np.array([0.2, -0.1], dtype=float)
    # gamma is a phase PER COORDINATE. Its sum is the packet's total phase.
    gamma1 = np.array([0.0, 0.0], dtype=float)
    # alp is positive, in bohr^-2. Larger alp means a narrower packet;
    # the probability density has variance 1/(4*alp) along each coordinate.
    alp1 = np.array([1.0, 0.8], dtype=float)

    # Packet 2 has different centers, momenta, phases, and widths.
    q2 = np.array([0.3, -0.2], dtype=float)
    p2 = np.array([-0.1, 0.4], dtype=float)
    gamma2 = np.array([0.1, 0.0], dtype=float)
    alp2 = np.array([0.7, 1.2], dtype=float)

    # One mass per DOF, measured in electron masses; the API takes 1/mass.
    mass = np.array([1.0, 2.0], dtype=float)
    inv_mass = 1.0 / mass

    # The first packet is conjugated (the bra); the second is the ket.
    # Overlap and kinetic energy are complex scalars after integrating both DOFs.
    overlap = integrals.gwp_overlap(q1, p1, gamma1, alp1, q2, p2, gamma2, alp2)
    kinetic = integrals.gwp_kinetic(q1, p1, gamma1, alp1, q2, p2, gamma2, alp2, inv_mass)
    # Position and spatial derivative matrix elements have one component per DOF.
    dipole = integrals.gwp_dipole(q1, p1, gamma1, alp1, q2, p2, gamma2, alp2)
    coupling = integrals.gwp_coupling(q1, p1, gamma1, alp1, q2, p2, gamma2, alp2)
    print("Overlap <G1|G2>:", overlap)
    print("Kinetic energy <G1|T|G2>:", kinetic)
    print("Position <G1|r|G2>:", dipole)
    print("Spatial derivative <G1|gradient|G2>:", coupling)

    # r is ONE point in 2D, not two separate 1D sample points.
    r = np.array([0.1, 0.2], dtype=float)
    # This four-argument overload assumes zero phase, matching gamma1 above.
    value = integrals.gwp_value(r, q1, p1, alp1)
    gradient = integrals.gwp_deriv(r, q1, p1, alp1)
    # Omitting p requests |G1|^2, not a zero-momentum wavefunction.
    density = integrals.gwp_value(r, q1, alp1)
    print("G1(r):", value)
    print("Gradient of G1 at r:", gradient)
    print("Probability density at r:", density)
    np.testing.assert_allclose(density, abs(value)**2)

    # The pointwise product is scale * G_product. Neither factor is conjugated.
    scale, q_product, p_product, gamma_product, alp_product = integrals.gwp_product_decomposition(
        q1, p1, gamma1, alp1, q2, p2, gamma2, alp2
    )
    print("Product scale:", scale)
    print("Product center:", q_product)
    print("Product momentum:", p_product)
    print("Product phase per DOF:", gamma_product)
    print("Product width:", alp_product)
    # A twelve-argument overlap conjugates BOTH first packets. Here G3 = G1.
    triple = integrals.gwp_overlap(
        q1, p1, gamma1, alp1, q2, p2, gamma2, alp2, q1, p1, gamma1, alp1
    )
    print("Integral conj(G1*G2)*G1:", triple)

    # Differentiate the real squared overlap with respect to packet 1's
    # centers and widths. These parameter derivatives differ from dG/dr above.
    # NumPy has no autograd: use centered differences, changing one DOF at a
    # time while holding all other parameters fixed. No optimization is done.
    overlap_probability = abs(overlap)**2
    step = 1e-6
    q_gradient = np.zeros_like(q1)
    alp_gradient = np.zeros_like(alp1)
    for dof in range(q1.size):
        q_plus = q1.copy()
        q_minus = q1.copy()
        q_plus[dof] += step
        q_minus[dof] -= step
        overlap_plus = integrals.gwp_overlap(q_plus, p1, gamma1, alp1, q2, p2, gamma2, alp2)
        overlap_minus = integrals.gwp_overlap(q_minus, p1, gamma1, alp1, q2, p2, gamma2, alp2)
        q_gradient[dof] = (abs(overlap_plus)**2 - abs(overlap_minus)**2) / (2*step)

        alp_plus = alp1.copy()
        alp_minus = alp1.copy()
        alp_plus[dof] += step
        alp_minus[dof] -= step
        overlap_plus = integrals.gwp_overlap(q1, p1, gamma1, alp_plus, q2, p2, gamma2, alp2)
        overlap_minus = integrals.gwp_overlap(q1, p1, gamma1, alp_minus, q2, p2, gamma2, alp2)
        alp_gradient[dof] = (abs(overlap_plus)**2 - abs(overlap_minus)**2) / (2*step)

    print("Squared overlap:", overlap_probability)
    print("d|S12|^2/dq1:", q_gradient)
    print("d|S12|^2/dalp1:", alp_gradient)
    assert np.isfinite(q_gradient).all()
    assert np.isfinite(alp_gradient).all()

    # Build each basis matrix explicitly from the named packet parameters.
    # axis=1 puts packets in COLUMNS: shape (ndof=2, ntraj=2).
    q_basis = np.stack([q1, q2], axis=1)
    p_basis = np.stack([p1, p2], axis=1)
    gamma_basis = np.stack([gamma1, gamma2], axis=1)
    alp_basis = np.stack([alp1, alp2], axis=1)
    # Reuse the same basis on both sides to get 2x2 Hermitian matrices.
    overlap_matrix = integrals.gwp_overlap_matrix(
        q_basis, p_basis, gamma_basis, alp_basis,
        q_basis, p_basis, gamma_basis, alp_basis
    )
    kinetic_matrix = integrals.gwp_kinetic_matrix(
        q_basis, p_basis, gamma_basis, alp_basis,
        q_basis, p_basis, gamma_basis, alp_basis, inv_mass
    )
    print("Basis overlap matrix:\n", overlap_matrix)
    print("Basis kinetic matrix:\n", kinetic_matrix)

    # Electronic labels are one per PACKET, not per DOF. Orthogonal states
    # remove overlaps between packets with different labels (off-diagonal here).
    states = np.array([0, 1], dtype=int)
    electronic_overlap = integrals.gwp_overlap_matrix(
        q_basis, p_basis, gamma_basis, alp_basis, states,
        q_basis, p_basis, gamma_basis, alp_basis, states
    )
    print("Overlap including orthogonal electronic states:\n", electronic_overlap)
    np.testing.assert_allclose(np.diag(overlap_matrix), 1.0)
    np.testing.assert_allclose(overlap_matrix, overlap_matrix.conj().T)
    np.testing.assert_allclose(kinetic_matrix, kinetic_matrix.conj().T)


if __name__ == "__main__":
    main()
