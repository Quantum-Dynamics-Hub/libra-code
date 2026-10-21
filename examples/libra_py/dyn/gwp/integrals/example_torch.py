"""Construct two 2D tensor packets and evaluate the same integrals as NumPy.

Run from the repository root:
    PYTHONPATH=src python examples/libra_py/dyn/gwp/integrals/example_torch.py

All quantities use atomic units (hbar = 1). Each vector has one entry per
spatial degree of freedom (DOF). Use --device cuda for a supported GPU.
Neither SymPy nor the C++ extension is required.
"""
import argparse

import torch

from libra_py.dyn.gwp import integrals


def main(device="cpu"):
    """Evaluate integrals, parameter gradients, and two-packet basis matrices."""
    # Float64 inputs give complex128 integrals. All tensors share a device.
    dtype = torch.float64

    # Packet 1: q is position in bohr; p is momentum in atomic units.
    # Request gradients for q1 and alp1; other parameters stay fixed.
    q1 = torch.tensor([0.0, 0.5], dtype=dtype, device=device, requires_grad=True)
    p1 = torch.tensor([0.2, -0.1], dtype=dtype, device=device)
    # gamma is a phase PER COORDINATE. Its sum is the packet's total phase.
    gamma1 = torch.tensor([0.0, 0.0], dtype=dtype, device=device)
    # alp is positive, in bohr^-2. Larger alp means a narrower packet;
    # the probability density has variance 1/(4*alp) along each coordinate.
    alp1 = torch.tensor([1.0, 0.8], dtype=dtype, device=device, requires_grad=True)

    # Packet 2 has different centers, momenta, phases, and widths.
    q2 = torch.tensor([0.3, -0.2], dtype=dtype, device=device)
    p2 = torch.tensor([-0.1, 0.4], dtype=dtype, device=device)
    gamma2 = torch.tensor([0.1, 0.0], dtype=dtype, device=device)
    alp2 = torch.tensor([0.7, 1.2], dtype=dtype, device=device)

    # One mass per DOF, measured in electron masses; the API takes 1/mass.
    mass = torch.tensor([1.0, 2.0], dtype=dtype, device=device)
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
    r = torch.tensor([0.1, 0.2], dtype=dtype, device=device)
    # This four-argument overload assumes zero phase, matching gamma1 above.
    value = integrals.gwp_value(r, q1, p1, alp1)
    gradient = integrals.gwp_deriv(r, q1, p1, alp1)
    # Omitting p requests |G1|^2, not a zero-momentum wavefunction.
    density = integrals.gwp_value(r, q1, alp1)
    print("G1(r):", value)
    print("Gradient of G1 at r:", gradient)
    print("Probability density at r:", density)
    assert torch.allclose(density, abs(value)**2, rtol=1e-12, atol=1e-13)

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
    # PyTorch uses autograd; the NumPy example uses centered differences.
    # A real scalar objective lets backward() use an implicit unit seed,
    # without importing SymPy. This computes gradients, not an optimizer step.
    overlap_probability = overlap.abs().square()
    overlap_probability.backward()
    q_gradient = q1.grad
    alp_gradient = alp1.grad
    print("Squared overlap:", overlap_probability)
    print("d|S12|^2/dq1:", q_gradient)
    print("d|S12|^2/dalp1:", alp_gradient)
    assert torch.isfinite(q_gradient).all()
    assert torch.isfinite(alp_gradient).all()

    # Build each basis matrix explicitly from the named packet parameters.
    # dim=1 puts packets in COLUMNS: shape (ndof=2, ntraj=2).
    # Stacking preserves the connection to q1/alp1 in the autograd graph.
    q_basis = torch.stack([q1, q2], dim=1)
    p_basis = torch.stack([p1, p2], dim=1)
    gamma_basis = torch.stack([gamma1, gamma2], dim=1)
    alp_basis = torch.stack([alp1, alp2], dim=1)
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
    states = torch.tensor([0, 1], dtype=torch.int64, device=device)
    electronic_overlap = integrals.gwp_overlap_matrix(
        q_basis, p_basis, gamma_basis, alp_basis, states,
        q_basis, p_basis, gamma_basis, alp_basis, states
    )
    print("Overlap including orthogonal electronic states:\n", electronic_overlap)
    assert torch.allclose(overlap_matrix.diag(), torch.ones(2, dtype=overlap_matrix.dtype, device=device))
    assert torch.allclose(overlap_matrix, overlap_matrix.mH, rtol=1e-12, atol=1e-13)
    assert torch.allclose(kinetic_matrix, kinetic_matrix.mH, rtol=1e-12, atol=1e-13)
    assert overlap_matrix.device == q1.device
    assert overlap_matrix.requires_grad
    assert kinetic_matrix.requires_grad


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--device", default="cpu", help="PyTorch device, e.g. cpu or cuda")
    main(parser.parse_args().device)
