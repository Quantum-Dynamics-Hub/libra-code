"""
Tutorial-style PySCF/Libra example using the formal adapter architecture.

This is an alternative to the inline ``pyscf_compute_adi`` function in
``tutorial.ipynb``. The tutorial builds Libra-native return objects manually
inside the callback. Here, the responsibilities are split explicitly:

* ``CASSCF`` computes backend-neutral energies / gradients / time overlaps
* ``LibraESAdapter`` converts Libra ``MATRIX`` input into ``MolecularGeometry``
* ``LibraESAdapter`` packs the PySCF results back into Libra-native
  ``CMATRIX`` structures

The example includes a direct callback suitable for Libra TSH/NAMD drivers.

Run:

    python src/libra_py/packages/pyscf/examples/tutorial_adapter_casscf.py
"""

from __future__ import annotations

from typing import Any

from liblibra_core import MATRIX

import libra_py.units as units
from libra_py.packages.pyscf.adapter import LibraESAdapter
from libra_py.packages.pyscf.implementations.casscf import CASSCF


ATOM_LABELS = ["Li", "H"]
NSTATES = 2


def build_strategy() -> CASSCF:
    """Construct the backend-neutral PySCF strategy."""
    return CASSCF(
        norbcas=2,
        nelecas=2,
        nroots=NSTATES,
        basis="sto-3g",
        charge=0,
    )


def build_compute_model() -> Any:
    """Return a Libra callback with the standard ``compute_model`` signature."""
    adapter = LibraESAdapter(build_strategy(), ATOM_LABELS, NSTATES)

    def compute_model(q: Any, params: dict[str, Any], full_id: Any) -> Any:
        return adapter.compute_model(q, params, full_id)

    return compute_model


def make_demo_q(bond_angstrom: float) -> MATRIX:
    """Construct a single-trajectory Libra coordinate matrix in Bohr."""
    q = MATRIX(3 * len(ATOM_LABELS), 1)
    bond_bohr = bond_angstrom * units.Angst

    coords = [
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        bond_bohr,
    ]
    for i, value in enumerate(coords):
        q.set(i, 0, value)
    return q


def run_callback_demo() -> None:
    """Mimic the tutorial's direct callback usage with Libra-native ``q``."""
    compute_model = build_compute_model()
    q = make_demo_q(1.60)
    params: dict[str, Any] = {}

    first = compute_model(q, params, [0, 0])
    second = compute_model(q, params, [0, 0])

    print("Direct callback demo")
    print("ham_adi diagonal (Hartree):")
    for i in range(NSTATES):
        print(f"  state {i}: {second.ham_adi.get(i, i).real:.10f}")
    print("time-overlap matrix from second call:")
    second.time_overlap_adi.real().show_matrix()
    print("first gradient component dE/dq_0:")
    for i in range(NSTATES):
        print(f"  state {i}: {second.d1ham_adi[0].get(i, i).real:.10f}")


def main() -> None:
    run_callback_demo()


if __name__ == "__main__":
    main()
