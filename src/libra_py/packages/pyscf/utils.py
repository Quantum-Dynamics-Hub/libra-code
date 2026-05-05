"""
Shared helpers for PySCF-backed electronic-structure strategies.
"""

from __future__ import annotations

from typing import Any, Optional

import numpy as np
from pyscf import gto, scf

from .interfaces import MolecularGeometry, ElectronicStructureStrategy


def geometry_to_pyscf_atom_spec(geom: MolecularGeometry) -> str:
    """Convert MolecularGeometry into the atom-spec string expected by PySCF."""
    return ";".join(
        f"{label} {coord[0]} {coord[1]} {coord[2]}"
        for label, coord in zip(geom.atom_labels, geom.coords_angstrom)
    )


def run_rhf_for_geometry(
    geom: MolecularGeometry,
    *,
    basis: str,
    unit: str,
    charge: int,
    prev_mol: Optional[Any] = None,
    overlap_integral: str = "int1e_ovlp",
    spin: int = 0,
    verbose: int = 0,
    dm0: Optional[np.ndarray] = None,
) -> tuple[Any, Any, Optional[np.ndarray]]:
    """Build a PySCF molecule, optionally compute cross-geometry AO overlap, and run RHF."""
    mol = gto.M(
        atom=geometry_to_pyscf_atom_spec(geom),
        basis=basis,
        unit=unit,
        charge=charge,
        spin=spin,
    )

    ao_overlap = None
    if prev_mol is not None:
        ao_overlap = gto.intor_cross(overlap_integral, prev_mol, mol)

    mf = scf.RHF(mol)
    if dm0 is not None:
        mf.run(dm0=dm0, verbose=verbose)
    else:
        mf.run(verbose=verbose)
    return mol, mf, ao_overlap


class PySCFBasedStrategy(ElectronicStructureStrategy):
    def __init__(self, mol=None, nroots=1, basis="sto-3g", unit="Angstrom", charge=0, overlap_integral="int1e_ovlp"):
        super().__init__(mol=mol, nroots=nroots, basis=basis, unit=unit, charge=charge)
        self._solver = None
        self._prev_mol = None
        self._prev_mf = None
        self._prev_solver = None
        self._overlap_integral = overlap_integral

    def save_cache(self) -> None:
        self._prev_mol = self._mol
        self._prev_mf = self._mf
        self._prev_solver = self._solver

    def run_hf(self) -> None:
        self._solver = None
        self._ao_overlap = None

        dm0 = None
        if self._prev_mf is not None:
            try:
                dm0 = self._prev_mf.make_rdm1()
            except Exception:
                dm0 = None

        self._mol, self._mf, self._ao_overlap = run_rhf_for_geometry(
            self._geom,
            basis=self._basis,
            unit=self._unit,
            charge=self._charge,
            prev_mol=self._prev_mol,
            overlap_integral=self._overlap_integral,
            spin=0,
            verbose=0,
            dm0=dm0,
        )
