"""Libra PySCF package interface.

Keep this module lightweight so importing ``libra_py.packages.pyscf.interfaces``
does not eagerly pull in concrete backends or adapter code.
"""

from __future__ import annotations

from libra_py.packages.pyscf.interfaces import ElectronicStructureStrategy, MolecularGeometry

__all__ = [
    "interfaces",
    "implementations",
    "methods",
    "ElectronicStructureStrategy",
    "MolecularGeometry",
    "CISD",
    "CASSCF",
    "TDDFT",
    "TDDFT_States",
]


def __getattr__(name: str):
    if name in {"CASSCF", "TDDFT", "TDDFT_States"}:
        from .implementations import CASSCF, TDDFT, TDDFT_States

        return { "CASSCF": CASSCF, "TDDFT": TDDFT, "TDDFT_States": TDDFT_States }[name]
    if name == "CISD":
        raise AttributeError( "CISD backend is not available yet; use CASSCF or TDDFT" )
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
