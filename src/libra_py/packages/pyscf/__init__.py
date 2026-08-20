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
    "CISD_States",
    "CASSCF",
    "TDDFT",
    "TDDFT_States",
]


def __getattr__(name: str):
    if name in {"CASSCF", "CISD", "CISD_States", "TDDFT", "TDDFT_States"}:
        from .implementations import CASSCF, CISD, CISD_States, TDDFT, TDDFT_States

        return {
            "CASSCF": CASSCF,
            "CISD": CISD,
            "CISD_States": CISD_States,
            "TDDFT": TDDFT,
            "TDDFT_States": TDDFT_States,
        }[name]
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
