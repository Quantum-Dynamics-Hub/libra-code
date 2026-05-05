"""Libra PySCF package interface.

Keep this module lightweight so importing ``libra_py.packages.pyscf.interfaces``
does not eagerly pull in concrete backends or adapter code.
"""

from __future__ import annotations

from .interfaces import ElectronicStructureStrategy, MolecularGeometry

__all__ = [
    "interfaces",
    "implementations",
    "methods",
    "ElectronicStructureStrategy",
    "MolecularGeometry",
    "CISD",
    "CASSCF",
]


def __getattr__(name: str):
    if name in {"CISD", "CASSCF"}:
        from .implementations import CASSCF, CISD

        return {"CISD": CISD, "CASSCF": CASSCF}[name]
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
