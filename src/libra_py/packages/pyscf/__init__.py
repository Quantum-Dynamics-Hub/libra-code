"""Libra PySCF package interface.

Keep this module lightweight so importing ``libra_py.packages.pyscf.interfaces``
does not eagerly pull in concrete backends or adapter code.
"""

from __future__ import annotations

from .interfaces import ElectronicStructureStrategy, MolecularGeometry

__all__ = [
    "ElectronicStructureStrategy",
    "MolecularGeometry",
    "CISD",
    "CASSCF",
    "build_strategy",
    "build_strategy_factory",
    "load_config",
    "LibraESAdapter",
    "LibraNAMDAdapter",
    "MultiTrajNAMDAdapter",
    "NAMDRunner",
]


def __getattr__(name: str):
    if name in {"CISD", "CASSCF"}:
        from .implementations import CASSCF, CISD

        return {"CISD": CISD, "CASSCF": CASSCF}[name]
    if name in {"build_strategy", "build_strategy_factory", "load_config"}:
        from .factory import build_strategy, build_strategy_factory, load_config

        return {
            "build_strategy": build_strategy,
            "build_strategy_factory": build_strategy_factory,
            "load_config": load_config,
        }[name]
    if name in {"LibraESAdapter", "LibraNAMDAdapter", "MultiTrajNAMDAdapter", "NAMDRunner"}:
        from .adapter import LibraESAdapter, LibraNAMDAdapter, MultiTrajNAMDAdapter, NAMDRunner

        return {
            "LibraESAdapter": LibraESAdapter,
            "LibraNAMDAdapter": LibraNAMDAdapter,
            "MultiTrajNAMDAdapter": MultiTrajNAMDAdapter,
            "NAMDRunner": NAMDRunner,
        }[name]
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
