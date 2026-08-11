"""PySCF implementation backends."""

from __future__ import annotations

from .casscf import CASSCF

__all__ = ["CASSCF", "CISD"]

try:
    from .cisd import CISD
except ModuleNotFoundError:  # pragma: no cover - optional backend
    CISD = None
