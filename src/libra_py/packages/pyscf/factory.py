# *********************************************************************************
# * Copyright (C) 2026 Jieyang Gu <jieyanggu792@gmail.com>
# * This file is distributed under the terms of the GNU General Public License
# * as published by the Free Software Foundation, either version 3 of
# * the License, or (at your option) any later version.
# * See the file LICENSE in the root directory of this distribution
# * or <http://www.gnu.org/licenses/>.
# *
# *********************************************************************************/
"""
Construction helpers for config-driven electronic-structure strategies.
"""

from __future__ import annotations

import importlib
import json
from pathlib import Path
from typing import Any

from libra_py.packages.pyscf.implementations import CASSCF, CISD
from libra_py.packages.pyscf.interfaces import ElectronicStructureStrategy


_STRATEGY_REGISTRY: dict[str, type[ElectronicStructureStrategy]] = {
    "casscf": CASSCF,
    "cisd": CISD,
}


def load_config(config_path: str | Path) -> dict[str, Any]:
    """Load a JSON configuration file."""
    path = Path(config_path)
    with path.open() as fh:
        return json.load(fh)


def build_strategy(es_cfg: dict[str, Any]) -> ElectronicStructureStrategy:
    """Instantiate an :class:`ElectronicStructureStrategy` from config.

    Supported config forms:

    1. Registry-based:

       {
         "kind": "casscf",
         "params": {...}
       }

    2. Dynamic import:

       {
         "module": "libra_py.packages.pyscf.implementations.casscf",
         "class": "CASSCF",
         "params": {...}
       }
    """
    params = dict(es_cfg.get("params", {}))

    kind = es_cfg.get("kind")
    if kind is not None:
        try:
            cls = _STRATEGY_REGISTRY[str(kind).lower()]
        except KeyError as exc:
            supported = ", ".join(sorted(_STRATEGY_REGISTRY))
            raise ValueError(
                f"Unknown strategy kind '{kind}'. Supported kinds: {supported}"
            ) from exc
    else:
        if "module" not in es_cfg or "class" not in es_cfg:
            raise ValueError(
                "Electronic-structure config must provide either 'kind' or "
                "both 'module' and 'class'."
            )
        module = importlib.import_module(es_cfg["module"])
        cls = getattr(module, es_cfg["class"])

    strategy = cls(**params)
    if not isinstance(strategy, ElectronicStructureStrategy):
        raise TypeError(
            f"{cls.__name__} does not implement ElectronicStructureStrategy"
        )
    return strategy


def build_strategy_factory(es_cfg: dict[str, Any]):
    """Return a zero-argument factory that builds fresh strategy instances."""

    def _factory() -> ElectronicStructureStrategy:
        return build_strategy(es_cfg)

    return _factory


__all__ = [
    "build_strategy",
    "build_strategy_factory",
    "load_config",
]
