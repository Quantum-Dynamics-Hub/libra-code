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
Compatibility wrapper for the shared Libra NAMD adapter.

The real implementation lives in :mod:`libra_py.namd_adapter`.  This module
keeps the older PySCF-local import path stable while ensuring the PySCF package
reuses the same generic adapter used by other backends.
"""

from __future__ import annotations

from libra_py.namd_adapter import (
    LibraNAMDAdapter,
    MultiTrajNAMDAdapter,
    NAMDRunner,
)
from libra_py.packages.pyscf.interfaces import ElectronicStructureStrategy


class LibraESAdapter(LibraNAMDAdapter):
    """Backward-compatible single-trajectory adapter.

    Older PySCF examples instantiated ``LibraESAdapter(strategy, atom_labels,
    nstates)``. The shared adapter derives the state count from the strategy, so
    ``nstates`` is now validated and otherwise ignored.
    """

    def __init__(
        self,
        strategy: ElectronicStructureStrategy,
        atom_labels: list[str],
        nstates: int | None = None,
        compute_gradients: bool = True,
        use_nac_vectors: bool = False,
    ) -> None:
        if nstates is not None and int(nstates) != int(strategy.nstates):
            raise ValueError(
                f"Requested nstates={nstates}, but strategy reports "
                f"{strategy.nstates} states."
            )
        super().__init__(
            strategy=strategy,
            atom_labels=atom_labels,
            compute_gradients=compute_gradients,
            use_nac_vectors=use_nac_vectors,
        )


__all__ = [
    "LibraESAdapter",
    "LibraNAMDAdapter",
    "MultiTrajNAMDAdapter",
    "NAMDRunner",
]
