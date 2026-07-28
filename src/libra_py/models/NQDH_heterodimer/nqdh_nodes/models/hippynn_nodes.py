"""Small hippy-nn graph node helpers used by project scripts."""

from __future__ import annotations

from collections.abc import Sequence

import torch
from hippynn.graphs import IdxType
from hippynn.graphs.nodes.base import SingleNode


class CatModule(torch.nn.Module):
    """Concatenate graph tensors along one dimension."""

    def __init__(self, dim: int = 1) -> None:
        super().__init__()
        self.dim = dim

    def forward(self, *values: torch.Tensor) -> torch.Tensor:
        return torch.cat(values, dim=self.dim)


class MoleculeCatNode(SingleNode):
    """Concatenate molecule-indexed tensors, preserving molecule index state."""

    _index_state = IdxType.Molecules

    def __init__(self, name: str, parents: Sequence, dim: int = 1) -> None:
        main_outputs = tuple(parent.main_output for parent in parents)
        super().__init__(name, main_outputs, module=CatModule(dim=dim))
