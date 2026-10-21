# Copyright (C) 2026 Alexey V. Akimov
# SPDX-License-Identifier: GPL-3.0-or-later

"""Gaussian wavepacket integrals and helpers."""

from . import integrals
from .integrals import (
    gwp_coupling,
    gwp_deriv,
    gwp_dipole,
    gwp_kinetic,
    gwp_kinetic_matrix,
    gwp_overlap,
    gwp_overlap_matrix,
    gwp_product_decomposition,
    gwp_value,
)

__all__ = [
    "integrals",
    "gwp_coupling",
    "gwp_deriv",
    "gwp_dipole",
    "gwp_kinetic",
    "gwp_kinetic_matrix",
    "gwp_overlap",
    "gwp_overlap_matrix",
    "gwp_product_decomposition",
    "gwp_value",
]
