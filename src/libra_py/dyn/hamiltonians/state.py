
from __future__ import annotations
from dataclasses import dataclass, replace
from typing import Any, Optional, Dict


@dataclass(frozen=True)
class HamiltonianState:
    """
    Immutable snapshot of all Hamiltonian-related quantities
    for a batch of TBFs at a given time step.

    This replaces the data-storage role of nHamiltonian.
    """

    # =========================================================
    # Representation metadata
    # =========================================================
    rep: str = "adiabatic"  # "adiabatic" | "diabatic"

    basis_transform: Optional[Any] = None   # (N, S, S)
    ordering: Optional[Any] = None          # permutation info

    # =========================================================
    # Electronic Hamiltonians
    # =========================================================
    H_adi: Optional[Any] = None             # (N, S, S)
    H_dia: Optional[Any] = None             # (N, S, S)

    Hvib_adi: Optional[Any] = None          # (N, S, S)
    Hvib_dia: Optional[Any] = None          # (N, S, S)

    # =========================================================
    # Overlaps & couplings
    # =========================================================
    S_adi: Optional[Any] = None             # (N, S, S)
    S_dia: Optional[Any] = None             # (N, S, S)

    NAC_adi: Optional[Any] = None           # (N, S, S)
    NAC_dia: Optional[Any] = None           # (N, S, S)

    DC1_adi: Optional[Any] = None           # (N, D, S, S)
    DC1_dia: Optional[Any] = None           # (N, D, S, S)

    # =========================================================
    # Derivatives of Hamiltonian
    # =========================================================
    dH_adi: Optional[Any] = None            # (N, D, S, S)
    dH_dia: Optional[Any] = None            # (N, D, S, S)

    d2H_adi: Optional[Any] = None           # (N, D, D, S, S)
    d2H_dia: Optional[Any] = None           # (N, D, D, S, S)

    # =========================================================
    # Time propagation continuity
    # =========================================================
    time_overlap_adi: Optional[Any] = None  # (N, S, S)
    time_overlap_dia: Optional[Any] = None  # (N, S, S)

    phase_correction: Optional[Any] = None  # (N, S, S)

    # =========================================================
    # Cached observables (optional)
    # =========================================================
    energies: Optional[Any] = None          # (N, S) or (N,)
    forces: Optional[Any] = None            # (N, D)

    # =========================================================
    # Arbitrary metadata (SSY, LD flags, provenance, etc.)
    # =========================================================
    metadata: Optional[Dict[str, Any]] = None

    # =========================================================
    # Convenience properties
    # =========================================================
    @property
    def H(self):
        """Active Hamiltonian depending on representation."""
        return self.H_adi if self.rep == "adiabatic" else self.H_dia

    @property
    def Hvib(self):
        """Active vibronic Hamiltonian."""
        return self.Hvib_adi if self.rep == "adiabatic" else self.Hvib_dia

    @property
    def S(self):
        """Active overlap matrix."""
        return self.S_adi if self.rep == "adiabatic" else self.S_dia

    @property
    def NAC(self):
        return self.NAC_adi if self.rep == "adiabatic" else self.NAC_dia

    @property
    def DC1(self):
        return self.DC1_adi if self.rep == "adiabatic" else self.DC1_dia

    # =========================================================
    # Immutable update helper
    # =========================================================
    def replace(self, **kwargs) -> "HamiltonianState":
        """
        Functional-style update (keeps immutability semantics).
        """
        return replace(self, **kwargs)

    # =========================================================
    # Representation switching helpers
    # =========================================================
    def to_adiabatic(self) -> "HamiltonianState":
        return self.replace(rep="adiabatic")

    def to_diabatic(self) -> "HamiltonianState":
        return self.replace(rep="diabatic")

    # =========================================================
    # Metadata helper
    # =========================================================
    def set_meta(self, key: str, value: Any) -> "HamiltonianState":
        md = dict(self.metadata or {})
        md[key] = value
        return self.replace(metadata=md)



