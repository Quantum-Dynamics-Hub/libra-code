from typing import Any
import numpy as np

from .interfaces import ElectronicStructureStrategy, MolecularGeometry


class LibraAdapter:
    # Delegation only: adapter receives an already-constructed strategy.
    def __init__(self, interface: ElectronicStructureStrategy) -> None: ...

    # Metadata delegated from the strategy.
    @property
    def _atom_labels(self) -> list[str]: ...

    @property
    def _nstates(self) -> int: ...

    @property
    def _natoms(self) -> int: ...

    # Main Libra callback sequencing:
    # 1. extract trajectory index
    # 2. convert Libra q (Bohr) -> MolecularGeometry (Angstrom)
    # 3. call strategy.set_geom_and_run_hf(...)
    # 4. call strategy.compute_energy(root) for all roots
    # 5. call strategy.compute_gradient(root) for all roots
    # 6. call strategy.time_overlap_matrix(nstates)
    # 7. convert NumPy -> Libra CMATRIX/list[CMATRIX]
    def compute_model(
        self,
        q: "MATRIX",                      # Libra coordinates, flat Cartesian, Bohr
        params: dict[str, Any],
        full_id: Any,
    ) -> Any: ...

    # Sequencing helper: extract trajectory id for swarm dispatch.
    @staticmethod
    def _traj_index(full_id: Any) -> int: ...

    # Type + unit conversion:
    # Libra MATRIX column (Bohr) -> MolecularGeometry((natoms,3), Angstrom)
    def _libra_q_to_geometry(
        self,
        q: "MATRIX",
        traj_idx: int,
    ) -> MolecularGeometry: ...

    # Interface sequencing helper:
    # calls strategy.set_geom_and_run_hf(geom)
    def _update_strategy(
        self,
        geom: MolecularGeometry,
    ) -> None: ...

    # Interface sequencing helper:
    # calls strategy.compute_energy(root) for root = 0..nstates-1
    def _collect_energies(self) -> np.ndarray: ...        # shape: (nstates,), Hartree

    # Interface sequencing helper:
    # calls strategy.compute_gradient(root) for root = 0..nstates-1
    def _collect_gradients(self) -> np.ndarray: ...       # shape: (nstates, natoms, 3), Hartree/Bohr

    # Interface sequencing helper:
    # calls strategy.time_overlap_matrix(nstates)
    def _collect_time_overlap(self) -> np.ndarray: ...    # shape: (nstates, nstates), unitless

    # Type conversion:
    # NumPy energies -> Libra ham_adi
    def _build_ham_adi(
        self,
        energies: np.ndarray,
    ) -> "CMATRIX": ...

    # Type conversion:
    # NumPy gradients -> Libra d1ham_adi
    def _build_d1ham_adi(
        self,
        gradients: np.ndarray,
    ) -> list["CMATRIX"]: ...

    # Type conversion:
    # NumPy overlap -> Libra time_overlap_adi
    def _build_time_overlap_adi(
        self,
        time_overlap: np.ndarray,
    ) -> "CMATRIX": ...

    # FSSH packing helper:
    # zero explicit derivative couplings when using time-overlap NACs
    def _build_dc1_adi(self) -> list["CMATRIX"]: ...

    # FSSH packing helper:
    # identity fixed-geometry adiabatic overlap
    def _build_ovlp_adi(self) -> "CMATRIX": ...

    # Final assembly of the Python object Libra inspects.
    def _pack_result(
        self,
        energies: np.ndarray,
        gradients: np.ndarray,
        time_overlap: np.ndarray,
    ) -> Any: ...


class _LibraResult:
    ham_adi: "CMATRIX"
    d1ham_adi: list["CMATRIX"]
    time_overlap_adi: "CMATRIX"
    dc1_adi: list["CMATRIX"]
    ovlp_adi: "CMATRIX"


# Type conversion helper:
# Libra MATRIX column -> NumPy 1D vector
def _matrix_col_to_numpy(
    mat: "MATRIX",
    col: int,
) -> np.ndarray: ...


# Type conversion helper:
# NumPy 2D array -> Libra CMATRIX
def _cmatrix_from_numpy(
    arr: np.ndarray,
) -> "CMATRIX": ...


# Construction helper:
# Libra identity overlap matrix
def _identity_cmatrix(
    size: int,
) -> "CMATRIX": ...
