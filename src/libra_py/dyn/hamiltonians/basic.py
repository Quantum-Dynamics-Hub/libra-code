from __future__ import annotations

from typing import Any


CORE_HAMILTONIAN_FIELDS = (
    "ovlp_dia",
    "ham_dia",
    "nac_dia",
    "hvib_dia",
    "ham_adi",
    "nac_adi",
    "hvib_adi",
    "basis_transform",
    "time_overlap_adi",
    "time_overlap_dia",
    "cum_phase_corr",
    "ordering_adi",
)

DERIVATIVE_LEVEL_1_FIELDS = (
    "dc1_dia",
    "d1ham_dia",
    "dc1_adi",
    "d1ham_adi",
)

DERIVATIVE_LEVEL_2_FIELDS = (
    "d2ham_dia",
    "d2ham_adi",
)


def init_hamiltonian_storage(storage: Any, der_lvl: int = 0) -> Any:
    """
    Allocate nHamiltonian-equivalent TensorStorage fields.

    TensorStorage already allocates the core Hamiltonian fields during
    construction. This function mirrors nHamiltonian::init_all by ensuring the
    optional derivative tensors are allocated when requested.
    """

    if storage.ham_adi is None or storage.ham_dia is None:
        storage.allocate_hamiltonian_vars()
    if der_lvl >= 1 and _missing_level_1_derivatives(storage):
        _allocate_derivatives_preserving(storage, der_lvl=der_lvl)
    elif der_lvl >= 2 and _missing_level_2_derivatives(storage):
        _allocate_derivatives_preserving(storage, der_lvl=2)
    return storage


def reset_hamiltonian_storage(storage: Any, der_lvl: int = 2) -> Any:
    """
    Zero Hamiltonian tensors that are currently allocated.

    This replaces C++ memory re-initialization with value reset semantics. It
    does not deallocate storage and does not touch nuclear/electronic amplitudes.
    """

    for name in _allocated_field_names(storage, der_lvl):
        getattr(storage, name)[...] = 0
    if storage.ordering_adi is not None:
        storage.ordering_adi[...] = _identity_ordering(storage)
    return storage


def copy_hamiltonian_content(dst: Any, src: Any, der_lvl: int = 2) -> Any:
    """
    Copy allocated Hamiltonian tensors from one TensorStorage object to another.

    Only fields allocated on both source and destination are copied, matching
    the conservative behavior of nHamiltonian::copy_level_content.
    """

    for name in _allocated_field_names(dst, der_lvl):
        source = getattr(src, name, None)
        target = getattr(dst, name, None)
        if source is not None and target is not None:
            target[...] = source
    dst.gs_kinetic_energy = getattr(src, "gs_kinetic_energy", dst.gs_kinetic_energy)
    return dst


def hamiltonian_memory_status(storage: Any, der_lvl: int = 2) -> dict[str, int]:
    """
    Return allocation status for Hamiltonian-related TensorStorage fields.

    Values follow the useful part of the old convention: 0 means unallocated and
    1 means allocated in TensorStorage. External pointer ownership does not
    exist in the Python storage model.
    """

    return {
        name: int(getattr(storage, name, None) is not None)
        for name in _field_names_for_level(der_lvl)
    }


def _allocated_field_names(storage: Any, der_lvl: int) -> tuple[str, ...]:
    return tuple(
        name
        for name in _field_names_for_level(der_lvl)
        if getattr(storage, name, None) is not None
    )


def _field_names_for_level(der_lvl: int) -> tuple[str, ...]:
    names = CORE_HAMILTONIAN_FIELDS
    if der_lvl >= 1:
        names += DERIVATIVE_LEVEL_1_FIELDS
    if der_lvl >= 2:
        names += DERIVATIVE_LEVEL_2_FIELDS
    return names


def _identity_ordering(storage: Any):
    ordering = storage.backend.zeros(
        storage.ordering_adi.shape,
        dtype=int,
    )
    ordering[...] = storage.backend.asarray(range(storage.nstates))
    return ordering


def _missing_level_1_derivatives(storage: Any) -> bool:
    return any(
        getattr(storage, name, None) is None
        for name in DERIVATIVE_LEVEL_1_FIELDS
    )


def _missing_level_2_derivatives(storage: Any) -> bool:
    return any(
        getattr(storage, name, None) is None
        for name in DERIVATIVE_LEVEL_2_FIELDS
    )


def _allocate_derivatives_preserving(storage: Any, der_lvl: int) -> None:
    names = DERIVATIVE_LEVEL_1_FIELDS + DERIVATIVE_LEVEL_2_FIELDS
    saved = {
        name: getattr(storage, name).copy()
        for name in names
        if getattr(storage, name, None) is not None
    }

    storage.allocate_hamiltonian_derivatives(der_lvl=der_lvl)

    for name, value in saved.items():
        target = getattr(storage, name, None)
        if target is not None and target.shape == value.shape:
            target[...] = value
