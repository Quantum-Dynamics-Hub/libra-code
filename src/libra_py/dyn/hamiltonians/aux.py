from __future__ import annotations

from collections.abc import Mapping
from typing import Any


def as_result_mapping(result: Any) -> Mapping[str, Any]:
    """
    Return a mapping view of a model result.

    C++ `compute_diabatic` and `compute_adiabatic` accepted Python objects with
    named attributes. The Python dyn layer accepts either dictionaries or
    attribute-style objects so model code can stay lightweight.
    """

    if isinstance(result, Mapping):
        return result
    return _AttributeMapping(result)


def result_get(result: Any, name: str, default=None):
    """Read a field from a dict-like or attribute-style model result."""

    if isinstance(result, Mapping):
        return result.get(name, default)
    return getattr(result, name, default)


def require_allocated(storage: Any, *names: str) -> None:
    """Raise AttributeError if any requested TensorStorage field is None."""

    missing = [
        name
        for name in names
        if getattr(storage, name) is None
    ]
    if missing:
        joined = ", ".join(missing)
        raise AttributeError(f"TensorStorage field(s) not allocated: {joined}")


def check_matrix_shape(value: Any, shape: tuple[int, ...], name: str) -> None:
    """
    Validate the trailing matrix shape of a scalar or batched tensor.

    This is the NumPy/Python analogue of nHamiltonian::check_cmatrix.
    """

    actual = getattr(value, "shape", None)
    if actual is None:
        raise TypeError(f"{name} must be an array-like object with a shape")
    if tuple(actual[-len(shape):]) != shape:
        raise ValueError(
            f"{name} has trailing shape {tuple(actual[-len(shape):])}; "
            f"expected {shape}"
        )


def check_tensor_shape(value: Any, shape: tuple[int, ...], name: str) -> None:
    """Validate a full tensor shape."""

    actual = getattr(value, "shape", None)
    if actual is None:
        raise TypeError(f"{name} must be an array-like object with a shape")
    if tuple(actual) != shape:
        raise ValueError(f"{name} has shape {tuple(actual)}; expected {shape}")


def active_slice(storage: Any, traj: Any, name: str):
    """Return a writable trajectory/TBF slice from TensorStorage."""

    require_allocated(storage, name)
    return getattr(storage, name)[traj.id, traj.tbf_ids]


def write_active_slice(storage: Any, traj: Any, name: str, value: Any) -> None:
    """Assign one active trajectory/TBF slice in TensorStorage."""

    require_allocated(storage, name)
    getattr(storage, name)[traj.id, traj.tbf_ids] = value


class _AttributeMapping(Mapping):
    """Small adapter for attribute-style model result objects."""

    def __init__(self, obj: Any):
        self._obj = obj

    def __getitem__(self, key: str):
        value = getattr(self._obj, key)
        if value is None:
            raise KeyError(key)
        return value

    def __iter__(self):
        return (
            name
            for name in dir(self._obj)
            if not name.startswith("_") and getattr(self._obj, name) is not None
        )

    def __len__(self):
        return sum(1 for _ in self)
