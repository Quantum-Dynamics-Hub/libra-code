from __future__ import annotations

from typing import Any, Callable

from .aux import as_result_mapping, check_matrix_shape, write_active_slice
from .basic import init_hamiltonian_storage


DIABATIC_RESULT_FIELDS = {
    "ham_dia": "ham_dia",
    "H_dia": "ham_dia",
    "ovlp_dia": "ovlp_dia",
    "S_dia": "ovlp_dia",
    "nac_dia": "nac_dia",
    "NAC_dia": "nac_dia",
    "hvib_dia": "hvib_dia",
    "Hvib_dia": "hvib_dia",
    "dc1_dia": "dc1_dia",
    "DC1_dia": "dc1_dia",
    "d1ham_dia": "d1ham_dia",
    "dH_dia": "d1ham_dia",
    "d2ham_dia": "d2ham_dia",
    "d2H_dia": "d2ham_dia",
    "time_overlap_dia": "time_overlap_dia",
}


def compute_diabatic(
    storage: Any,
    traj: Any,
    model_fn: Callable,
    params: Any = None,
    der_lvl: int = 1,
) -> Any:
    """
    Evaluate diabatic model data and write it into TensorStorage.

    This is the storage-backed translation of nHamiltonian::compute_diabatic.
    The model function receives keyword arguments `R`, `P`, `storage`, `traj`,
    and `params`, and may return a dict or attribute object containing any of
    the fields listed in DIABATIC_RESULT_FIELDS.
    """

    init_hamiltonian_storage(storage, der_lvl=der_lvl)
    idx = traj.tbf_ids
    result = as_result_mapping(
        model_fn(
            R=storage.q[traj.id, idx],
            P=storage.p[traj.id, idx],
            storage=storage,
            traj=traj,
            params=params,
        )
    )

    for source, target in DIABATIC_RESULT_FIELDS.items():
        value = result.get(source)
        if value is None:
            continue
        _validate_diabatic_field(storage, value, source, target)
        write_active_slice(storage, traj, target, value)

    if "hvib_dia" not in result and "Hvib_dia" not in result:
        build_diabatic_hvib(storage, traj)

    return storage


def build_diabatic_hvib(storage: Any, traj: Any) -> Any:
    """Build `hvib_dia = ham_dia - i * nac_dia` for active TBF slices."""

    idx = traj.tbf_ids
    if storage.ham_dia is None or storage.hvib_dia is None:
        return storage
    value = storage.ham_dia[traj.id, idx]
    if storage.nac_dia is not None:
        value = value - 1j * storage.nac_dia[traj.id, idx]
    storage.hvib_dia[traj.id, idx] = value
    return storage


def _validate_diabatic_field(storage: Any, value: Any, source: str, target: str) -> None:
    ns = storage.ndia
    nd = storage.nnucl
    if target in ("dc1_dia", "d1ham_dia"):
        check_matrix_shape(value, (nd, ns, ns), source)
    elif target == "d2ham_dia":
        check_matrix_shape(value, (nd, nd, ns, ns), source)
    else:
        check_matrix_shape(value, (ns, ns), source)
