from __future__ import annotations

from typing import Any, Callable, Iterable, Optional

from .aux import as_result_mapping


class HamiltonianEngine:
    """
    Stateless Hamiltonian construction and transformation engine.

    The engine writes Hamiltonian quantities directly into TensorStorage. This
    keeps TensorStorage as the single owner of dynamical tensors and avoids a
    parallel HamiltonianState container drifting out of sync.
    """

    _RESULT_TO_STORAGE = {
        "H_adi": "ham_adi",
        "H_dia": "ham_dia",
        "Hvib_adi": "hvib_adi",
        "Hvib_dia": "hvib_dia",
        "NAC_adi": "nac_adi",
        "NAC_dia": "nac_dia",
        "DC1_adi": "dc1_adi",
        "DC1_dia": "dc1_dia",
        "dH_adi": "d1ham_adi",
        "dH_dia": "d1ham_dia",
        "d2H_adi": "d2ham_adi",
        "d2H_dia": "d2ham_dia",
        "S_dia": "ovlp_dia",
        "time_overlap_adi": "time_overlap_adi",
        "time_overlap_dia": "time_overlap_dia",
        "basis_transform": "basis_transform",
        "phase_correction": "cum_phase_corr",
        "ordering": "ordering_adi",
    }

    _REP_FIELDS = {
        "adiabatic": {
            "hamiltonian": "ham_adi",
            "vibronic": "hvib_adi",
            "nac": "nac_adi",
        },
        "diabatic": {
            "hamiltonian": "ham_dia",
            "vibronic": "hvib_dia",
            "nac": "nac_dia",
        },
    }

    def __init__(self, backend: Any):
        self.backend = backend

    def evaluate(
        self,
        traj: Any,
        storage: Any,
        model_fn: Callable,
        rep: str = "adiabatic",
    ) -> Any:
        """
        Evaluate a model and write available Hamiltonian tensors to storage.

        model_fn receives R, P, storage, and traj and should return a mapping
        with keys such as H_adi, H_dia, NAC_adi, dH_adi, or their diabatic
        counterparts. Only keys present in the result are written.
        """

        self._validate_rep(rep)
        idx = traj.tbf_ids
        result = as_result_mapping(
            model_fn(
                R=storage.q[traj.id, idx],
                P=storage.p[traj.id, idx],
                storage=storage,
                traj=traj,
            )
        )

        self._ensure_derivative_storage(storage, result)

        for result_name, storage_name in self._RESULT_TO_STORAGE.items():
            value = result.get(result_name)
            if value is None:
                continue
            self._write(storage, traj, storage_name, value)

        self._derive_missing_representation_data(storage, traj, result, rep)

        if result.get(f"Hvib_{self._rep_suffix(rep)}") is None:
            self.build_hvib(storage, traj, reps=(rep,))
        return storage

    def build_state(
        self,
        traj: Any,
        storage: Any,
        model_fn: Callable,
        rep: str = "adiabatic",
        T: Optional[Any] = None,
        apply_ssy: bool = False,
    ) -> Any:
        """
        Full Hamiltonian construction pipeline.

        The name is retained as a compatibility alias, but the returned object
        is TensorStorage, not a separate HamiltonianState snapshot.
        """

        self.evaluate(traj, storage, model_fn, rep)

        if apply_ssy:
            self.apply_ssy(storage, traj)
            self.build_hvib(storage, traj, reps=(rep,))

        if T is not None:
            self.apply_ld(storage, traj, T, reps=(rep,))
            self.build_hvib(storage, traj, reps=(rep,))

        return storage

    def active_matrix(
        self,
        storage: Any,
        traj: Any,
        rep: str = "adiabatic",
        kind: str = "vibronic",
    ) -> Any:
        """Return the active Hamiltonian matrix slice for a trajectory."""

        field = self._field_for(rep, kind)
        return getattr(storage, field)[traj.id, traj.tbf_ids]

    def apply_ld(
        self,
        storage: Any,
        traj: Any,
        T: Any,
        reps: Iterable[str] = ("adiabatic", "diabatic"),
    ) -> Any:
        """
        Apply a local-diabatization rotation to stored Hamiltonian tensors.

        H' = T† H T and NAC' = T† NAC T are written back to TensorStorage.
        Vibronic Hamiltonians are rebuilt by build_state() after this rotation.
        """

        bd = self.backend
        Tdag = bd.conjugate_transpose(T)

        for rep in reps:
            self._validate_rep(rep)
            for kind in ("hamiltonian", "nac"):
                field = self._field_for(rep, kind)
                matrix = getattr(storage, field)
                if matrix is None:
                    continue
                current = matrix[traj.id, traj.tbf_ids]
                matrix[traj.id, traj.tbf_ids] = bd.einsum(
                    "...ij,...jk,...kl->...il",
                    Tdag,
                    current,
                    T,
                )

        self._write(storage, traj, "basis_transform", T)
        return storage

    def apply_ssy(self, storage: Any, traj: Any) -> Any:
        """
        Placeholder for SSY correction logic.

        Future implementations should modify TensorStorage slices in place.
        """

        return storage

    def build_hvib(
        self,
        storage: Any,
        traj: Any,
        reps: Iterable[str] = ("adiabatic", "diabatic"),
    ) -> Any:
        """Construct stored vibronic Hamiltonians as Hvib = H - i * NAC."""

        for rep in reps:
            self._validate_rep(rep)
            ham = getattr(storage, self._field_for(rep, "hamiltonian"))
            hvib = getattr(storage, self._field_for(rep, "vibronic"))
            nac = getattr(storage, self._field_for(rep, "nac"))
            if ham is None or hvib is None:
                continue

            idx = traj.tbf_ids
            value = ham[traj.id, idx]
            if nac is not None:
                value = value - 1j * nac[traj.id, idx]
            hvib[traj.id, idx] = value

        return storage

    def _write(self, storage: Any, traj: Any, storage_name: str, value: Any) -> None:
        target = getattr(storage, storage_name)
        if target is None:
            raise AttributeError(
                f"TensorStorage field '{storage_name}' has not been allocated"
            )
        target[traj.id, traj.tbf_ids] = value

    def _ensure_derivative_storage(self, storage: Any, result: dict) -> None:
        needs_d1 = any(
            result.get(name) is not None
            for name in ("DC1_adi", "DC1_dia", "dH_adi", "dH_dia")
        )
        needs_d2 = any(
            result.get(name) is not None
            for name in ("d2H_adi", "d2H_dia")
        )
        if (needs_d1 or needs_d2) and (
            storage.dc1_adi is None
            or storage.dc1_dia is None
            or storage.d1ham_adi is None
            or storage.d1ham_dia is None
            or (needs_d2 and (storage.d2ham_adi is None or storage.d2ham_dia is None))
        ):
            storage.allocate_hamiltonian_derivatives(der_lvl=2 if needs_d2 else 1)

    def _derive_missing_representation_data(
        self,
        storage: Any,
        traj: Any,
        result: dict,
        rep: str,
    ) -> None:
        suffix = self._rep_suffix(rep)

        if rep == "adiabatic" and result.get("H_adi") is None and result.get("H_dia") is not None:
            from .adiabatic import compute_adiabatic_from_diabatic

            der_lvl = 1 if _has_derivative_data(result, "dia") else 0
            compute_adiabatic_from_diabatic(storage, traj, der_lvl=der_lvl)

        if result.get(f"NAC_{suffix}") is None:
            self._build_velocity_projected_nac(storage, traj, rep)

    def _build_velocity_projected_nac(self, storage: Any, traj: Any, rep: str) -> None:
        dc1_name = f"dc1_{self._rep_suffix(rep)}"
        nac_name = self._field_for(rep, "nac")
        dc1 = getattr(storage, dc1_name, None)
        nac = getattr(storage, nac_name, None)
        if dc1 is None or nac is None or storage.p is None:
            return

        idx = traj.tbf_ids
        velocity = storage.p[traj.id, idx]
        if getattr(storage, "iM", None) is not None:
            velocity = velocity * storage.iM[traj.id, idx]
        nac[traj.id, idx] = self.backend.einsum(
            "...d,...dij->...ij",
            velocity,
            dc1[traj.id, idx],
        )

    def _field_for(self, rep: str, kind: str) -> str:
        self._validate_rep(rep)
        try:
            return self._REP_FIELDS[rep][kind]
        except KeyError as exc:
            raise ValueError(f"Unknown Hamiltonian kind: {kind}") from exc

    @staticmethod
    def _rep_suffix(rep: str) -> str:
        return "adi" if rep == "adiabatic" else "dia"

    @staticmethod
    def _validate_rep(rep: str) -> None:
        if rep not in ("adiabatic", "diabatic"):
            raise ValueError("rep must be 'adiabatic' or 'diabatic'")


def _has_derivative_data(result: dict, suffix: str) -> bool:
    return (
        result.get(f"DC1_{suffix}") is not None
        or result.get(f"dH_{suffix}") is not None
        or result.get(f"d2H_{suffix}") is not None
    )
