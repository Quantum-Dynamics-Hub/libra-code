from __future__ import annotations

from typing import Any, Callable, Optional

from .state import HamiltonianState


class HamiltonianEngine:
    """
    Stateless Hamiltonian construction and transformation engine.

    Responsibilities:
    ----------------------------------------
    1. Build electronic Hamiltonians (dia / adi)
    2. Build vibronic Hamiltonians
    3. Apply local diabatization (LD)
    4. Apply corrections (SSY hooks)
    5. Provide unified pipeline for TDSE propagation setup

    This replaces:
        - nHamiltonian::compute_adiabatic
        - nHamiltonian::compute_diabatic
        - transform_all()
        - hvib reconstruction logic
    """

    def __init__(self, backend: Any):
        self.backend = backend

    # ============================================================
    # 1. MAIN ENTRY: BUILD FULL HAMILTONIAN STATE
    # ============================================================

    def evaluate(
        self,
        traj: Any,
        storage: Any,
        model_fn: Callable,
        rep: str = "adiabatic",
    ) -> HamiltonianState:
        """
        Build HamiltonianState from a physical model.

        Parameters
        ----------
        traj : Trajectory
        storage : TensorStorage
        model_fn : callable
            User-defined electronic structure model
        rep : str
            "adiabatic" or "diabatic"
        """

        idx = traj.tbf_ids

        # --------------------------------------------------------
        # Load nuclear DOFs
        # --------------------------------------------------------
        R = storage.q[traj.id, idx]
        P = storage.p[traj.id, idx]

        # --------------------------------------------------------
        # Call electronic structure model
        # --------------------------------------------------------
        result = model_fn(R=R, P=P, storage=storage, traj=traj)

        # --------------------------------------------------------
        # Extract Hamiltonian components
        # --------------------------------------------------------
        H_adi = result.get("H_adi")
        H_dia = result.get("H_dia")

        S_adi = result.get("S_adi")
        S_dia = result.get("S_dia")

        NAC_adi = result.get("NAC_adi")
        NAC_dia = result.get("NAC_dia")

        DC1_adi = result.get("DC1_adi")
        DC1_dia = result.get("DC1_dia")

        dH_adi = result.get("dH_adi")
        dH_dia = result.get("dH_dia")

        d2H_adi = result.get("d2H_adi")
        d2H_dia = result.get("d2H_dia")

        Hvib_adi = result.get("Hvib_adi")
        Hvib_dia = result.get("Hvib_dia")

        basis_transform = result.get("basis_transform")

        # --------------------------------------------------------
        # Build immutable state
        # --------------------------------------------------------
        return HamiltonianState(
            rep=rep,

            H_adi=H_adi,
            H_dia=H_dia,

            Hvib_adi=Hvib_adi,
            Hvib_dia=Hvib_dia,

            S_adi=S_adi,
            S_dia=S_dia,

            NAC_adi=NAC_adi,
            NAC_dia=NAC_dia,

            DC1_adi=DC1_adi,
            DC1_dia=DC1_dia,

            dH_adi=dH_adi,
            dH_dia=dH_dia,

            d2H_adi=d2H_adi,
            d2H_dia=d2H_dia,

            basis_transform=basis_transform,
        )

    # ============================================================
    # 2. LOCAL DIABATIZATION TRANSFORMATION
    # ============================================================

    def apply_ld(self, state: HamiltonianState, T: Any) -> HamiltonianState:
        """
        Apply local diabatization transform:

            H' = T† H T
            C' = T† C   (handled in propagation layer)
        """

        bd = self.backend
        Tdag = bd.conjugate_transpose(T)

        H_adi_new = state.H_adi
        H_dia_new = state.H_dia
        NAC_adi_new = state.NAC_adi
        NAC_dia_new = state.NAC_dia

        # --------------------------------------------------------
        # Transform adiabatic Hamiltonian
        # --------------------------------------------------------
        if H_adi_new is not None:
            H_adi_new = bd.einsum(
                "nij,njk,nkl->nil",
                Tdag, H_adi_new, T
            )

        # --------------------------------------------------------
        # Transform diabatic Hamiltonian
        # --------------------------------------------------------
        if H_dia_new is not None:
            H_dia_new = bd.einsum(
                "nij,njk,nkl->nil",
                Tdag, H_dia_new, T
            )

        # --------------------------------------------------------
        # Transform NACs if present
        # --------------------------------------------------------
        if NAC_adi_new is not None:
            NAC_adi_new = bd.einsum(
                "nij,njk,nkl->nil",
                Tdag, NAC_adi_new, T
            )

        if NAC_dia_new is not None:
            NAC_dia_new = bd.einsum(
                "nij,njk,nkl->nil",
                Tdag, NAC_dia_new, T
            )

        return state.replace(
            H_adi=H_adi_new,
            H_dia=H_dia_new,
            NAC_adi=NAC_adi_new,
            NAC_dia=NAC_dia_new,
            basis_transform=T
        )

    # ============================================================
    # 3. SSY / EMPIRICAL CORRECTIONS (HOOK)
    # ============================================================

    def apply_ssy(self, state: HamiltonianState, traj: Any) -> HamiltonianState:
        """
        Placeholder for SSY correction logic.

        In C++ this was:
            SSY_correction(H, dyn_var, ham, itraj)

        Here we isolate it as a hook.
        """

        # Example (user-defined later):
        # modify H based on population, decoherence, etc.
        return state

    # ============================================================
    # 4. VIBRONIC HAMILTONIAN CONSTRUCTION
    # ============================================================

    def build_hvib(self, state: HamiltonianState) -> HamiltonianState:
        """
        Construct vibronic Hamiltonian:

            Hvib = H - i * NAC
        """

        Hvib_adi = state.Hvib_adi
        Hvib_dia = state.Hvib_dia

        # --------------------------------------------------------
        # Adiabatic vibronic Hamiltonian
        # --------------------------------------------------------
        if Hvib_adi is None and state.H_adi is not None:
            Hvib_adi = state.H_adi
            if state.NAC_adi is not None:
                Hvib_adi = state.H_adi - 1j * state.NAC_adi

        # --------------------------------------------------------
        # Diabatic vibronic Hamiltonian
        # --------------------------------------------------------
        if Hvib_dia is None and state.H_dia is not None:
            Hvib_dia = state.H_dia
            if state.NAC_dia is not None:
                Hvib_dia = state.H_dia - 1j * state.NAC_dia

        return state.replace(
            Hvib_adi=Hvib_adi,
            Hvib_dia=Hvib_dia
        )

    # ============================================================
    # 5. FULL PIPELINE (MOST USEFUL FUNCTION)
    # ============================================================

    def build_state(
        self,
        traj: Any,
        storage: Any,
        model_fn: Callable,
        rep: str = "adiabatic",
        T: Optional[Any] = None,
        apply_ssy: bool = False,
    ) -> HamiltonianState:
        """
        Full Hamiltonian construction pipeline.

        This replaces:
            - compute_adiabatic()
            - compute_diabatic()
            - hvib construction
            - LD transforms
            - SSY hooks
        """

        state = self.evaluate(traj, storage, model_fn, rep)

        state = self.build_hvib(state)

        if apply_ssy:
            state = self.apply_ssy(state, traj)

        if T is not None:
            state = self.apply_ld(state, T)

        return state
