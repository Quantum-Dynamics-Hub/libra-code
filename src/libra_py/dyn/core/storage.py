# storage/tensor_storage.py

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Any


# ============================================================
# Tensor storage
# ============================================================

@dataclass
class TensorStorage:
    """
    Central backend-agnostic storage.

    TensorStorage is the tensor-backed Python analogue of the storage roles
    split across the C++ dyn_variables and nHamiltonian classes. It stores
    arrays for many independent trajectories, each of which may contain a
    dynamic list of trajectory basis functions (TBFs).

    Dimensions:

        ntraj   -> trajectories
        ntbf    -> trajectory basis functions (dynamic)
        ndof    -> nuclear degrees of freedom
        nstates -> electronic states

        C++ naming guide:

        dyn_variables.h:
            ndia, nadi, ndof, ntraj
            q, p, iM, f
            ampl_adi, ampl_dia
            q_mm, p_mm, proj_adi
            dm_adi, dm_dia
            act_states, act_states_dia
            basis_transform
            dR, dP
            reversal_events, coherence_time
            dm_adi_prev, dm_dia_prev, fssh3_errors
            SHXF/MQCXF auxiliary trajectory fields
            TCNBRA, QTSH, KC-RPMD, and simple-decoherence fields

        nHamiltonian.h:
            nnucl, ndia, nadi
            ham_adi, ham_dia
            hvib_adi, hvib_dia
            nac_adi, nac_dia
            dc1_adi, dc1_dia
            d1ham_adi, d1ham_dia
            d2ham_adi, d2ham_dia
            ovlp_dia
            time_overlap_adi, time_overlap_dia
            basis_transform
            ordering_adi, cum_phase_corr
            eigen_algo, phase_corr_ovlp_tol, gs_kinetic_energy

    Notes:

        The C++ API distinguishes ndia and nadi and stores many quantities
        separately in diabatic and adiabatic representations. This first
        Python storage pass uses one nstates dimension for both spaces but
        keeps the C++ names so adiabatic and diabatic data can coexist.

        Non-array fields such as tcnbra_thermostats remain Python lists. Status
        flags mirror the C++ allocation-status fields and are set to 1 by the
        allocation function that initializes the corresponding tensor group.

    Shapes:

        q            (ntraj, ntbf, ndof)
        p            (ntraj, ntbf, ndof)
        iM           (ntraj, ntbf, ndof)
        f            (ntraj, ntbf, ndof)

        ampl_adi     (ntraj, ntbf, nstates)
        ampl_dia     (ntraj, ntbf, nstates)
        q_mm         (ntraj, ntbf, nstates)
        p_mm         (ntraj, ntbf, nstates)
        proj_adi     (ntraj, ntbf, nstates, nstates)
        dm_adi       (ntraj, ntbf, nstates, nstates)
        dm_dia       (ntraj, ntbf, nstates, nstates)
        act_states   (ntraj, ntbf)
        act_states_dia (ntraj, ntbf)

        ham_adi      (ntraj, ntbf, nstates, nstates)
        ham_dia      (ntraj, ntbf, nstates, nstates)
        hvib_adi     (ntraj, ntbf, nstates, nstates)
        hvib_dia     (ntraj, ntbf, nstates, nstates)
        nac_adi      (ntraj, ntbf, nstates, nstates)
        nac_dia      (ntraj, ntbf, nstates, nstates)
        basis_transform (ntraj, ntbf, nstates, nstates)
        ovlp_dia     (ntraj, ntbf, nstates, nstates)
        time_overlap_adi (ntraj, ntbf, nstates, nstates)
        time_overlap_dia (ntraj, ntbf, nstates, nstates)
        dc1_*        (ntraj, ntbf, ndof, nstates, nstates)
        d1ham_*      (ntraj, ntbf, ndof, nstates, nstates)
        d2ham_*      (ntraj, ntbf, ndof, ndof, nstates, nstates)

        dR, dP       (ntraj, ntbf, ndof, nstates, nstates)
        q_aux, p_aux (ntraj, ntbf, nstates, ndof)
        simple-decoherence tensors (ntraj, ntbf, nstates, nstates)

    """

    # Array backend, e.g. numpy, torch, jax, or a facade implementing zeros
    # and the linear algebra operations used by the dyn prototype.
    backend: Any

    # Dimension numbers. These correspond to dyn_variables::ntraj,
    # dyn_variables::ndof, and nHamiltonian::nadi/ndia in this single
    # electronic-state-count prototype.
    ntraj: int
    ndof: int
    nstates: int

    # C++ dimension aliases. In this prototype ndia == nadi == nstates.
    ndia: int = field(init=False)
    nadi: int = field(init=False)
    nnucl: int = field(init=False)

    # Initial and allocated TBF slots. ntbf is the logical count;
    # ntbf_capacity is the currently allocated tensor width.
    ntbf_initial: int = 1
    ntbf_capacity: Optional[int] = None

    # Backend/device placement hint for CPU/GPU array implementations.
    device: str = "cpu"

    # Dynamic TBF count currently available to all trajectories.
    ntbf: int = field(init=False)

    # Activity mask for TBF slots. True means a trajectory owns an active TBF
    # in that slot; False means the slot is empty or has been pruned.
    alive: Any = field(init=False)

    # Current MD time step.
    timestep: int = 0

    # nHamiltonian control parameters and scalar metadata.
    eigen_algo: int = 0
    phase_corr_ovlp_tol: float = 0.0
    gs_kinetic_energy: float = 0.0

    # Allocation/status flags following dyn_variables.h conventions.
    electronic_vars_status: int = 0
    nuclear_vars_status: int = 0
    afssh_vars_status: int = 0
    bcsh_vars_status: int = 0
    dish_vars_status: int = 0
    fssh2_vars_status: int = 0
    shxf_vars_status: int = 0
    mqcxf_vars_status: int = 0
    tcnbra_vars_status: int = 0
    qtsh_vars_status: int = 0
    kcrpmd_vars_status: int = 0
    simple_decoherence_vars_status: int = 0

    # nHamiltonian memory-status flags. The C++ convention is:
    # 0 - not allocated, 1 - allocated internally, 2 - allocated externally.
    ovlp_dia_mem_status: int = 0
    dc1_dia_mem_status: int = 0
    ham_dia_mem_status: int = 0
    nac_dia_mem_status: int = 0
    hvib_dia_mem_status: int = 0
    d1ham_dia_mem_status: int = 0
    d2ham_dia_mem_status: int = 0
    dc1_adi_mem_status: int = 0
    ham_adi_mem_status: int = 0
    nac_adi_mem_status: int = 0
    hvib_adi_mem_status: int = 0
    d1ham_adi_mem_status: int = 0
    d2ham_adi_mem_status: int = 0
    basis_transform_mem_status: int = 0
    time_overlap_adi_mem_status: int = 0
    time_overlap_dia_mem_status: int = 0
    cum_phase_corr_mem_status: int = 0

    # -------------------------------------------------
    # Nuclear
    # -------------------------------------------------

    # Nuclear coordinates.
    q: Optional[Any] = None

    # Nuclear momenta.
    p: Optional[Any] = None

    # Inverse nuclear masses.
    iM: Optional[Any] = None

    # Active nuclear forces.
    f: Optional[Any] = None

    # -------------------------------------------------
    # Electronic
    # -------------------------------------------------

    # Electronic amplitudes in adiabatic representation.
    ampl_adi: Optional[Any] = None

    # Electronic amplitudes in diabatic representation.
    ampl_dia: Optional[Any] = None

    # Meyer-Miller electronic coordinate variables:
    # ampl_adi = (q_mm + i * p_mm) / sqrt(2).
    q_mm: Optional[Any] = None

    # Meyer-Miller electronic momentum variables.
    p_mm: Optional[Any] = None

    # Cumulative projection matrices for adiabatic-state ordering.
    proj_adi: Optional[Any] = None

    # Electronic density matrix in adiabatic representation.
    dm_adi: Optional[Any] = None

    # Electronic density matrix in diabatic representation.
    dm_dia: Optional[Any] = None

    # Electronic density matrix in adiabatic representation at previous step.
    dm_adi_prev: Optional[Any] = None

    # Electronic density matrix in diabatic representation at previous step.
    dm_dia_prev: Optional[Any] = None

    # Active adiabatic state index for each TBF.
    act_states: Optional[Any] = None

    # Active diabatic state index for each TBF.
    act_states_dia: Optional[Any] = None

    # -------------------------------------------------
    # Hamiltonians
    # -------------------------------------------------

    # Hamiltonian in diabatic representation: <psi_dia|H|psi_dia>,
    # generally non-diagonal.
    ham_dia: Optional[Any] = None

    # Hamiltonian in adiabatic representation: <psi_adi|H|psi_adi>,
    # diagonal by construction.
    ham_adi: Optional[Any] = None

    # Vibronic Hamiltonian in diabatic representation:
    # hvib_dia = ham_dia - i * hbar * nac_dia.
    hvib_dia: Optional[Any] = None

    # Vibronic Hamiltonian in adiabatic representation:
    # hvib_adi = ham_adi - i * hbar * nac_adi.
    hvib_adi: Optional[Any] = None

    # Nonadiabatic/time-derivative couplings in diabatic representation:
    # <psi_dia|d/dt|psi_dia>.
    nac_dia: Optional[Any] = None

    # Nonadiabatic/time-derivative couplings in adiabatic representation:
    # <psi_adi|d/dt|psi_adi>.
    nac_adi: Optional[Any] = None

    # Basis transformation, the eigenvectors U in |psi_adi> = |psi_dia> * U.
    basis_transform: Optional[Any] = None

    # Overlap of the diabatic basis, <psi_dia|psi_dia>.
    ovlp_dia: Optional[Any] = None

    # First-order derivative coupling matrices in diabatic representation.
    dc1_dia: Optional[Any] = None

    # First-order derivative coupling matrices in adiabatic representation.
    dc1_adi: Optional[Any] = None

    # First derivatives of ham_dia with respect to nuclear coordinates.
    d1ham_dia: Optional[Any] = None

    # First derivatives of ham_adi with respect to nuclear coordinates.
    d1ham_adi: Optional[Any] = None

    # Second derivatives of ham_dia with respect to nuclear coordinates.
    d2ham_dia: Optional[Any] = None

    # Second derivatives of ham_adi with respect to nuclear coordinates.
    d2ham_adi: Optional[Any] = None

    # Time-overlap matrix <psi_i_adi(t)|psi_j_adi(t+dt)>.
    time_overlap_adi: Optional[Any] = None

    # Time-overlap matrix <psi_i_dia(t)|psi_j_dia(t+dt)>.
    time_overlap_dia: Optional[Any] = None

    # Cumulative phase corrections already applied to basis_transform.
    cum_phase_corr: Optional[Any] = None

    # Adiabatic-state ordering/permutation metadata.
    ordering_adi: Optional[Any] = None

    # -------------------------------------------------
    # A-FSSH
    # -------------------------------------------------

    # Moments of coordinates in adiabatic representation.
    dR: Optional[Any] = None

    # Moments of momenta in adiabatic representation.
    dP: Optional[Any] = None

    # -------------------------------------------------
    # BCSH / DISH / FSSH3
    # -------------------------------------------------

    # BCSH reversal event matrix.
    reversal_events: Optional[Any] = None

    # DISH coherence times.
    coherence_time: Optional[Any] = None

    # Various kinds of errors in FSSH3 approach.
    fssh3_errors: Optional[Any] = None

    # -------------------------------------------------
    # SHXF / MQCXF
    # -------------------------------------------------

    # Whether an adiabatic state interacts with others.
    is_mixed: Optional[Any] = None

    # Whether decoherence is turned on for the first time.
    is_first: Optional[Any] = None

    # Whether to fix an auxiliary trajectory.
    is_fixed: Optional[Any] = None

    # Whether to keep auxiliary momenta.
    is_keep: Optional[Any] = None

    # State-wise auxiliary nuclear coordinates.
    q_aux: Optional[Any] = None

    # State-wise auxiliary nuclear momenta.
    p_aux: Optional[Any] = None

    # Auxiliary momenta from the previous step.
    p_aux_old: Optional[Any] = None

    # Spatial derivative of the phase of coefficients.
    nab_phase: Optional[Any] = None

    # Previous-step spatial derivative of the phase of coefficients.
    nab_phase_old: Optional[Any] = None

    # XF Hamiltonian.
    ham_xf: Optional[Any] = None

    # Wave-packet widths based on Gaussian approximation.
    wp_width: Optional[Any] = None

    # Quantum momenta defined as (-1) * grad_nuc |chi| / |chi|.
    p_quant: Optional[Any] = None

    # Exact vector potential.
    VP: Optional[Any] = None

    # Decoherence force in MQCXF.
    f_xf: Optional[Any] = None

    # -------------------------------------------------
    # TCNBRA / QTSH / KC-RPMD
    # -------------------------------------------------

    # Alpha parameters used to scale NACs.
    thermal_correction_factors: Optional[Any] = None

    # Auxiliary thermostats for each trajectory.
    tcnbra_thermostats: list[Any] = field(default_factory=list)

    # Kinetic energies for each trajectory.
    tcnbra_ekin: Optional[Any] = None

    # Nonclassical force in QTSH.
    qtsh_f_nc: Optional[Any] = None

    # Classical auxiliary electronic variable mass.
    m_aux_var: Optional[Any] = None

    # Classical auxiliary electronic variable coordinate.
    y_aux_var: Optional[Any] = None

    # Classical auxiliary electronic variable momentum.
    p_aux_var: Optional[Any] = None

    # Classical auxiliary electronic variable force.
    f_aux_var: Optional[Any] = None

    # -------------------------------------------------
    # Simple decoherence
    # -------------------------------------------------

    # Integrated decoherence factors, reset at accepted hops.
    coherence_factors: Optional[Any] = None

    # Pair-wise coherence clocks.
    coherence_clocks: Optional[Any] = None

    # Previous and current state-pair energy gaps.
    gaps_prev: Optional[Any] = None
    gaps_curr: Optional[Any] = None

    # First-time flags for gap running averages.
    is_first_gap: Optional[Any] = None

    # Running averages for gaps and gap squares.
    mean_gap: Optional[Any] = None
    mean_gap2: Optional[Any] = None

    # Number of averaging steps for each state pair.
    averaging_steps: Optional[Any] = None

    # Gap fluctuations and previous-step fluctuations.
    gap_fluctuations: Optional[Any] = None
    gap_fluctuations_prev: Optional[Any] = None

    # Gap correlation functions.
    gap_correlations: Optional[Any] = None

    # Decoherence rates averaged over trajectories.
    ave_decoherence_rates: Optional[Any] = None

    # -------------------------------------------------

    def __post_init__(self):
        """
        Initialize dynamic TBF bookkeeping and allocate tensor buffers.

        If ntbf_capacity is not supplied, a small spare capacity is chosen so
        early spawning does not immediately reallocate all tensors.
        """

        self.ndia = self.nstates
        self.nadi = self.nstates
        self.nnucl = self.ndof
        self.ntbf = self.ntbf_initial

        if self.ntbf_capacity is None:
            self.ntbf_capacity = max(
                16,
                2 * self.ntbf_initial
            )

        self.allocate()

    # ==================================================
    # Allocation
    # ==================================================

    def allocate(self):
        """
        Allocate all tensor fields using the configured backend.

        All arrays are allocated with shape convention
        (ntraj, ntbf_capacity, ...). The first ntbf_initial slots are marked
        alive for every trajectory; remaining capacity is inactive until
        spawn() marks a slot alive for a particular trajectory.

        Method-specific arrays are left unallocated until the corresponding
        allocate_* method is called.
        """

        self.alive = self.backend.zeros(
            (self.ntraj, self.ntbf_capacity),
            dtype=bool
        )
        self.alive[:, :self.ntbf] = True

        self.allocate_nuclear_vars()
        self.allocate_electronic_vars()
        self.allocate_hamiltonian_vars()

    # ==================================================
    # Allocation groups
    # ==================================================

    def _zeros(self, shape, dtype=float):
        """Allocate a zero-filled backend array."""

        return self.backend.zeros(shape, dtype=dtype)

    def _shape(self, *tail):
        """Return a shape prefixed by (ntraj, ntbf_capacity)."""

        return (self.ntraj, self.ntbf_capacity, *tail)

    def allocate_nuclear_vars(self):
        """
        Allocate general nuclear variables.

        Mirrors dyn_variables::allocate_nuclear_vars() for q, p, iM, and f.
        """

        nd = self.ndof

        self.q = self._zeros(self._shape(nd))
        self.p = self._zeros(self._shape(nd))
        self.iM = self._zeros(self._shape(nd))
        self.f = self._zeros(self._shape(nd))

        self.nuclear_vars_status = 1

    def allocate_electronic_vars(self):
        """
        Allocate general electronic variables.

        Mirrors dyn_variables::allocate_electronic_vars() for amplitudes,
        Meyer-Miller variables, density matrices, active states, projection
        matrices, and basis_transform.
        """

        ns = self.nstates

        self.ampl_adi = self._zeros(
            self._shape(ns),
            dtype=complex
        )

        self.ampl_dia = self._zeros(
            self._shape(ns),
            dtype=complex
        )

        self.q_mm = self._zeros(self._shape(ns))
        self.p_mm = self._zeros(self._shape(ns))

        self.proj_adi = self._zeros(
            self._shape(ns, ns),
            dtype=complex
        )

        self.dm_adi = self._zeros(
            self._shape(ns, ns),
            dtype=complex
        )

        self.dm_dia = self._zeros(
            self._shape(ns, ns),
            dtype=complex
        )

        self.act_states = self._zeros(
            self._shape(),
            dtype=int
        )

        self.act_states_dia = self._zeros(
            self._shape(),
            dtype=int
        )

        self.basis_transform = self._zeros(
            self._shape(ns, ns),
            dtype=complex
        )

        self.electronic_vars_status = 1
        self.basis_transform_mem_status = 1

    def allocate_hamiltonian_vars(self):
        """
        Allocate core nHamiltonian variables.

        This includes Hamiltonians, vibronic Hamiltonians, NACs, overlaps,
        time-overlaps, phase corrections, and ordering metadata. Derivative
        tensors are allocated separately by allocate_hamiltonian_derivatives().
        """

        ns = self.nstates

        self.ham_dia = self._zeros(
            self._shape(ns, ns),
            dtype=complex
        )

        self.ham_adi = self._zeros(
            self._shape(ns, ns),
            dtype=complex
        )

        self.hvib_dia = self._zeros(
            self._shape(ns, ns),
            dtype=complex
        )

        self.hvib_adi = self._zeros(
            self._shape(ns, ns),
            dtype=complex
        )

        self.nac_dia = self._zeros(
            self._shape(ns, ns),
            dtype=complex
        )

        self.nac_adi = self._zeros(
            self._shape(ns, ns),
            dtype=complex
        )

        self.ovlp_dia = self._zeros(
            self._shape(ns, ns),
            dtype=complex
        )

        self.time_overlap_adi = self._zeros(
            self._shape(ns, ns),
            dtype=complex
        )

        self.time_overlap_dia = self._zeros(
            self._shape(ns, ns),
            dtype=complex
        )

        self.cum_phase_corr = self._zeros(
            self._shape(ns, ns),
            dtype=complex
        )

        self.ordering_adi = self._zeros(
            self._shape(ns),
            dtype=int
        )

        self.ovlp_dia_mem_status = 1
        self.ham_dia_mem_status = 1
        self.nac_dia_mem_status = 1
        self.hvib_dia_mem_status = 1
        self.ham_adi_mem_status = 1
        self.nac_adi_mem_status = 1
        self.hvib_adi_mem_status = 1
        self.basis_transform_mem_status = 1
        self.time_overlap_adi_mem_status = 1
        self.time_overlap_dia_mem_status = 1
        self.cum_phase_corr_mem_status = 1

    def allocate_hamiltonian_derivatives(self, der_lvl=1):
        """
        Allocate derivative-coupling and Hamiltonian-derivative tensors.

        Parameters
        ----------
        der_lvl : int, optional
            Highest Hamiltonian derivative level to allocate. Level 1
            allocates dc1_* and d1ham_*. Level 2 also allocates d2ham_*.
        """

        nd = self.ndof
        ns = self.nstates

        self.dc1_dia = self._zeros(
            self._shape(nd, ns, ns),
            dtype=complex
        )

        self.dc1_adi = self._zeros(
            self._shape(nd, ns, ns),
            dtype=complex
        )

        self.d1ham_dia = self._zeros(
            self._shape(nd, ns, ns),
            dtype=complex
        )

        self.d1ham_adi = self._zeros(
            self._shape(nd, ns, ns),
            dtype=complex
        )

        self.dc1_dia_mem_status = 1
        self.dc1_adi_mem_status = 1
        self.d1ham_dia_mem_status = 1
        self.d1ham_adi_mem_status = 1

        if der_lvl >= 2:
            self.d2ham_dia = self._zeros(
                self._shape(nd, nd, ns, ns),
                dtype=complex
            )

            self.d2ham_adi = self._zeros(
                self._shape(nd, nd, ns, ns),
                dtype=complex
            )

            self.d2ham_dia_mem_status = 1
            self.d2ham_adi_mem_status = 1

    def allocate_afssh(self):
        """
        Allocate A-FSSH coordinate and momentum moment tensors.

        Mirrors dyn_variables::allocate_afssh() for dR and dP.
        """

        self.dR = self._zeros(
            self._shape(self.ndof, self.nstates, self.nstates),
            dtype=complex
        )

        self.dP = self._zeros(
            self._shape(self.ndof, self.nstates, self.nstates),
            dtype=complex
        )

        self.afssh_vars_status = 1

    def allocate_bcsh(self):
        """
        Allocate BCSH reversal event data.

        Mirrors dyn_variables::allocate_bcsh().
        """

        self.reversal_events = self._zeros(
            self._shape(self.nstates)
        )

        self.bcsh_vars_status = 1

    def allocate_dish(self):
        """
        Allocate DISH coherence-time data.

        Mirrors dyn_variables::allocate_dish().
        """

        self.coherence_time = self._zeros(
            self._shape(self.nstates)
        )

        self.dish_vars_status = 1

    def allocate_fssh2(self):
        """
        Allocate FSSH2 previous-step density matrices.

        Mirrors dyn_variables::allocate_fssh2().
        """

        ns = self.nstates

        self.dm_dia_prev = self._zeros(
            self._shape(ns, ns),
            dtype=complex
        )

        self.dm_adi_prev = self._zeros(
            self._shape(ns, ns),
            dtype=complex
        )

        self.fssh2_vars_status = 1

    def allocate_fssh3(self, nerrors=5):
        """
        Allocate FSSH3 diagnostic error storage.

        Parameters
        ----------
        nerrors : int, optional
            Number of error channels to keep for each trajectory/TBF.
        """

        self.fssh3_errors = self._zeros(
            self._shape(nerrors)
        )

    def allocate_shxf(self):
        """
        Allocate independent-trajectory SHXF auxiliary variables.

        Mirrors dyn_variables::allocate_shxf().
        """

        ns = self.nstates
        nd = self.ndof

        self.is_mixed = self._zeros(self._shape(ns), dtype=int)
        self.is_first = self._zeros(self._shape(ns), dtype=int)
        self.is_fixed = self._zeros(self._shape(ns), dtype=int)
        self.is_keep = self._zeros(self._shape(ns), dtype=int)

        self.q_aux = self._zeros(self._shape(ns, nd))
        self.p_aux = self._zeros(self._shape(ns, nd))
        self.p_aux_old = self._zeros(self._shape(ns, nd))
        self.nab_phase = self._zeros(self._shape(ns, nd))
        self.nab_phase_old = self._zeros(self._shape(ns, nd))

        self.ham_xf = self._zeros(
            self._shape(ns, ns),
            dtype=complex
        )

        self.wp_width = self._zeros(self._shape(nd))
        self.p_quant = self._zeros(self._shape(nd))
        self.VP = self._zeros(self._shape(nd))

        self.shxf_vars_status = 1

    def allocate_mqcxf(self):
        """
        Allocate MQCXF-specific decoherence force data.

        Mirrors dyn_variables::allocate_mqcxf().
        """

        self.f_xf = self._zeros(self._shape(self.ndof))
        self.mqcxf_vars_status = 1

    def allocate_tcnbra(self, thermostats=None):
        """
        Allocate thermally-corrected NBRA data.

        Parameters
        ----------
        thermostats : list, optional
            Optional thermostat objects to associate with trajectories.
        """

        self.thermal_correction_factors = self._zeros(self.nstates)
        self.tcnbra_thermostats = list(thermostats or [])
        self.tcnbra_ekin = self._zeros(self._shape())

        self.tcnbra_vars_status = 1

    def allocate_qtsh(self):
        """
        Allocate QTSH nonclassical force data.

        Mirrors dyn_variables::allocate_qtsh().
        """

        self.qtsh_f_nc = self._zeros(self._shape(self.ndof))
        self.qtsh_vars_status = 1

    def allocate_kcrpmd(self):
        """
        Allocate KC-RPMD auxiliary electronic variables.

        Mirrors dyn_variables::allocate_kcrpmd().
        """

        ns = self.nstates

        self.m_aux_var = self._zeros(self._shape(ns))
        self.y_aux_var = self._zeros(self._shape(ns))
        self.p_aux_var = self._zeros(self._shape(ns))
        self.f_aux_var = self._zeros(self._shape(ns))

        self.kcrpmd_vars_status = 1

    def allocate_simple_decoherence(self):
        """
        Allocate pair-wise simple-decoherence history tensors.

        Mirrors dyn_variables::allocate_simple_decoherence().
        """

        ns = self.nstates

        self.coherence_factors = self._zeros(self._shape(ns, ns))
        self.coherence_clocks = self._zeros(self._shape(ns, ns))
        self.gaps_prev = self._zeros(self._shape(ns, ns))
        self.gaps_curr = self._zeros(self._shape(ns, ns))
        self.is_first_gap = self._zeros(self._shape(ns, ns), dtype=int)
        self.mean_gap = self._zeros(self._shape(ns, ns))
        self.mean_gap2 = self._zeros(self._shape(ns, ns))
        self.averaging_steps = self._zeros(self._shape(ns, ns))
        self.gap_fluctuations = self._zeros(self._shape(ns, ns))
        self.gap_fluctuations_prev = self._zeros(self._shape(ns, ns))
        self.gap_correlations = self._zeros(self._shape(ns, ns))
        self.ave_decoherence_rates = self._zeros((ns, ns))

        self.simple_decoherence_vars_status = 1

    # ==================================================
    # Active TBF view
    # ==================================================

    def active_mask(self):
        """
        Return the TBF activity mask.

        Returns
        -------
        alive : array-like, shape (ntraj, ntbf_capacity)
            Boolean mask identifying which TBF slots are active for each
            trajectory.
        """

        return self.alive

    # ==================================================
    # Spawn
    # ==================================================

    def spawn(self, traj):
        """
        Activate a new TBF slot for one trajectory.

        Parameters
        ----------
        traj : int
            Trajectory index in the first tensor dimension.

        Returns
        -------
        int
            The newly activated TBF slot index.
        """

        if self.ntbf >= self.ntbf_capacity:
            self.expand()

        idx = self.ntbf

        self.alive[traj, idx] = True

        self.ntbf += 1

        return idx

    # ==================================================
    # Kill
    # ==================================================

    def prune(self, traj, tbf):
        """
        Mark one TBF slot inactive for one trajectory.

        Parameters
        ----------
        traj : int
            Trajectory index.
        tbf : int
            TBF slot index to deactivate.
        """

        self.alive[traj, tbf] = False

    # ==================================================
    # Expand
    # ==================================================

    def expand(self):
        """
        Double the allocated TBF capacity of every tensor field.

        Existing data are copied into the leading TBF slots of the new arrays.
        Non-array attributes, scalar metadata, and backend objects are skipped.
        """

        old = self.ntbf_capacity
        self.ntbf_capacity *= 2

        for name, value in vars(self).items():

            shape = getattr(value, "shape", None)
            try:
                shape = list(shape)
            except TypeError:
                continue

            if len(shape) < 2:
                continue

            if shape[0] != self.ntraj or shape[1] != old:
                continue

            shape[1] = self.ntbf_capacity

            new = self.backend.zeros(
                tuple(shape),
                dtype=value.dtype
            )

            new[:, :old] = value

            setattr(self, name, new)

    # ==================================================
    # Compact
    # ==================================================

    def compact(self):
        """
        Remove TBF slots that are inactive for all trajectories.

        The TBF dimension is filtered by alive.any(axis=0). This preserves any
        slot used by at least one trajectory and drops slots that are globally
        inactive.
        """

        keep = self.alive.any(axis=0)

        for name, value in vars(self).items():

            shape = getattr(value, "shape", None)
            try:
                shape = list(shape)
            except TypeError:
                continue

            if len(shape) < 2:
                continue

            if shape[0] != self.ntraj or shape[1] != len(keep):
                continue

            setattr(
                self,
                name,
                value[:, keep]
            )

        self.ntbf = int(keep.sum())
        self.ntbf_capacity = self.ntbf
