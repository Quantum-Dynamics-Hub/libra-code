"""
Views into TensorStorage.

TensorView gives one trajectory/TBF a lightweight window into the shared
TensorStorage arrays. It never copies data; properties return writable slices.
"""

from .storage import TensorStorage


class TensorView:
    """
    Lightweight views into TensorStorage for a single TBF.
    """

    def __init__(self, storage: TensorStorage, traj_id: int, tbf_id: int):
        self._s = storage
        self._traj = traj_id
        self._tbf = tbf_id

    # ==================================================
    # Generic access
    # ==================================================

    @property
    def storage(self):
        """The backing TensorStorage object."""

        return self._s

    @property
    def traj_id(self):
        """Trajectory index addressed by this view."""

        return self._traj

    @property
    def tbf_id(self):
        """TBF slot index addressed by this view."""

        return self._tbf

    def get(self, name):
        """
        Return one storage field sliced to this trajectory/TBF.

        Parameters
        ----------
        name : str
            TensorStorage member name.

        Raises
        ------
        AttributeError
            If the field does not exist or has not been allocated.
        """

        value = getattr(self._s, name)
        if value is None:
            raise AttributeError(
                f"TensorStorage field '{name}' has not been allocated"
            )
        return value[self._traj, self._tbf]

    def set(self, name, value):
        """
        Assign one storage field slice for this trajectory/TBF.
        """

        target = getattr(self._s, name)
        if target is None:
            raise AttributeError(
                f"TensorStorage field '{name}' has not been allocated"
            )
        target[self._traj, self._tbf] = value

    def maybe(self, name):
        """
        Return a sliced field or None if that field is unallocated.
        """

        value = getattr(self._s, name)
        if value is None:
            return None
        return value[self._traj, self._tbf]

    def as_dict(self, names):
        """
        Return selected allocated fields as a dictionary of writable slices.
        """

        return {name: self.get(name) for name in names}

    # ==================================================
    # Core nuclear variables
    # ==================================================

    @property
    def q(self): return self._s.q[self._traj, self._tbf]
    @property
    def p(self): return self._s.p[self._traj, self._tbf]
    @property
    def iM(self): return self._s.iM[self._traj, self._tbf]
    @property
    def f(self): return self._s.f[self._traj, self._tbf]

    # ==================================================
    # Core electronic variables
    # ==================================================

    @property
    def ampl_adi(self): return self._s.ampl_adi[self._traj, self._tbf]
    @property
    def ampl_dia(self): return self._s.ampl_dia[self._traj, self._tbf]
    @property
    def q_mm(self): return self._s.q_mm[self._traj, self._tbf]
    @property
    def p_mm(self): return self._s.p_mm[self._traj, self._tbf]
    @property
    def proj_adi(self): return self._s.proj_adi[self._traj, self._tbf]
    @property
    def dm_adi(self): return self._s.dm_adi[self._traj, self._tbf]
    @property
    def dm_dia(self): return self._s.dm_dia[self._traj, self._tbf]
    @property
    def act_states(self): return self._s.act_states[self._traj, self._tbf]
    @property
    def act_states_dia(self): return self._s.act_states_dia[self._traj, self._tbf]

    @property
    def dm_adi_prev(self): return self.get("dm_adi_prev")
    @property
    def dm_dia_prev(self): return self.get("dm_dia_prev")

    # ==================================================
    # Core Hamiltonian variables
    # ==================================================

    @property
    def ham_adi(self): return self._s.ham_adi[self._traj, self._tbf]
    @property
    def ham_dia(self): return self._s.ham_dia[self._traj, self._tbf]
    @property
    def hvib_adi(self): return self._s.hvib_adi[self._traj, self._tbf]
    @property
    def hvib_dia(self): return self._s.hvib_dia[self._traj, self._tbf]
    @property
    def nac_adi(self): return self._s.nac_adi[self._traj, self._tbf]
    @property
    def nac_dia(self): return self._s.nac_dia[self._traj, self._tbf]
    @property
    def basis_transform(self): return self._s.basis_transform[self._traj, self._tbf]
    @property
    def ovlp_dia(self): return self._s.ovlp_dia[self._traj, self._tbf]
    @property
    def time_overlap_adi(self): return self._s.time_overlap_adi[self._traj, self._tbf]
    @property
    def time_overlap_dia(self): return self._s.time_overlap_dia[self._traj, self._tbf]
    @property
    def cum_phase_corr(self): return self._s.cum_phase_corr[self._traj, self._tbf]
    @property
    def ordering_adi(self): return self._s.ordering_adi[self._traj, self._tbf]

    @property
    def dc1_dia(self): return self.get("dc1_dia")
    @property
    def dc1_adi(self): return self.get("dc1_adi")
    @property
    def d1ham_dia(self): return self.get("d1ham_dia")
    @property
    def d1ham_adi(self): return self.get("d1ham_adi")
    @property
    def d2ham_dia(self): return self.get("d2ham_dia")
    @property
    def d2ham_adi(self): return self.get("d2ham_adi")

    # ==================================================
    # Method-specific variables
    # ==================================================

    @property
    def dR(self): return self.get("dR")
    @property
    def dP(self): return self.get("dP")

    @property
    def reversal_events(self): return self.get("reversal_events")
    @property
    def coherence_time(self): return self.get("coherence_time")
    @property
    def fssh3_errors(self): return self.get("fssh3_errors")

    @property
    def is_mixed(self): return self.get("is_mixed")
    @property
    def is_first(self): return self.get("is_first")
    @property
    def is_fixed(self): return self.get("is_fixed")
    @property
    def is_keep(self): return self.get("is_keep")
    @property
    def q_aux(self): return self.get("q_aux")
    @property
    def p_aux(self): return self.get("p_aux")
    @property
    def p_aux_old(self): return self.get("p_aux_old")
    @property
    def nab_phase(self): return self.get("nab_phase")
    @property
    def nab_phase_old(self): return self.get("nab_phase_old")
    @property
    def ham_xf(self): return self.get("ham_xf")
    @property
    def wp_width(self): return self.get("wp_width")
    @property
    def p_quant(self): return self.get("p_quant")
    @property
    def VP(self): return self.get("VP")
    @property
    def f_xf(self): return self.get("f_xf")

    @property
    def tcnbra_ekin(self): return self.get("tcnbra_ekin")
    @property
    def qtsh_f_nc(self): return self.get("qtsh_f_nc")
    @property
    def m_aux_var(self): return self.get("m_aux_var")
    @property
    def y_aux_var(self): return self.get("y_aux_var")
    @property
    def p_aux_var(self): return self.get("p_aux_var")
    @property
    def f_aux_var(self): return self.get("f_aux_var")

    @property
    def coherence_factors(self): return self.get("coherence_factors")
    @property
    def coherence_clocks(self): return self.get("coherence_clocks")
    @property
    def gaps_prev(self): return self.get("gaps_prev")
    @property
    def gaps_curr(self): return self.get("gaps_curr")
    @property
    def is_first_gap(self): return self.get("is_first_gap")
    @property
    def mean_gap(self): return self.get("mean_gap")
    @property
    def mean_gap2(self): return self.get("mean_gap2")
    @property
    def averaging_steps(self): return self.get("averaging_steps")
    @property
    def gap_fluctuations(self): return self.get("gap_fluctuations")
    @property
    def gap_fluctuations_prev(self): return self.get("gap_fluctuations_prev")
    @property
    def gap_correlations(self): return self.get("gap_correlations")

    # ==================================================
    # Grouped views
    # ==================================================

    def nuclear(self):
        """Return the core nuclear slices."""

        return self.as_dict(["q", "p", "iM", "f"])

    def electronic(self):
        """Return the core electronic slices."""

        return self.as_dict([
            "ampl_adi",
            "ampl_dia",
            "q_mm",
            "p_mm",
            "proj_adi",
            "dm_adi",
            "dm_dia",
            "act_states",
            "act_states_dia",
        ])

    def hamiltonian(self):
        """Return the core Hamiltonian slices."""

        return self.as_dict([
            "ham_adi",
            "ham_dia",
            "hvib_adi",
            "hvib_dia",
            "nac_adi",
            "nac_dia",
            "basis_transform",
            "ovlp_dia",
            "time_overlap_adi",
            "time_overlap_dia",
            "cum_phase_corr",
            "ordering_adi",
        ])

    def optional(self):
        """Return allocated optional/method-specific slices."""

        names = [
            "dm_adi_prev",
            "dm_dia_prev",
            "dc1_dia",
            "dc1_adi",
            "d1ham_dia",
            "d1ham_adi",
            "d2ham_dia",
            "d2ham_adi",
            "dR",
            "dP",
            "reversal_events",
            "coherence_time",
            "fssh3_errors",
            "is_mixed",
            "is_first",
            "is_fixed",
            "is_keep",
            "q_aux",
            "p_aux",
            "p_aux_old",
            "nab_phase",
            "nab_phase_old",
            "ham_xf",
            "wp_width",
            "p_quant",
            "VP",
            "f_xf",
            "tcnbra_ekin",
            "qtsh_f_nc",
            "m_aux_var",
            "y_aux_var",
            "p_aux_var",
            "f_aux_var",
            "coherence_factors",
            "coherence_clocks",
            "gaps_prev",
            "gaps_curr",
            "is_first_gap",
            "mean_gap",
            "mean_gap2",
            "averaging_steps",
            "gap_fluctuations",
            "gap_fluctuations_prev",
            "gap_correlations",
        ]
        return {
            name: self.maybe(name)
            for name in names
            if getattr(self._s, name) is not None
        }
