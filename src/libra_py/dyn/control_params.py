# *********************************************************************************
# * Copyright (C) 2026 Alexey V. Akimov
# *
# * This file is distributed under the terms of the GNU General Public License
# * as published by the Free Software Foundation, either version 3 of
# * the License, or (at your option) any later version.
# * See the file LICENSE in the root directory of this distribution
# * or <http://www.gnu.org/licenses/>.
# *
# *********************************************************************************/

"""
.. module:: control_params
   :platform: Unix, Windows
   :synopsis: Implements various dynamics control options. Replaces dyn_control_params class
   initially written in C++ and exposed to Python
.. moduleauthor:: Alexey V. Akimov

"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Dict, Any
import numpy as np


@dataclass
class DynControlParams:
    """
    Full Python equivalent of the C++ dyn_control_params class.

    This class is a *pure control container* for nonadiabatic dynamics.
    No physics is implemented here — only algorithmic switches and parameters.

    All C++ MATRIX* / CMATRIX* pointers are mapped to:
        numpy.ndarray | None
    """

    # ==========================================================================
    # Hamiltonian / representation control
    # ==========================================================================

    #: Representation for TDSE propagation
    #: Selects the representation in which nuclear/electronic (Ehrenfest core) dynamics 
    #: is executed. This is the representation chosen as the main one to represent the 
    #: time-dependent wavefunction and the one used to integrate the TD-SE
    #: 0 diabatic wfc
    #: 1 adiabatic wfc [default]
    #: 2 diabatic density matrix
    #: 3 adiabatic density matrix
    #: 4 MMST adiabatic wfc
    rep_tdse: int = 1

    #: Hamiltonian update method
    #: How to update Hamiltonian and which type of Hamiltonian to update
    #: 0 none
    #: 1 recompute diabatic [default]
    #: 2 recompute adiabatic
    ham_update_method: int = 1

    #: Hamiltonian transformation method
    #: How to transform the Hamiltonians between representations
    #: 0 none
    #: 1 dia → adi by diagonalization [default]
    #: 2 dia → adi by stored transform
    #: 3 adi → dia by stored transform
    #: 4 adi → dia by local diabatization
    ham_transform_method: int = 1

    #: Representation for surface hopping
    #: 0 diabatic
    #: 1 adiabatic [default]
    rep_sh: int = 1

    #: Representation for Landau–Zener probabilities
    #: 0 diabatic, diabatic, Eq. 1 of the Belyaev-Lebedev paper [default]
    #: 1 adiabatic, Eq. 3 of the Belyaev-Lebedev paper, 
    #: crossing point is determined by the sign change of the diabatic gap
    #: 2 adiabatic, Eq. 3 of the Belyaev-Lebedev paper, 
    #: crossing point is determined by the sign change of the NAC
    rep_lz: int = 0

    #: Representation for force evaluation
    #: In which representation to compute forces. To clarify - the forces
    #: in both representations shall be equivalent, so this flag actually
    #: selects the type of the properties needed to compute such forces.
    #: For instance, if it is set to 0, we may be using the derivatives
    #: of the diabatic Hamiltonians, diabatic states overlaps, etc.
    #: 0 diabatic
    #: 1 adiabatic [default]
    rep_force: int = 1

    #: Force computation method
    #: 0 none, e.g. we do not really need them
    #: 1 state-specific  as in the TSH or adiabatic (including adiabatic excited states) [default]
    #: 2 Ehrenfest / MMST (SQC)
    #: 3 QTSH
    #: 4 KC-RPMD
    force_method: int = 1

    #: MMST / SQC gamma parameter
    #: Gamma parameter in the Meyer-Miller Hamiltonian of form:
    #:  H = sum_{i,j} H_{ij} [ 1/(2 hbar) * (q_i - i * p_j) * (q_j + i * p_i) - gamma * delta_ij ]
    #: gamma = 0 → Ehrenfest [default]
    sqc_gamma: float = 0.0

    #: Enforce nuclear motion on a fixed electronic state
    #: 0 no [default]
    #: 1 yes
    enforce_state_following: int = 0

    #: Index of enforced state (>= 0)
    #: If we enforce the nuclear dynamics to be on a given state, 
    #: what is the index of that state [any integer >- 0, default = 0 ]
    #: The default value of 0 enforces the nuclear dynamics to be on the ground state. 
    #: This is a convenient way of doing NBRA calculations with model systems without 
    #: the need for pre-computing the trajectories 
    enforced_state_index: int = 0

    #: Time overlap computation
    #: 0 do not compute, perhaps because it was already pre-computed or read
    #: 1 explicitly compute it from the wavefunctions (the Hamiltonian shall have 
    #: the basis_transform variables updated) [default]
    time_overlap_method: int = 1  ### NOTE: In C++ constructor it was actually set to 0 even though
                                  ### NOTE: the docstring said it was 1

    #: NAC update method
    #: 0 none, don't update them (e.g. for simplest NAC)
    #: 1 update according to changed momentum and existing derivative couplings [default]
    #: 2 update according to time-overlaps (only time-derivative NACs)
    nac_update_method: int = 1

    #: Algorithm for time-derivative NACs
    #: -1 none, e.g. we use NACs from somewhere else [default]
    #:  0 HST formula (if nac_update_method==2)
    #:  1 NPI of Meek–Levine (if nac_update_method==2)
    nac_algo: int = -1

    #: Hvib update method
    #: 0 none, don't update them (e.g. if it is read externally)
    #: 1 Hvib = H − iħNAC [default]
    hvib_update_method: int = 1

    #: Apply SSY Hamiltonian modification
    #: Whether to modify the Hamiltonian in the dynamics according the Shenvi-Subotnik-Yang (SSY)
    #: method, see my Chapter Eq. 3.27 Note that this is only applied in the adiabatic representation
    #: 0 no [default]
    #: 1 yes
    do_ssy: int = 0

    #: Adiabatic phase correction
    #: 0 no [default]
    #: 1 yes 
    do_phase_correction: int = 0  ### NOTE: since LD is used by default, we should turn this correction off
                                  ### NOTE: this is different from the C++ prototype

    #: Phase correction overlap threshold
    #: The minimal magnutude of the matrix element for which we'll be computing the phase correction
    #: If the overlap is zero, then we don't really care about the phase, but if it is not, then this
    #: parameter sets out threshold for when we do.  [ default: 1e-3 ]
    phase_correction_tol: float = 1e-3

    #: Phase correction applied directly to NACs
    #: New phase correction, directly applied to NACs. Intended to be used mostly with state_tracking_algo == 4,
    #: although can be useful with other state treacking algorithms. Should not be used together with 
    #: `do_phase_correction`
    #: 0 no [default]
    #: 1 yes
    do_nac_phase_correction: int = 0

    #: State tracking algorithm
    #: -1 use local diabatization (LD) approach, it includes phase correction too [default]
    #:  0 none
    #:  1 method of Kosuke Sato (may fail by getting trapped into an infinite loop)
    #:  2 Munkres-Kuhn (Hungarian)
    #: 21 ChatGPT-generated Munkres-Kuhn (Hungarian) algorithm
    #:  3 experimental stochastic algorithm, the original version with elimination (known problems)
    #: 32 experimental stochastic algorithms with all permutations (too expensive)
    #: 33 the improved stochastic algorithm with good scaling and performance, on par with the mincost
    #:  4 new, experimental force-based tracking
    state_tracking_algo: int = -1

    #: Munkres–Kuhn alpha (selects the range of orbitals included in reordering) [default: 0.0]
    MK_alpha: float = 0.0

    #: Munkres–Kuhn verbosity
    #: 0 no extra output [ default ]
    #: 1 prints extra details on what the algorithm is doing, for debugging
    MK_verbosity: int = 0

    #: Convergence behavior for stochastic reordering
    #: A swtich for stochastic reordering algorithm 3 to choose what happens when an 
    #: acceptable permutation isn't generated in the set number of attempts:
    #: 0 returns the identity permutation (does not require convergence) [default] 
    #: 1 exits and prints an error (requires convergence)
    convergence: int = 0

    #: Maximum attempts in stochastic reordering
    #: The maximum number of hops that an be attempted before either choosing the 
    #: identity or exiting in stochastic reordering algorithm 3 [default: 100]
    max_number_attempts: int = 100

    #: Probability cutoff in stochastic reordering [default: 0.0]
    #: The probability threshold for stochastic state reordering algorithm. 
    #: If a probability for a multi-state stransition is below this value, it will be disregarded and set to 0
    #: The rest of the probabilities will be renormalized
    min_probability_reordering: float = 0.0

    # ==========================================================================
    # Surface hopping
    # ==========================================================================

    #: Surface hopping method
    #: -1 adiabatic dynamics, no hops [default]
    #: 0 Fewest Switches Surface Hopping (FSSH)
    #: 1 Global Flux Surface Hopping (GFSH)
    #: 2 Markov-State Surface Hopping (MSSH)
    #: 3 Landau-Zener (LZ) options
    #: 4 Zhu-Nakamura (ZN) options
    #: 5 DISH
    #: 6 MASH
    #: 7 FSSH2
    #: 8 FSSH3
    #: 9 original GFSH
    tsh_method: int = -1  ### NOTE: again, this is different from C++ prototype where it is set to 0
                          ### NOTE: even though the docstring says it should be -1 

    #: Hop acceptance algorithm
    #: Options to control the acceptance of the proposed hops.
    #:  0 accept all proposed hops  [ default ]
    #: 10 based on adiabatic energy - accept only those hops that can obey the energy conservation with 
    #:    adiabatic potential energies
    #: 11 based on diabatic energy - same as 10, but we use diabatic potential energies
    #: 20 based on derivative coupling vectors - accept only those hops that can obey the energy conservation
    #:    by rescaling nuclear velocities along the directions of derivative couplings for the quantum nuclear DOFs                   
    #: 21 based on difference of state-specific forces - same as 20, but the rescaling is done along the vector
    #:    parallel to the difference of adiabatic forces on initial and target states
    #: 31 accept hops with the probability taken from the quantum Boltzmann distribution
    #: 32 accept hops with the probability taken from the classical Maxwell-Boltzmann distribution
    #: 33 accept hops with the probability taken from the updated quantum Boltzmann distribution (experimental)
    #: 40 based on possibility to conserve energy using tcnbra_ekin variables (for TC-NBRA)
    hop_acceptance_algo: int = 0

    #: Momentum rescaling algorithm
    #: Options to control nuclear momenta changes upon successful or frustrated hops.
    #:   0 don't rescale [ default ]
    #: 100 based on adiabatic energy, don't reverse on frustrated hops
    #: 101 based on adiabatic energy, reverse on frustrated hops
    #: 110 based on diabatic energy, don't reverse on frustrated hops
    #: 111 based on diabatic energy, reverse on frustrated hops
    #: 200 along derivative coupling vectors, don't reverse on frustrated hops
    #: 201 along derivative coupling vectors, reverse on frustrated hops
    #: 210 along difference of state-specific forces, don't reverse on frustrated hops
    #: 211 along difference of state-specific forces, reverse on frustrated hops
    #:  40 does not rescale velocities, but rescales  tcnbra_ekin variables
    momenta_rescaling_algo: int = 0

    #: Jasper–Truhlar frustrated-hop criterion
    #: A flag to turn on/off the Jasper-Truhlar criterion for reversal of velocities on frustrated hops.
    #: According to: Jasper, A. W.; Truhlar, D. G. Chem. Phys. Lett. 2003, 369, 60− 67
    #: 0: don't use this criterion (naive handling)
    #: 1: use it [ default ] - the velocities are reversed along the direction d_{a,j} if
    #:    a) (F_a * d_{a,j}) * (F_j * d_{a,j}) < 0 and 
    #:    b) (v * d_{a,j}) * (F_j * d_{a,j}) < 0 
    #:    where a - is the active state index;  Only in effect, if `momenta_rescaling_algo == 201`
    use_Jasper_Truhlar_criterion: int = 1

    #: Scale hopping probabilities by Boltzmann factor
    #: This is needed for the libra_py/workflows/nbra calculations, where the hop proposal 
    #: probability also includes the factor to account for the hop acceptance probabilities
    #: 0 no [default]
    #: 1 yes
    use_boltz_factor: int = 0

    # ==========================================================================
    # FSSH2 / FSSH3
    # ==========================================================================

    #: Use revised FSSH2
    #: 0 no, use the FSSH2 as it was formulated in the original paper [default]
    #: 1 use the revised version (experimental)
    fssh2_revision: int = 0

    #: FSSH3 density matrix size option
    #: The size of the vectorized density matrix in equations to determine hopping probabilites/fluxes
    #: 0: N elements - only populations; the matrices are overdetermined
    #: 1: N^2 elements - first N elements are populations, then Re and Im parts of upper-triangular coherences
    #:    that is rho_{0,1}, rho_{0,2}, ... rho_{0,N-1}, rho_{1,2}, ... rho_{1,N-1}, ... rho_{N-2,N-1}  [ default ]
    fssh3_size_option: int = 1

    #: FSSH3 approach
    #: The approach to determine the hopping probabilities:
    #: 0: based on master equation, rho(t+dt) = J * rho(t);  J matrix contains hopping probabilities directly  [ default ]
    #: 1: based on kinetic approach, drho/dt = J * rho; J matrix contains fluxes
    fssh3_approach_option: int = 0

    #: FSSH3 matrix decomposition method
    #: The matrix decomposition method for solving the least-squares problem.
    #: In the present implementation, is not used.
    #: 0: bdcSvd
    #: 1: fullPivLu
    #: 2: fullPivHouseholderQr
    #: 3: completeOrthogonalDecomposition [ default ]
    fssh3_decomp_option: int = 3

    #: FSSH3 optimization timestep (a.u.)
    # The time-step of the optimization procedure in the FSSH3 calculations
    # [default: 0.001 a.u.]
    fssh3_dt: float = 0.001

    #: FSSH3 max optimization steps [default: 1000]
    fssh3_max_steps: int = 1000

    #: FSSH3 error tolerance [default: 1e-7]
    fssh3_err_tol: float = 1e-7

    # ==========================================================================
    # QTSH
    # ==========================================================================

    #: Use QTSH
    #: 0 no [default]
    #: 1 yes 
    use_qtsh: int = 0

    #: How to compute nonclassical force in the adiabatic QTSH. Only used with `use_qtsh == 1`
    #: 0 Considering only the first-order force, i.e., off-diagonal Ehrenfest force
    #: 1: The whole force including the second-order term is used [default]
    qtsh_force_option: int = 1

    # ==========================================================================
    # KC-RPMD
    # ==========================================================================

    #: Whether to use KC-RPMD 
    #: 0: don't apply [ default ]
    #: 1: use it 
    use_kcrpmd: int = 0

    #: KC-RPMD free energy parameter eta [default: 6.28]
    kcrpmd_eta: float = 6.28

    #: KC-RPMD kinetic constraint parameter a [default: 0.1]
    kcrpmd_a: float = 0.1

    #: KC-RPMD Heaviside parameter b [default: 1000]
    kcrpmd_b: float = 1000.0

    #: KC-RPMD kinetic constraint switching parameter c [default: 0.5]
    kcrpmd_c: float = 0.5

    #: KC-RPMD free energy switching parameter d [default: 3.0]
    kcrpmd_d: float = 3.0

    #: KC-RPMD Langevin frictional coefficient within donor-acceptor basin [default: 0.0]
    kcrpmd_gamma: float = 0.0

    #: KC-RPMD Langevin frictional coefficient within kinked-pair regime [default: 0.0]
    kcrpmd_gammaKP: float = 0.0

    # ==========================================================================
    # Decoherence
    # ==========================================================================

    #: Decoherence algorithm selector
    #: -1 no decoherence [ default ]
    #:  0 SDM and alike
    #:  1 instantaneous decoherence options (ID-S, ID-A, ID-C)
    #:  2 AFSSH
    #:  3 BCSH of Linjun Wang
    #:  4 MF-SD of Bedard-Hearn, Larsen, Schwartz
    #:  5 SHXF of Min
    #:  6 MQCXF
    #:  7 DISH, rev2023
    #:  8 diabatic IDA, experimental
    #:  9 simple decoherence, experimental
    decoherence_algo: float = -1

    #: SDM norm tolerance
    #: Corresponds to the `tol` parameter in the sdm function. It controls 
    #: how much the norm of the old state can be larger than 1.0  before the 
    #: code stops with the error message [ default: 0.0 ]
    #: Note: only matters if `decoherence_algo == 0`
    sdm_norm_tolerance: float = 0.0

    #: DISH decoherence event sampling
    #: Selects the how to sample decoherence events in the DISH.
    #: 0 compare the coherence time counter with the decoherence time (simplified DISH) 
    #: 1 compare the coherence time counter with the time drawn from the exponential distribution
    #:   with the parameter lambda = 1/decoherence time - this distribution corresponds to 
    #:   the statistics of wait times between the Poisson-distributed events (decoherence)
    #:   This is what the original DISH meant to do [ default ]
    #: Note: only matters if `tsh_method == 5`
    dish_decoherence_event_option: int = 1

    #: Dephasing times calculation type
    #: Type of dephasing times/rates calculation:
    #: -1 set all dephasing rates to zero [ default ]
    #:  0 use the rates read out from the input 
    #:  1 use the energy-based decoherence method (EDC)    
    #:  2 Schwartz - mean-field Force-based decoherence (Schwartz 1), using inv_alpha
    #:  3 Schwartz - pair-wise-based decoherence, (Schwartz 2), using inv_alpha
    #:  4 Schwartz - mean-field Force-based decoherence (Schwartz 1), but using interaction width
    #:  5 Gu-Franco 
    decoherence_times_type: int = -1

    #: Schwartz 1/alpha parameters (ndof, 1)
    #: the parameters used in GWP in computing decoherence rates [default: None]
    schwartz_decoherence_inv_alpha: Optional[np.ndarray] = None

    #: Schwartz interaction widths (ndof, 1)
    #: the parameters for the spatial extent of NAC in computing decoherence rates [default: None]
    schwartz_interaction_width: Optional[np.ndarray] = None

    #: Empirical EDC parameter C [default: 1.0]
    decoherence_C_param: float = 1.0

    #: Empirical EDC parameter eps [default: 0.1 Ha]
    decoherence_eps_param: float = 0.1

    #: Apply dephasing-informed correction of Sifain et al to correct dephasing times:
    #: 0 do not apply [default]
    #: 1 apply
    dephasing_informed: int = 0

    #: Instantaneous decoherence variant
    #: Option to control the instantaneous decoherence methodology, only used with `decoherence_algo == 1`
    #: 0 ID-S (ID at successful hops)
    #: 1 ID-A (ID at attempted hops) [default]
    #:   if the proposed hop is not successful, we project back to the initial state 
    #:   if the proposed hop is accepted - we project onto that state
    #: 2 ID-C - consistent ID - an experimental algorithm
    #: 3 ID-A, new: 
    #:   if the proposed hop is not successful, we project out the proposed states
    #:   if the proposed hop is accepted - we project onto that state
    #: 4 ID-F (ID at frustrated hops), new:
    #:   if the proposed hop is not successful, we project out the proposed states
    #:   but we don't do anything if the hop is successful
    instantaneous_decoherence_variant: int = 1

    #: Collapse option
    #: How to collapse wavefunction amplitudes in the decoherence schemes:
    #: 0 by rescaling the magnitude of the amplitude vector elements, but preserving "phase" [default]
    #: 1 by resetting the amplitudes to 1.0+0.0j. This option changes phase 
    collapse_option: int = 0

    #: User-provided decoherence rates (nstates, nstates) [default: None]
    decoherence_rates: Optional[np.ndarray] = None

    #: A matrix of averaged energy gaps ⟨|Ei − Ej|⟩ (nstates, nstates)
    #: it is needed when `dephasing_informed` option is used [default : None]
    ave_gaps: Optional[np.ndarray] = None

    #: Gaussian wavepacket widths (ndof, 1)
    #: Widths of all DOFs for Gaussian function as an approximation to adiabatic wave packets. 
    #: According to the choice of the Gaussian width approximation,
    #: this parameter has different meanings:         
    #: - A constant width in the fixed-width approximation, that is, `use_td_width == 0`
    #: - The initial width in the free-particle Gaussian wave packet approximation, that is, `use_td_width == 1`
    #: - The interaction width in the Schwarz scheme, that is, `use_td_width == 2`
    #: - No influence on the dynamics since the width will be determined by internal variables in the Subotnik scheme, 
    #:   that is, `use_td_width == 3`
    #: Only used with independent-trajectory XF methods, that is, `decoherence_algo == 5 or 6` [default: None]
    wp_width: Optional[np.ndarray] = None

    #: Gaussian wavepacket velocities (ndof, 1)
    #: The velocity of Gaussian wave packet in the free-particle Gaussian  approximation, 
    #: that is, `use_td_width == 1` 
    #: Only used with independent-trajectory XF methods, that is, with `decoherence_algo == 5 or 6` [default : None]
    wp_v: Optional[np.ndarray] = None

    #: Population threshold for creating/destroying auxiliary trajectories [default: 0.01]
    #: Only used with independent-trajectory XF methods, that is, with `decoherence_algo == 5 or 6`
    coherence_threshold: float = 0.01

    #: Masking parameter for nabla-phase vectors [default: 1e-4 Ha]
    #: Only used with the MQCXF method, that is, `decoherence_algo == 5`
    e_mask: float = 1e-4

    #: Use XF force
    #: Whether to use the decoherence force in MQCXF. The corresponding electronic propagation 
    #: is adjusted for the energy conservation. Only used with `decoherence_algo == 6`
    #: 0 don't use it, so for XF-methods this is only Ehrenfest-like force; EhXF [ default ]
    #: 1 The whole force including the XF-correction; MQCXF 
    use_xf_force: int = 0

    #: Project out forbidden auxiliary trajectories
    #: Whether to project out the density on an auxiliary trajectory when its motion is classically 
    #: forbidden. Only used with independent-trajectory XF methods, that is, `decoherence_algo == 5 or 6`
    #: 0 don't [default]
    #: 1 do
    project_out_aux: int = 0

    #: Turning-point algorithm for auxiliary trajectories. 
    #: Only used with independent-trajectory XF methods, that is, `decoherence_algo == 5 or 6`
    #: 0 no treatment of a turning point
    #: 1 collapse to the active state [default]
    #: 2 fix auxiliary positions of adiabatic states except for the active state
    #: 3 keep auxiliary momenta of adiabatic states except for the active state
    tp_algo: int = 1

    #: Time-dependent Gaussian width option
    #: Whether to use the td Gaussian width for the nuclear wave packet approximation
    #: This option can be considered when it comes to unbounded systems.
    #: This approximation is based on a nuclear wave packet on a free surface:
    #: \sigma_x(t)=\sqrt[\sigma_x(0)^2 + (wp_v * t)^2]
    #: Only used with independent-trajectory XF methods, that is, `decoherence_algo == 5 or 6`
    #: 0 no td width; use the fixed-width Gaussian approximation [ default ]
    #: 1 the td Gaussian width from a free particle Gaussian wave packet, \sigma(t)=\sqrt[\sigma(0)^2 + (wp_v * t)^2]
    #: 2 the Schwarz scheme where the width depends on the instantaneous de Broglie 
    #:   wavelength, \sigma(t)^(-2) = [\sigma(0)^2 * P/ (4 * PI) ]^2
    #: 3 the Subotnik scheme where the width is given as a sum of pairwise widths depending on the auxiliary 
    #:   variables, \sigma_ij(t)^2 = |R_i - R_j| / |P_i - P_j| 
    use_td_width: int = 0

    # ==========================================================================
    # NBRA
    # ==========================================================================

    #: Use NBRA optimization
    #: A flag for NBRA calculations. Since in NBRA, the Hamiltonian is the same for all the trajectories
    #: we can only compute the Hamiltonian related properties once for one trajectory and increase the speed of calculations.
    #: If we set the value to 1 it will consider the NBRA type calculations and other integers the Hamiltonian related properties
    #: are computed for all trajectories.
    #: 0 no NBRA - Hamiltonians for all trajectories are computed explicitly [ default ]
    #: 1: the NBRA is involved - the calculations of the Hamiltonian are conducted for only 1 trajectory, 
    #:    and re-used by all other SH trajectories.  
    isNBRA: int = 0

    # ==========================================================================
    # Entanglement
    # ==========================================================================

    #: Trajectory entanglement option
    #: A selector of a method to couple the trajectories in this ensemble.
    #:  0 no coupling [ default ]
    #:  1 ETHD
    #:  2 ETHD3 (experimental)
    #: 22 another flavor of ETHD3 (experimental)
    #:  3 RPMD
    entanglement_opt: int = 0

    #: ETHD3 coordinate-space exponent
    #: Gaussian exponents that dresses up the trajectories in the ETHD3 method
    #: in the coordinate space, that is   ~exp(-alpha*(R-R0)^2 ) [ default: 1.0 ]
    ETHD3_alpha: float = 1.0

    #: ETHD3 momentum-space exponent
    #: Gaussian exponents that dresses up the trajectories in the ETHD3 method
    #: in the momentum space, that is   ~exp(-beta*(P-P0)^2 ) [ default: 1.0 ]
    ETHD3_beta: float = 1.0

    # ==========================================================================
    # QTAG
    # ==========================================================================

    #: QTAG potential approximation - how to approximate the Hamiltonian matrix 
    #: elements for trajectories that belong to different (or same) surfaces
    #: 0 BAT [default]
    #: 1 LHA
    #: 2 LHAe
    #: 3 BATe
    qtag_pot_approx_method: int = 0

    # ==========================================================================
    # Thermostat, constraints, integration
    # ==========================================================================

    #: Temperature of the system. This parameter could be used even in the NVE simulations
    #: e.g. as a parameters to compute hop acceptance probabilities based on Boltzmann factors [default: 300 K]
    Temperature: float = 300.0

    #: Which ensemble to use in the dynamics.
    #: 0 NVE [default]
    #: 1 NVT
    ensemble: int = 0

    # Thermostat parameters
    thermostat_params: Dict[str, Any] = field(default_factory=dict)

    #: Thermostat DOFs - this list contains the indices of nuclear DOFs which 
    #: shall be coupled to a thermostat directly [ default: [] ]
    thermostat_dofs: List[int] = field(default_factory=list)

    #: Quantum-classical partitioning
    #: This list of integers contains the indices of nuclear DOFs which chall be 
    #: treated "quantum-mechanically", well including with TSH that is. 
    #: These DOFs will determine the velocity-rescaling-based acceptance of the hops,
    #: and these DOFs will be rescaled when the transition is accepted [ default: [0] ]
    quantum_dofs: List[int] = field(default_factory=lambda: [0])

    #: Constrained DOFs
    #: This list of integers contains the indices of the nuclear DOFs to be constrained - their momenta will be constantly 
    #: reset to zero, so the corresponding coordinates will stay fixed [ default: [] ]
    constrained_dofs: List[int] = field(default_factory=list)

    #: Nuclear timestep (a.u.)
    #: the nuclear and electronic integration timesteps [ units: a.u. of time, default: 41.0 a.u. = 1 fs ]
    dt: float = 41.0

    #: Electronic substeps
    #: the number of electronic integration substeps per a nuclear step, such that dt_el = dt_nucl / num_electronic_substeps
    num_electronic_substeps: int = 1

    #: Electronic TDSE integrator selector
    #: for `rep_tdse = 0` (diabatic): 1** - with NBRA
    #: -1 No propagation
    #:  0 Lowdin exp_ with 2-point Hvib_dia 
    #:  1 based on QTAG propagator
    #:  2 based on modified QTAG propagator (Z at two times)
    #:  3 non-Hermitian integrator with 2-point Hvib_dia
    #: 
    #: for `rep_tdse = 1` (adiabatic):  1** - with NBRA
    #: -1 No propagation
    #:  0 LD, with crude splitting,  with exp_  [ default ]
    #:  1 LD, with symmetric splitting, with exp_
    #:  2 LD, original, with exp_
    #:  3 1-point, Hvib integration, with exp_
    #:  4 2-points, Hvib integration, with exp_
    #:  5 3-points, Hvib, integration with the second-point correction of Hvib, with exp_
    #:  6 same as 4, but without projection matrices (T_new = I)
    #: 10 same as 0, but with rotations
    #: 11 same as 1, but with rotations
    #: 12 same as 2, but with rotations
    #: 13 same as 3, but with rotations
    #: 14 same as 4, but with rotations
    #: 15 same as 5, but with rotations
    #:
    #: for `rep_tdse = 2` ( diabatic, density matrix formalism): 1** - with NBRA
    #:  0 mid-point Hvib with the second-point correction of Hvib
    #: 
    #: for `rep_tdse = 3` ( adiabatic, density matrix formalism): 1** - with NBRA
    #:  0 mid-point Hvib with the second-point correction of Hvib
    #:  1 Zhu Liouvillian
    #: 10 same as 0, but with rotations
    electronic_integrator: int = 0

    #: Amplitude transformation option - whether transform the amplitudes by the T transformation matrix
    #: 0 do not transform by the T matrix (naive, but potentially correct approach) 
    #: 1 do transform it (as in LD, but maybe not needed if we directly transform basis)
    ampl_transformation_method: int = 1

    #: Assume trivial crossings always consistent
    #: If set to True (1), we will force the reprojection matrix T_new to be the identity matrix. 
    #: This effectively removes basis-reprojection (local diabatization) approach and turns on the 
    #: "naive" approach where no trivial crossings exist.
    #: 0 No - we do want to use the LD approaches by default. [default]
    #: 1 Yes - one may need to turn on additional state tracking and phase correction methods
    assume_always_consistent: int = 0

    #: Thermal correction to NBRA
    #: Flag setting to use the thermal correction to NBRA 
    #: 0 no [default]
    #: 1 rescale NACs according to TC-NBRA
    thermally_corrected_nbra: int = 0

    # Total energy of the system (a.u.)should be set according to the initial condition of the experiment
    # Used int the TC-NBRA [default: 0.01 a.u.]
    total_energy: float = 0.01

    #: TC-NBRA thermostat frequency
    #: Frequency of the auxiliary thermostats in the TC-NBRA method [default: 0.001]
    tcnbra_nu_therm: float = 0.001

    #: Length of the auxiliary NHC thermostat in the TC-NBRA method [default: 1]
    tcnbra_nhc_size: int = 1

    #: TC-NBRA NAC scaling
    #: Whether to rescale NACs to reflect the fact that the instantaneous velocities 
    #: may be different from those in the ground state sampling. This would rescale off-diagonal
    #: elements of NAC, Hvib, and time-overlap matrices. Does not rescale derivative coupling
    #: vectors cause the chances are - they won't be used in NBRA or if they are used, it is not
    #: NBRA anymore.
    #: 0 do not rescale NACs etc. 
    #: 1 do rescale them [default]
    tcnbra_do_nac_scaling: int = 1

    #: Properties to save (Python list of strings)
    properties_to_save: List[str] = field(default_factory=list)

    #: Bath reorganization energy (Ha) [default: 0.0 Ha]
    reorg_energy: float = 0.0

    # ==========================================================================
    # Utilities
    # ==========================================================================

    def sanity_check(self) -> None:
        """Light consistency checks."""
        if self.enforce_state_following and self.enforced_state_index < 0:
            raise ValueError("enforced_state_index must be >= 0")

    def set_parameters(self, params: Dict[str, Any]) -> None:
        """Update parameters from a dictionary (bp::dict analogue)."""
        for k, v in params.items():
            if not hasattr(self, k):
                pass
                #raise AttributeError(f"Unknown parameter '{k}'")
            setattr(self, k, v)
        self.sanity_check()


