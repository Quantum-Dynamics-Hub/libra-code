# *********************************************************************************
# * Copyright (C) 2019-2024 Alexey V. Akimov
# *
# * This file is distributed under the terms of the GNU General Public License
# * as published by the Free Software Foundation, either version 3 of
# * the License, or (at your option) any later version.
# * See the file LICENSE in the root directory of this distribution
# * or <http://www.gnu.org/licenses/>.
# ***********************************************************************************
"""
.. module:: dynamics
   :platform: Unix, Windows
   :synopsis: This module implements a wrapper function for doing Ehrenfest/TSH/Verlet/etc. dynamics
       List of functions:
           * init_nuclear_dyn_var(Q, P, M, params, rnd)
           * init_electronic_dyn_var(params, isNBRA, rnd)
           * init_amplitudes(q, Cdia, Cadi, dyn_params, compute_model, model_params, transform_direction=0)
           * run_dynamics(_q, _p, _iM, _Cdia, _Cadi, _states, _dyn_params, compute_model, _model_params, rnd)
           * generic_recipe(q, p, iM, _dyn_params, compute_model, _model_params, _init_elec, rnd)
           * run_multiple_sets(init_cond, _dyn_params, compute_model, _model_params, _init_nucl, _init_elec, rnd)

.. moduleauthor:: Alexey V. Akimov

"""

__author__ = "Alexey V. Akimov"
__copyright__ = "Copyright 2019 Alexey V. Akimov"
__credits__ = ["Alexey V. Akimov"]
__license__ = "GNU-3"
__version__ = "1.0"
__maintainer__ = "Alexey V. Akimov"
__email__ = "alexvakimov@gmail.com"
__url__ = "https://quantum-dynamics-hub.github.io/libra/index.html"


import os
import sys
import math
import copy
import time
if sys.platform == "cygwin":
    from cyglibra_core import *
elif sys.platform == "linux" or sys.platform == "linux2":
    from liblibra_core import *

import util.libutil as comn
import libra_py.units as units
import libra_py.data_outs as data_outs
import libra_py.data_savers as data_savers
import libra_py.tsh_stat as tsh_stat
# import libra_py.dynamics as dynamics_io

from . import save


def run_dynamics(dyn_var, _dyn_params, ham, compute_model, _model_params, rnd):
    """

    Args:
        dyn_var ( dynamic_variables ): contains all the dynamical variables (momenta, coordinates, forces,
            quantum amplitudes, active states, projectors, etc.)

        dyn_params ( dictionary ): parameters controlling the execution of the dynamics
            Can contain:

            ///===============================================================================
            ///================= Computing Hamiltonian-related properties ====================
            ///===============================================================================

            * **dyn_params["rep_tdse"]** ( int ): Selects the representation in which nuclear/electronic
                (Ehrenfest core) dynamics is executed. This is the representation chosen as the main
                one to represent the time-dependent wavefunction and the one used to integrate the TD-SE

                - 0: diabatic representation, wfc
                - 1: adiabatic representation, wfc [ default ]
                - 2: diabatic representation, density matrix (e.g. Liouville's picture)
                - 3: adiabatic representation, density matrix (e.g. Liouville's picture)
                - 4: MMST adiabatic representation, wfc

            * **dyn_params["ham_update_method"]** ( int ): How to update Hamiltonian and which type of
                Hamiltonian to update

                - 0: don't update any Hamiltonians
                - 1: recompute only diabatic Hamiltonian [ default ]
                - 2: recompute only adiabatic Hamiltonian

            * **dyn_params["ham_transform_method"]** ( int ): How to transform the Hamiltonians between
                representations

                - 0: don't do any transforms
                - 1: diabatic->adiabatic according to internal diagonalization [ default ]
                - 2: diabatic->adiabatic according to internally stored basis transformation matrix
                - 3: adiabatic->diabatic according to internally stored basis transformation matrix
                - 4: adiabatic->diabatic according to local diabatization method


            * **dyn_params["rep_sh"]** ( int ): selects the representation which is
                used to perform surface hopping

                - 0: diabatic representation
                - 1: adiabatic representation [ default ]


            * **dyn_params["rep_lz"]** ( int ): The representation to compute LZ probabilitieis:

                - 0: diabatic, Eq. 1 of the Belyaev-Lebedev paper [ default ]
                - 1: adiabatic, Eq. 3 of the Belyaev-Lebedev paper, crossing point is determined
                     by the sign change of the diabatic gap
                - 2: adiabatic, Eq. 3 of the Belyaev-Lebedev paper, crossing point is determined
                     by the sign change of the NAC


            * **dyn_params["rep_force"]** ( int ): In which representation to compute forces.
                To clarify - the forces in both representations shall be equivalent, so this flag actually
                selects the type of the properties needed to compute such forces.
                For instance, if it is set to 0, we may be using the derivatives
                of the diabatic Hamiltonians, diabatic states overlaps, etc.

                - 0: diabatic
                - 1: adiabatic [ default ]


            * **dyn_params["force_method"]** ( int ): How to compute forces in the dynamics:

                - 0: don't compute forces at all - e.g. we do not really need them
                - 1: state-specific  as in the TSH or adiabatic (including adiabatic excited states) [ default ]
                - 2: Ehrenfest
                - 3: QTSH force
                - 4: KC-RPMD force


            * **dyn_params["sqc_gamma"]** ( double ): Gamma parameter in the Meyer-Miller Hamiltonian of form:
                H = sum_{i,j} H_{ij} [ 1/(2 hbar) * (q_i - i * p_j) * (q_j + i * p_i) - gamma * delta_ij ]
 
                - 0.0: corresponds to Ehrenfest force [ default ]


            * **dyn_params["enforce_state_following"]** ( int ): Wheather we want to enforce nuclear dynamics
                to be on a given state, regardlenss of the TSH transitions. Note: only matters is `force_method == 1`

                - 0: no [ default ]
                - 1: yes


            * **dyn_params["enforced_state_index"]** (int): If we enforce the nuclear dynamics to be on a given state,
                what is the index of that state [any integer >- 0, default = 0 ]

                The default value of 0 enforces the nuclear dynamics to be on the ground state.
                This is a convenient way of doing NBRA calculations with model systems without the need for
                pre-computing the trajectories


            * **dyn_params["time_overlap_method"]** ( int ): the flag to select how the time-overlaps are computed:

                - 0: don't compute it (perhaps because it was already pre-computed or read)
                - 1: explicitly compute it from the wavefunctions (the Hamiltonian shall have the basis_transform
                     variables updated)  [ default ]


            * **dyn_params["nac_update_method"]** ( int ): How to update NACs and vibronic Hamiltonian before
               electronic TD-SE propagation:

                - 0: don't update them (e.g. for simplest NAC)
                - 1: update according to changed momentum and existing derivative couplings [ default ]
                - 2: update according to time-overlaps (only time-derivative NACs)


            * **dyn_params["nac_algo"]** ( int ): How to compute time-derivative NACs

                - -1: don't, e.g. we use NACs from somewhere else [ default ]
                -  0: use HST formula (if nac_update_method==2)
                -  1: use NPI of Meek and Levine (if nac_update_method==2)


            * **dyn_params["hvib_update_method"]** ( int ): How to update Hvib

                - 0: don't update them (e.g. if it is read externally)
                - 1: update according to regular formula: Hvib = Ham - i * hbar * NAC [ default ]


            * **dyn_params["do_ssy"]** ( int ): Whether to modify the Hamiltonian in the dynamics according
                the Shenvi-Subotnik-Yang (SSY) method, see my Chapter Eq. 3.27 Note, that this is only applied
                in the adiabatic representation

                - 0: don't [ default ]
                - 1: do


            * **dyn_params["do_phase_correction"]** ( int ): the algorithm to correct phases on adiabatic states

                - 0: no phase correction
                - 1: according to our phase correction algorithm [ default ]


            * **dyn_params["phase_correction_tol"]** ( double ): The minimal magnutude of the matrix element
               for which we'll be computing the phase correction. If the overlap is zero, then we don't really
               care about the phase, but if it is not, then this parameter sets out threshold for when we do. [ default: 1e-3 ]


            * **dyn_params["do_nac_phase_correction"]** ( int ):  New phase correction, directly applied to NACs.
               Intended to be used mostly with state_tracking_algo == 4, although can be useful with other state
               tracking algorithms. Should not be used together with `do_phase_correction`

                - 0: no correction [ default ]
                - 1: do this correction


            * **dyn_params["state_tracking_algo"]** ( int ): the algorithm to keep track of the states' identities

                - -1: use LD approach, it includes phase correction too [ default ]
                - 0: no state tracking
                - 1: method of Kosuke Sato (may fail by getting trapped into an infinite loop)
                - 2: Munkres-Kuhn (Hungarian) algorithm 
                - 21: ChatGPT-generated Munkres-Kuhn (Hungarian) algorithm
                - 3: experimental stochastic algorithm, the original version with elimination (known problems)
                - 32: experimental stochastic algorithms with all permutations (too expensive)
                - 33: the improved stochastic algorithm with good scaling and performance, on par with the mincost
                - 4: new, experimental force-based tracking


            * **dyn_params["MK_alpha"]** ( double ): Munkres-Kuhn alpha
                (selects the range of orbitals included in reordering) [default: 0.0]


            * **dyn_params["MK_verbosity"]** ( double ): Munkres-Kuhn verbosity:

                - 0: no extra output [ default ]
                - 1: print details for debugging


            * **dyn_params["convergence"]** ( int ):  A swtich for stochastic reordering algorithm 3 (and variants)
                to choose what happens when an acceptable permutation isn't generated in the set number of attempts:

                - 0: returns the identity permutation (does not require convergence) [ default ]
                - 1: exits and prints an error (requires convergence)


            * **dyn_params["max_number_attempts"]** ( int ): The maximum number of hops that an be
                attempted before either choosing the identity or exiting in stochastic reordering algorithm 3 (and variants).
                Default: 100


            * **dyn_params["min_probability_reordering"]** ( double ): The probability threshold for stochastic
                state reordering algorithm. If a probability for a multi-state stransition is below this value,
                it will be disregarded and set to 0
                The rest of the probabilities will be renormalized
                Default: 0.0


            ///===============================================================================
            ///================= Surface hopping: proposal, acceptance =======================
            ///===============================================================================

            * **dyn_params["tsh_method"]** ( int ): Surface hop proposal methodology:

                - [-1]: adiabatic dynamics, no hops [ default ]
                - 0: Fewest Switches Surface Hopping (FSSH)
                - 1: Global Flux Surface Hopping (GFSH)
                - 2: Markov-State Surface Hopping (MSSH)
                - 3: Landau-Zener (LZ) options
                - 4: Zhu-Nakamura (ZN) options
                - 5: DISH
                - 6: MASH
                - 7: FSSH2
                - 8: FSSH3
                - 9: GFSH (original)


            * **dyn_params["hop_acceptance_algo"]** ( int ): Options to control the acceptance of the proposed hops:

                - 0: accept all proposed hops  [ default ]

                - 10: based on adiabatic energy - accept only those hops that can obey the energy conservation with
                    adiabatic potential energies
                - 11: based on diabatic energy - same as 10, but we use diabatic potential energies

                - 20: based on derivative coupling vectors - accept only those hops that can obey the energy conservation
                    by rescaling nuclear velocities along the directions of derivative couplings for the quantum nuclear DOFs
                - 21: based on difference of state-specific forces - same as 20, but the rescaling is done along the vector
                    parallel to the difference of adiabatic forces on initial and target states

                - 31: accept hops with the probability taken from the quantum Boltzmann distribution
                - 32: accept hops with the probability taken from the classical Maxwell-Boltzmann distribution
                - 33: accept hops with the probability taken from the updated quantum Boltzmann distribution (experimental)

                - 40: based on possibility to conserve energy using tcnbra_ekin variables (for TC-NBRA)


            * **dyn_params["momenta_rescaling_algo"]** ( int ): Options to control nuclear momenta changes upon successful
                or frustrated hops.

                - 0: don't rescale [ default ]

                - 100: based on adiabatic energy, don't reverse on frustrated hops
                - 101: based on adiabatic energy, reverse on frustrated hops
                - 110: based on diabatic energy, don't reverse on frustrated hops
                - 111: based on diabatic energy, reverse on frustrated hops

                - 200: along derivative coupling vectors, don't reverse on frustrated hops
                - 201: along derivative coupling vectors, reverse on frustrated hops
                - 210: along difference of state-specific forces, don't reverse on frustrated hops
                - 211: along difference of state-specific forces, reverse on frustrated hops

                - 40: does not rescale velocities, but rescales  tcnbra_ekin variables


            * **dyn_params["use_Jasper_Truhlar_criterion"]** ( int ): A flag to turn on/off the Jasper-Truhlar criterion
                for reversal of velocities on frustrated hops. According to: Jasper, A. W.; Truhlar, D. G. Chem. Phys. Lett.
                2003, 369, 60− 67

                - 0: don't use this criterion (naive handling)
                - 1: use it [ default ] - the velocities are reversed along the direction d_{a,j} if
                    a) (F_a * d_{a,j}) * (F_j * d_{a,j}) < 0 and b) (v * d_{a,j}) * (F_j * d_{a,j}) < 0
                    where a - is the active state index;  Only in effect, if `momenta_rescaling_algo == 201`


            * **dyn_params["use_boltz_factor"]** ( int ): A flag to scale the proposed hopping probabilities by the
                Boltzmann factor. This is needed for the libra_py/workflows/nbra calculations, where the hop proposal
                probability also includes the factor to account for the hop acceptance probabilities

                - 0: don't scale [ default ]
                - 1: do scale


            //=========== FSSH2 options ==========

            * **dyn_params["fssh2_revision"]** ( int ): Whether to use the revised FSSH2

                - 0: use the FSSH2 as it was formulated in the original paper [ default ]
                - 1: apply the revised version


            //=========== FSSH3 options ==========

            * **dyn_params["fssh3_size_option"]** ( int ): The size of the vectorized density matrix in equations
                to determine hopping probabilites/fluxes

                - 0: N elements - only populations; the matrices are overdetermined
                - 1: N^2 elements - first N elements are populations, then Re and Im parts of upper-triangular coherences
                         that is rho_{0,1}, rho_{0,2}, ... rho_{0,N-1}, rho_{1,2}, ... rho_{1,N-1}, ... rho_{N-2,N-1}  [ default ]


            * **dyn_params["fssh3_approach_option"]** ( int ): The approach to determine the hopping probabilities

                - 0: based on master equation, rho(t+dt) = J * rho(t);  J matrix contains hopping probabilities directly  [ defualt ]
                - 1: based on kinetic approach, drho/dt = J * rho; J matrix contains fluxes


            * **dyn_params["fssh3_decomp_option"]** ( int ): The matrix decomposition method for solving the least-squares problem.
                In the present implementation, is not used

                - 0: bdcSvd
                - 1: fullPivLu
                - 2: fullPivHouseholderQr
                - 3: completeOrthogonalDecomposition [ default ]

            * **dyn_params["fssh3_dt"]** ( double ): The time-step of the optimization procedure in the FSSH3 calculations
                Default: 0.001 a.u.


            * **dyn_params["fssh3_max_steps"]** ( int ): The maximal number of steps in the FSSH3 optimization step
                Default: 1000


            * **dyn_params["fssh3_err_tol"]** ( double ): FSSH3 error tolerance
                Default: 1e-7

            //=========== QTSH options ==========

            * **dyn_params["use_qtsh"]** ( int ): Whether to use QTSH

                - 0: don't apply [ default ]
                - 1: use it

            * **dyn_params["qtsh_force_option"]** ( int ): Nonclassical force options in the adiabatic QTSH

                - 0: Considering only the first-order force, i.e., off-diagonal Ehrenfest force
                - 1: The whole force including the second-order term is used [default]


            //=========== KC-RPMD options ==========

            * **dyn_params["use_kcrpmd"]** ( int ): Whether to use KC-RPMD

                - 0: don't apply [ default ]
                - 1: use it


            ///===============================================================================
            ///================= Decoherence options =========================================
            ///===============================================================================


            * **dyn_params["decoherence_algo"]** ( int ): Selector of the method to incorporate decoherence.

                - [-1]: no decoherence [ default ]
                - 0: SDM and alike
                - 1: instantaneous decoherence options (ID-S, ID-A, ID-C)
                - 2: AFSSH
                - 3: BCSH of Linjun Wang
                - 4: MF-SD of Bedard-Hearn, Larsen, Schwartz
                - 5: SHXF of Min
                - 6: MQCXF
                - 7: DISH, rev2023
                - 8: diabatic IDA, experimental


            * **dyn_params["sdm_norm_tolerance"]** ( double ): Corresponds to the "tol" parameter in the sdm function.
                It controls how much the norm of the old state can be larger than 1.0  before the code stops
                with the error message [ default : 0.0 ]

                Note: only matters if decoherence_algo == 0


            * **dyn_params["dish_decoherence_event_option"]** ( int ):  Selects the how to sample decoherence
                events in the DISH. Note: only matters if dyn_params["tsh_method"] == 5

                - 0: compare the coherence time counter with the decoherence time (simplified DISH)
                - 1: compare the coherence time counter with the time drawn from the exponential distribution
                    with the parameter lambda = 1/decoherence time - this distribution corresponds to
                    the statistics of wait times between the Poisson-distributed events (decoherence)
                    This is what the original DISH meant to do  [ default ]


            * **dyn_params["decoherence_times_type"]** ( int ): Type of dephasing times/rates calculation:

                - -1: set all dephasing rates to zero [ default ]
                - 0: use the rates read out from the input
                - 1: use the energy-based decoherence method (EDC)
                - 2: Schwartz - mean-field Force-based decoherence
                - 3: Schwartz - pair-wise-based decoherences


            * **dyn_params["schwartz_decoherence_inv_alpha"]** ( MATRIX(ndof, 1) ): a matrix of 1/alpha - the parameters
                used in GWP in computing decoherence rates [ default: NULL ]


            * **dyn_params["schwartz_interaction_width"]** ( MATRIX(ndof, 1) ): the parameters for the spatial extent
                of NAC in computing decoherence rates [ default: NULL ]


            * **dyn_params["decoherence_C_param"]** ( double ): An empirical parameter used in the EDC method [ default: 1.0 Ha]


            * **dyn_params["decoherence_eps_param"]** ( double ): An empirical parameter used in the EDC method [ default: 0.1 Ha]


            * **dyn_params["dephasing_informed"]** ( int ): A flag to apply the dephasing-informed approach of Sifain et al.
                to correct dephasing times:

                - 0: don't apply [ default ]
                - 1: use it, will need also the `ave_gaps` input


            * **dyn_params["instantaneous_decoherence_variant"]** ( int ): Option to control the instantaneous
                decoherence methodology, only used with decoherence_algo == 1

                - 0: ID-S
                - 1: ID-A [default] - if the proposed hop is not successful, we project back to the initial state
                                      if the proposed hop is accepted - we project onto that state
                - 2: ID-C - consistent ID - an experimental algorithm
                - 3: ID-A, new: if the proposed hop is not successful, we project out the proposed states
                                if the proposed hop is accepted - we project onto that state
                - 4: ID-F, new: if the proposed hop is not successful, we project out the proposed states
                                but we don't do anything if the hop is successful


            * **dyn_params["collapse_option"]** ( int ): How to collapse wavefunction amplitudes in the decoherence schemes:

                - 0: by rescaling the magnitude of the amplitude vector elements, but preserving "phase" [ default ]
                - 1: by resetting the amplitudes to 1.0+0.0j. This option changes phase


            * **dyn_params["decoherence_rates"]** ( MATRIX(ntates, nstates) ): the matrix of dephasing rates [ units: a.u. of time ^-1 ]


            * **dyn_params["ave_gaps"]** ( MATRIX(ntates, nstates) ): A matrix that contains the averaged
                moduli of the energy gaps: E_ij = <|E_i - E_j|>.  It is needed when dephasing_informed option is used


            * **dyn_params["wp_width"]** ( MATRIX(ndof, 1) ): Widths of all DOFs for Gaussian function as an
                approximation to adiabatic wave packets. According to the choice of the Gaussian width approximation,
                this parameter has different meanings:

                - A constant width in the fixed-width approximation, that is, `use_td_width == 0`
                - The initial width in the free-particle Gaussian wave packet approximation, that is, `use_td_width == 1`
                - The interaction width in the Schwarz scheme, that is, `use_td_width == 2`
                - No influence on the dynamics since the width will be determined by internal variables in the Subotnik scheme,
                    that is, `use_td_width == 3`

                Only used with independent-trajectory XF methods, that is, `decoherence_algo == 5 or 6`


            * **dyn_params["wp_v"]** ( MATRIX(ndof,1) ): The velocity of Gaussian wave packet in the free-particle Gaussian
                approximation, that is, `use_td_width == 1` Only used with independent-trajectory XF methods, that is,
                with `decoherence_algo == 5 or 6`


            * **dyn_params["coherence_threshold"]** ( double ): A population threshold for creating/destroying auxiliary
                trajectories. [ default: 0.01 ]. Only used with independent-trajectory XF methods, that is, with
                `decoherence_algo == 5 or 6`


            * **dyn_params["e_mask"]** ( double ): The masking parameter for computing nabla phase vectors. [ default: 0.0001 Ha ]
                Only used with the MQCXF method, that is, `decoherence_algo == 5`


            * **dyn_params["use_xf_force"]** (int): Whether to use the decoherence force in MQCXF.
                The corresponding electronic propagation is adjusted for the energy conservation.
                Only used with `decoherence_algo == 6`

                - 0: don't use it, so for XF-methods this is only Ehrenfest-like force; EhXF [ default ]
                - 1: The whole force including the XF-correction; MQCXF


            * **dyn_params["project_out_aux"]** (int): Whether to project out the density on an auxiliary
                trajectory when its motion is classically forbidden. Only used with independent-trajectory XF methods,
                that is, `decoherence_algo == 5 or 6`

                - 0: don't [default]
                - 1: do


            * **dyn_params["tp_algo"]** (int): Turning-point algorithm for auxiliary trajectories
                Only used with independent-trajectory XF methods, that is, `decoherence_algo == 5 or 6`

               - 0: no treatment of a turning point
               - 1: collapse to the active state [default]
               - 2: fix auxiliary positions of adiabatic states except for the active state
               - 3: keep auxiliary momenta of adiabatic states except for the active state


            * **dyn_params["use_td_width"]** (int): Whether to use the td Gaussian width for the nuclear
                wave packet approximation. This option can be considered when it comes to unbounded systems.
                This approximation is based on a nuclear wave packet on a free surface:
                \\sigma_x(t)=\\sqrt[\\sigma_x(0)^2 + (wp_v * t)^2]
                Only used with independent-trajectory XF methods, that is, `decoherence_algo == 5 or 6`

                - 0: no td width; use the fixed-width Gaussian approximation [ default ]
                - 1: the td Gaussian width from a free particle Gaussian wave packet, \\sigma(t)=\\sqrt[\\sigma(0)^2 + (wp_v * t)^2]
                - 2: the Schwarz scheme where the width depends on the instantaneous de Broglie
                     wavelength, \\sigma(t)^(-2) = [\\sigma(0)^2 * P/ (4 * PI) ]^2
                - 3: the Subotnik scheme where the width is given as a sum of pairwise widths depending on the auxiliary
                     variables, \\sigma_ij(t)^2 = |R_i - R_j| / |P_i - P_j|

            ///===============================================================================
            ///================= NBRA options =========================================
            ///===============================================================================

            * **dyn_params["isNBRA"]** (int): A flag for NBRA calculations. Since in NBRA, the Hamiltonian is
              the same for all the trajectories we can only compute the Hamiltonian related properties once for one
              trajectory and increase the speed of calculations. If we set the value to 1 it will consider the NBRA
              type calculations and other integers the Hamiltonian related properties are computed for all trajectories.

              - 0: no NBRA - Hamiltonians for all trajectories are computed explicitly [ default ]
              - 1: the NBRA is involved - the calculations of the Hamiltonian are conducted for only 1 trajectory,
                   and re-used by all other SH trajectories.


            ///===============================================================================
            ///================= Entanglement of trajectories ================================
            ///===============================================================================

            * **dyn_params["entanglement_opt"]** ( int ): A selector of a method to couple the trajectories in this ensemble:

                - 0: no coupling across trajectories [ default ]
                - 1: ETHD
                - 2: ETHD3 (experimental)
                - 22: another flavor of ETHD3 (experimental)
                - 3: RPMD


            * **dyn_params["ETHD3_alpha"]** ( double ): Gaussian exponents that dresses up the trajectories in the ETHD3 method
                in the coordinate space, that is   ~exp(-alpha*(R-R0)^2 ) [ default: 0.0 ]


            * **dyn_params["ETHD3_beta"]** ( double ): Gaussian exponents that dresses up the trajectories in the ETHD3 method
                in the coordinate space, that is   ~exp(-alpha*(P-P0)^2 ) [ default: 0.0 ]


            ///===============================================================================
            ///================= QTAG parameters =============================================
            ///===============================================================================

            * **dyn_params["qtag_pot_approx_method"]** ( int ): How to approximate the Hamiltonian matrix elements
              for trajectories that belong to different (or same) surfaces

                - 0 : BAT [ default ]
                - 1 : LHA
                - 2 : LHAe
                - 3 : BATe

            ///===============================================================================
            ///================= Bath, Constraints, and Dynamical controls ===================
            ///===============================================================================

            * **dyn_params["Temperature"]** ( double ): Temperature of the system. This parameter
                could be used even in the NVE simulations e.g. as a parameters to compute hop
                acceptance probabilities based on Boltzmann factors [ default: 300 K]


            * **dyn_params["ensemble"]** ( int ): Which ensemble to use in the dynamics.

                - 0: NVE [ default ]
                - 1: NVT


            * **dyn_params["thermostat_params"]** ( dict ): Parameters controlling the thermostat,
                only relevant for `ensemble = 1`  [ default: {} ]


            * **dyn_params["thermostat_dofs"]** ( list of ints ): Thermostat DOFs
                This list contains the indices of nuclear DOFs which shall be coupled to a thermostat directly.
                [ default: [] ]


            * **dyn_params["quantum_dofs"]** ( list of ints ): Quantum-classical partitioning
                This list of integers contains the indices of nuclear DOFs which chall be treated "quantum-mechanically", well
                including with TSH that is. These DOFs will determine the velocity-rescaling-based acceptance of the hops,
                and these DOFs will be rescaled when the transition is accepted
                [ default: all DOFs ]

                Note: this default is different from the C++ default, where we just use [0]


            * **dyn_params["constrained_dofs"]** ( list of ints ):  Constrained DOFs
                This list of integers contains the indices of the nuclear DOFs to be constrained - their momenta will be constantly
                reset to zero, so the corresponding coordinates will stay fixed
                [ default: [] ]


            * **dyn_params["dt"]** ( double ): the nuclear and electronic integration
                timesteps [ units: a.u. of time, default: 41.0 a.u. = 1 fs ]


            * **dyn_params["num_electronic_substeps"]** ( int ): the number of electronic integration substeps per
                a nuclear step, such that dt_el = dt_nucl / num_electronic_substeps

            * **dyn_params["electronic_integrator"]** ( int ): the method for electronic TD-SE integration:

              rep_tdse = 0 (diabatic): 1** - with NBRA

                -1              - No propagation

                 0              - Lowdin exp_ with 2-point Hvib_dia
                 1              - based on QTAG propagator
                 2              - based on modified QTAG propagator (Z at two times)
                 3              - non-Hermitian integrator with 2-point Hvib_dia

              rep_tdse = 1 (adiabatic):  1** - with NBRA

                -1              -  No propagation

                 0              -  ld, with crude splitting,  with exp_  [ default ]
                 1              -  ld, with symmetric splitting, with exp_
                 2              -  ld, original, with exp_
                 3              -  1-point, Hvib integration, with exp_
                 4              -  2-points, Hvib integration, with exp_
                 5              -  3-points, Hvib, integration with the second-point correction of Hvib, with exp_
                 6              -  same as 4, but without projection matrices (T_new = I)

                10              -  same as 0, but with rotations
                11              -  same as 1, but with rotations
                12              -  same as 2, but with rotations
                13              -  same as 3, but with rotations
                14              -  same as 4, but with rotations
                15              -  same as 5, but with rotations

              rep_tdse = 2 ( diabatic, density matrix formalism): 1** - with NBRA

                 0              -  mid-point Hvib with the second-point correction of Hvib

              rep_tdse = 3 ( adiabatic, density matrix formalism): 1** - with NBRA

                 0              -  mid-point Hvib with the second-point correction of Hvib
                 1              -  Zhu Liouvillian

                10              -  same as 0, but with rotations


            * **dyn_params["ampl_transformation_method"]** (int): Whether transform the amplitudes by the T
              transformation matrix:

              - 0: do not transform by the T matrix (naive, but potentially correct approach)
              - 1: do transform it (as in LD, but maybe not needed if we directly transform basis)


            * **dyn_params["assume_always_consistent"]** (int): If set to True (1), we will force the reprojection matrix
              T_new to be the identity matrix. This effectively removes basis-reprojection (local diabatization) approach
              and turns on the "naive" approach where no trivial crossings exist.

              - 0: No - we do want to use the LD approaches by default. [ default]
              - 1: Yes - one may need to turn on additional state tracking and phase correction methods


            * **dyn_params["thermally_corrected_nbra"]** (int): Flag setting to use the thermal correction to NBRA:

              - 0: No [default]
              - 1: Yes, rescale NACs


            * **dyn_params["total_energy"]** (double): Total energy of the system - should be set according to the
              initial condition of the experiment. Used by the nbra rescaling approach (experimental method) [ units: a.u. ]
              [Default = 0.01 a.u.]


            * **dyn_params["tcnbra_nu_therm"]** (double): Frequency of the auxiliary thermostats in the TC-NBRA method
              [Default = 0.001]


            * **dyn_params["tcnbra_nhc_size"]** (int): Length of the auxiliary NHC thermostat in the TC-NBRA method [units: unitless]
              [Default = 1]


            * **dyn_params["tcnbra_do_nac_scaling"]** (int): Whether to rescale NACs to reflect the fact that the
              instantaneous velocities may be different from those in the ground state sampling. This would rescale off-diagonal
              elements of NAC, Hvib, and time-overlap matrices. Does not rescale derivative coupling vectors cause the chances
              are - they won't be used in NBRA or if they are used, it is not NBRA anymore.

              - 0: - do not rescale NACs etc.
              - 1: - do rescale them [ default ]

            ///===============================================================================
            ///================= Variables specific to Python version: saving ================
            ///===============================================================================


            * **dyn_params["nsteps"]** ( int ): the number of the dynamical steps to do [ default: 1 ]


            * **dyn_params["prefix"]** ( string ): the name of the folder to be created by this function.
                All the data files will be created in that folder [ default: "out" ]


            * **dyn_params["hdf5_output_level"]** ( int ): controls what info to save into HDF5 files (as we run)

                - -1: all returned values are None [ default ]
                - 0: also, all return values are none
                - 1: 1-D info - Relevant variables to request to save: ["timestep", "time", "Ekin_ave", "Epot_ave", "Etot_ave",
                    "Etherm", "E_NHC", "dEkin_ave", "dEpot_ave", "dEtot_ave", "tcnbra_thermostat_energy", "Ekin_ave_qtsh" ]
                - 2: 2-D info - Relevant variables to request to save: ["states", "states_dia", "se_pop_adi", "se_pop_dia",
                    "sh_pop_adi", "sh_pop_dia", "sh_pop_adi_TR", "sh_pop_dia_TR", "mash_pop_adi", "mash_pop_dia", 
                    "fssh3_average_errors"]
                - 3: 3-D info - Relevant variables to request to save: [ "SH_pop", "SH_pop_raw", "D_adi", "D_adi_raw", 
                    "D_dia", "D_dia_raw", "q", "p", "f", "Cadi", "Cdia", "coherence_adi", "coherence_dia", "q_mm", "p_mm",
                    "wp_width", "p_quant", "VP", "f_xf", "qtsh_f_nc" ]
                - 4: 4-D info - Relevant variables to request to save: ["hvib_adi", "hvib_dia", "St", "basis_transform", 
                    "projector", "q_aux", "p_aux", "nab_phase"]
                - 5: 5-D info - Relevat variables to request to save: ["dc1_adi"]


            * **dyn_params["mem_output_level"]** ( int ): controls what info to save into HDF5 files (all at the end)

                Same meaning and output as with hdf5_output_level, except all the variables are first stored in memory
                (while the calcs are running) and then they are written into the HDF5 file at the end of the calculations.
                This is a much faster version of hdf5 saver.


            * **dyn_params["txt_output_level"]** ( int ): controls what info to save into TXT files

                Same meaning and output as with hdf5_output_level, except all the variables are written as text files. [ default: -1 ]


            * **dyn_params["use_compression"]** ( int ): whether to use the compression when writing the info into HDF5 files

                - 0: don't use it [ default ]
                - 1: use the compression, with the level defined by the `compression_level` parameter. However,
                     we found that the compression slows down the processing significantly and sometimes
                     creates even larger files, so it is not recommended to use it

            * **dyn_params["compression_level"]** ( list of 3 ints ): the compression level for ints, doubles, and complex, respectively.
                The larger the value, the higher the compression. [ default: [0, 0, 0] ]


            * **dyn_params["progress_frequency"]** ( double ):  at what intervals print out some "progress" messages. For
                instance, if you have `nsteps = 100` and `progress_frequency = 0.1`, the code will notify you every 10 steps. [ default : 0.1 ]

                Note: this doesn't affect the frequency of the data saving to the file - at this point, we save the information for every
                timestep


            * **dyn_params["properties_to_save"]** ( list of string ): describes what properties to save to the HDF5 files. Note that
                if some properties are not listed in this variable, then they are not saved, even if `mem_output_level` suggests they may be
                saved. You need to BOTH set the appropriate `mem_output_level` AND `properties_to_save`

                [ default:  [ "timestep", "time", "Ekin_ave", "Epot_ave", "Etot_ave",
                           "dEkin_ave", "dEpot_ave", "dEtot_ave", "states", "SH_pop", "SH_pop_raw",
                           "D_adi", "D_adi_raw", "D_dia", "D_dia_raw", "q", "p", "Cadi", "Cdia",
                           "hvib_adi", "hvib_dia", "St", "basis_transform", "projector"
                ]

        ham (nHamiltonian): object controlling Hamiltonian-related calculations

        compute_model ( PyObject ): the pointer to the Python function that performs the Hamiltonian calculations

        _model_params ( dictionary ): contains the selection of a model and the parameters
            for that model Hamiltonian

        rnd ( Random ): random numbers generator object


    Returns:
        None: saves the data to files
    """

    # Create copies of the input dynamical variables, so we could run several such
    # functions with the same input variables without worries that they will be altered
    # inside of each other

    model_params = dict(_model_params)
    dyn_params = dict(_dyn_params)

    nstates = model_params["nstates"]

    # Parameters and dimensions
    critical_params = []
    default_params = {}
    # ================= Computing Hamiltonian-related properties ====================
    default_params.update({"rep_tdse": 1, "ham_update_method": 1, "ham_transform_method": 1,
                           "rep_sh": 1, "rep_lz": 0, "rep_force": 1,
                           "force_method": 1, "enforce_state_following": 0, "enforced_state_index": 0,
                           "time_overlap_method": 0, "nac_update_method": 1, "nac_algo": 0,
                           "hvib_update_method": 1, "do_ssy": 0,
                           "do_phase_correction": 1, "phase_correction_tol": 1e-3,
                           "do_nac_phase_correction": 0,
                           "state_tracking_algo": -1, "MK_alpha": 0.0, "MK_verbosity": 0,
                           "convergence": 0, "max_number_attempts": 100, "min_probability_reordering": 0.0,
                           "is_nbra": 0, "icond": 0, "nfiles": -1, "thermally_corrected_nbra": 0, "total_energy": 0.01,
                           "tcnbra_nu_therm": 0.001, "tcnbra_nhc_size": 1, "tcnbra_do_nac_scaling": 1
                           })

    # ================= Surface hopping: proposal, acceptance =======================
    default_params.update({"tsh_method": -1, "hop_acceptance_algo": 0, "momenta_rescaling_algo": 0,
                           "use_Jasper_Truhlar_criterion": 1,
                           "use_boltz_factor": 0
                           })

    # ================= FSSH2 specific ====================
    default_params.update({"fssh2_revision": 0})

    # ================= FSSH3 specific ====================
    default_params.update({"fssh3_size_option": 1, "fssh3_approach_option": 0, "fssh3_decomp_option": 3,
                           "fssh3_dt": 0.001, "fssh3_max_steps": 1000, "fssh3_err_tol": 1e-7
                           })

    # ================= QTSH specific ====================
    default_params.update({"use_qtsh": 0, "qtsh_force_option": 1})

    # ================= Decoherence options =========================================
    default_params.update(
        {
            "decoherence_algo": -1,
            "sdm_norm_tolerance": 0.0,
            "dish_decoherence_event_option": 1,
            "decoherence_times_type": -1,
            "decoherence_C_param": 1.0,
            "decoherence_eps_param": 0.1,
            "dephasing_informed": 0,
            "instantaneous_decoherence_variant": 1,
            "collapse_option": 0,
            "decoherence_rates": MATRIX(
                nstates,
                nstates),
            "ave_gaps": MATRIX(
                nstates,
                nstates),
            "schwartz_decoherence_inv_alpha": MATRIX(
                dyn_var.ndof,
                1),
            "schwartz_interaction_width": MATRIX(
                dyn_var.ndof,
                1),
            "wp_width": MATRIX(
                dyn_var.ndof,
                1),
            "wp_v": MATRIX(
                dyn_var.ndof,
                1),
            "coherence_threshold": 0.01,
            "e_mask": 0.0001,
            "use_xf_force": 0,
            "project_out_aux": 1,
            "tp_algo": 1,
            "use_td_width": 0})

    # ================= Entanglement of trajectories ================================
    default_params.update({"entanglement_opt": 0, "ETHD3_alpha": 0.0, "ETHD3_beta": 0.0})

    # ================= Bath, Constraints, and Dynamical controls ===================
    default_params.update({"Temperature": 300.0, "ensemble": 0, "thermostat_params": {},
                           "quantum_dofs": None, "thermostat_dofs": [], "constrained_dofs": [],
                           "dt": 1.0 * units.fs2au, "num_electronic_substeps": 1,
                           "electronic_integrator": 0
                           })

    # ================= Variables specific to Python version: saving ================
    default_params.update({"nsteps": 1,
                           "prefix": "out",
                           "prefix2": "out2",
                           "hdf5_output_level": -1,
                           "mem_output_level": -1,
                           "txt_output_level": -1,
                           "txt2_output_level": -1,
                           "use_compression": 0,
                           "compression_level": [0,
                                                 0,
                                                 0],
                           "progress_frequency": 0.1,
                           "properties_to_save": ["timestep",
                                                  "time",
                                                  "Ekin_ave",
                                                  "Epot_ave",
                                                  "Etot_ave",
                                                  "dEkin_ave",
                                                  "dEpot_ave",
                                                  "dEtot_ave",
                                                  "states",
                                                  "SH_pop",
                                                  "SH_pop_raw",
                                                  "se_pop_adi",
                                                  "se_pop_dia",
                                                  "sh_pop_adi",
                                                  "D_adi",
                                                  "D_adi_raw",
                                                  "D_dia",
                                                  "D_dia_raw",
                                                  "q",
                                                  "p",
                                                  "f",
                                                  "Cadi",
                                                  "Cdia",
                                                  "q_mm",
                                                  "p_mm",
                                                  "hvib_adi",
                                                  "hvib_dia",
                                                  "St",
                                                  "basis_transform",
                                                  "projector"]})

    comn.check_input(dyn_params, default_params, critical_params)

    prefix = dyn_params["prefix"]
    prefix2 = dyn_params["prefix2"]
    rep_tdse = dyn_params["rep_tdse"]
    # rep_ham = dyn_params["rep_ham"]
    nsteps = dyn_params["nsteps"]
    dt = dyn_params["dt"]
    phase_correction_tol = dyn_params["phase_correction_tol"]
    hdf5_output_level = dyn_params["hdf5_output_level"]
    mem_output_level = dyn_params["mem_output_level"]
    txt_output_level = dyn_params["txt_output_level"]
    txt2_output_level = dyn_params["txt2_output_level"]
    do_phase_correction = dyn_params["do_phase_correction"]
    state_tracking_algo = dyn_params["state_tracking_algo"]
    force_method = dyn_params["force_method"]
    properties_to_save = dyn_params["properties_to_save"]
    use_compression = dyn_params["use_compression"]
    compression_level = dyn_params["compression_level"]
    ensemble = dyn_params["ensemble"]
    time_overlap_method = dyn_params["time_overlap_method"]
    decoherence_algo = dyn_params["decoherence_algo"]
    is_nbra = dyn_params["is_nbra"]
    icond = dyn_params["icond"]    # Setting the initial geoemtry in the dynamics
    nfiles = dyn_params["nfiles"]  # The number of loaded Ham files
    tsh_method = dyn_params["tsh_method"]
    thermally_corrected_nbra = dyn_params["thermally_corrected_nbra"]
    total_energy = dyn_params["total_energy"]
    use_qtsh = dyn_params["use_qtsh"]
    qtsh_force_option = dyn_params["qtsh_force_option"]

    states = intList()

    ndia = dyn_var.ndia  # Cdia.num_of_rows
    nadi = dyn_var.nadi  # Cadi.num_of_rows
    nnucl = dyn_var.ndof  # q.num_of_rows
    ntraj = dyn_var.ntraj  # q.num_of_cols

    if (dyn_params["quantum_dofs"] is None):
        dyn_params["quantum_dofs"] = list(range(nnucl))

    # Initialize savers
    _savers = save.init_tsh_savers(dyn_params, model_params, nsteps, ntraj, nnucl, nadi, ndia)

    # Open and close the output files for further writing
    if _savers["txt_saver"] is not None:
        _savers["txt_saver"].save_data_txt(F"{prefix}", properties_to_save, "w", 0)

    if _savers["txt2_saver"] is not None:
        _savers["txt2_saver"].save_data_txt(F"{prefix2}", properties_to_save, "w", 0)

    model_params.update({"timestep": icond})
    update_Hamiltonian_variables(dyn_params, dyn_var, ham, ham, compute_model, model_params, 0)
    update_Hamiltonian_variables(dyn_params, dyn_var, ham, ham, compute_model, model_params, 1)

    therm = ThermostatList()
    if ensemble == 1:
        for traj in range(ntraj):
            therm.append(Thermostat(dyn_params["thermostat_params"]))
            therm[traj].set_Nf_t(len(dyn_params["thermostat_dofs"]))
            therm[traj].init_nhc()

    if decoherence_algo == 2:
        dyn_var.allocate_afssh()
    elif decoherence_algo == 3:  # BCSH
        dyn_var.allocate_bcsh()
    elif decoherence_algo == 5:  # SHXF
        dyn_var.allocate_shxf()
    elif decoherence_algo == 6:  # MQCXF
        dyn_var.allocate_mqcxf()
    elif decoherence_algo == 9:  # simple decoherence
        dyn_var.allocate_simple_decoherence()

    if tsh_method == 5 or decoherence_algo == 7:  # DISH or DISH_rev2023
        dyn_var.allocate_dish()
    if tsh_method == 7 or tsh_method == 8 or tsh_method == 9:  # FSSH2 or FSSH3 or original GFSH
        dyn_var.allocate_fssh2()
        dyn_var.save_curr_dm_into_prev()

    if thermally_corrected_nbra == 1:
        dyn_var.allocate_tcnbra()

    if use_qtsh == 1:  # QTSH
        dyn_var.allocate_qtsh()

    ham_aux = nHamiltonian(ham)

    # Do the propagation
    for i in range(nsteps):
        # Energies
        Ekin, Epot, Etot, dEkin, dEpot, dEtot = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        Etherm, E_NHC = 0.0, 0.0
        save.save_tsh_data_1234_new(_savers, dyn_params, i, dyn_var, ham)

        # ============ Propagate ===========
        index = i + icond
        if nfiles > 0:
            index = index % nfiles
        model_params.update({"timestep": index})

        compute_dynamics(dyn_var, dyn_params, ham, ham_aux, compute_model, model_params, rnd, therm)

        if _savers["txt_saver"] is not None:
            _savers["txt_saver"].save_data_txt(F"{prefix}", properties_to_save, "a", i)

        if _savers["txt2_saver"] is not None:
            _savers["txt2_saver"].save_data_txt(F"{prefix2}", properties_to_save, "a", 0)

    if _savers["mem_saver"] is not None:
        _savers["mem_saver"].save_data(F"{prefix}/mem_data.hdf", properties_to_save, "w")
        return _savers["mem_saver"]


# def generic_recipe(q, p, iM, _dyn_params, compute_model, _model_params, _init_elec, rnd):
def generic_recipe(_dyn_params, compute_model, _model_params, _init_elec, _init_nucl, rnd):
    """
    This function initializes electronic variables and Hamiltonians for a given set of
    nuclear variables (those are needed, not initialized here) and then runs the calculations
    on multiple trajectories

    Args:
        _dyn_params ( dictionary ): control parameters for running the dynamics
            see the documentation of the :func:`libra_py.dynamics.run_dynamics` function

        compute_model ( PyObject ): the pointer to the Python function that performs the Hamiltonian calculations

        _model_params ( dictionary ): contains the selection of a model and the parameters
            for that model Hamiltonian. In addition to all the parameters, should contain the
            key

            * **_model_params["model0"]** ( int ): the selection of the model to be used to compute diabatic-to-adiabatic
                transformation. The particular value chosen here depends on the way the `compute_model` functions is defined

               :Note: the function selected should be able to generate the transformation matrix.

        _init_elec ( dictionary ): control parameters for initialization of electronic variables
            see the documentation of the :func:`libra_py.dynamics.init_electronic_dyn_var` function

        _init_nucl (dict): controls how to sample nuclear DOFs from what the user provides. Can contain:
          * **_init_nucl["init_type"]** (int) : the type of the sampling:
            - 0 : Coords and momenta are set exactly to the given value
            - 1 : Coords are set, momenta are sampled
            - 2 : Coords are sampled, momenta are set
            - 3 : Both coords and momenta are samples [ default ]
          * **_init_nucl["ndof"]** (int) : the number of nuclear DOFs [default: 1]
          * **_init_nucl["q"]** ( list of `ndof` doubles ): average coordinate of all trajectories, for each nuclear DOF,
                  a.u. of lenth [ default: [0.0] ]
          * **_init_nucl["p"]** ( list of `ndof` doubles ): average momentum of all trajectories, for each nuclear DOF,
                  a.u. of lenth [ default: [0.0] ]
          * **_init_nucl["mass"]** (list of `ndof` doubles): the mass for each nuclear DOF, in a.u. of mass [ default: [2000.0] ]
          * **_init_nucl["force_constant"]** (list of `ndof` doubles) : the force constant of the harmonic oscillator in each
                  nuclear DOF, determins the width of the coordinate and momentum samplings [ default: [0.001] ]

        rnd ( Random ): random numbers generator object - controls the initialization of nuclear DOFs

    Returns:
        tuple : tuple returned by the :func:`libra_py.dynamics.run_dynamics` function



    """

    model_params = dict(_model_params)
    dyn_params = dict(_dyn_params)
    init_elec = dict(_init_elec)
    init_nucl = dict(_init_nucl)

    # ========== Model params ==========
    comn.check_input(model_params, {}, ["model0"])

    # ========== Nuclear params ==========
    # First check the inputs for the initialization of nuclear variables
    comn.check_input(init_nucl, {"init_type": 3, "ndof": 1, "q": [0.0], "p": [0.0], "mass": [2000.0],
                                 "force_constant": [0.01], "q_width": [1.0], "p_width": [1.0]}, [])
    ndof = len(init_nucl["mass"])

    # ========== Dynamics params ==========
    # Now, we can use the information of ndof to define the default quantum_dofs
    comn.check_input(dyn_params, {"rep_tdse": 1, "rep_sh": 1, "is_nbra": 0, "direct_init": 0, "ntraj": 1,
                                  "quantum_dofs": list(range(ndof))
                                  }, [])

    # Extract info from the updated dyn_params:
    is_nbra = dyn_params["isNBRA"]
    ntraj = dyn_params["ntraj"]
    direct_init = dyn_params["direct_init"]

    # ========== Electronic params =========
    # Setup the electronic variables
    init_elec.update({"is_nbra": is_nbra})
    comn.check_input(init_elec, {"ndia": 1, "nadi": 1}, [])

    # Extract information from the updated electronic variables
    ndia = init_elec["ndia"]
    nadi = init_elec["nadi"]

    # sys.exit(0)

    # Setup the dynamical variables object
    dyn_var = dyn_variables(ndia, nadi, ndof, ntraj)

    # Initialize nuclear variables
    dyn_var.init_nuclear_dyn_var(init_nucl, rnd)

    # print("Initial coordinates")
    # dyn_var.get_coords().show_matrix()
    # sys.exit(0)

    # Initialize electronic variables
    dyn_var.init_amplitudes(init_elec, rnd)
    dyn_var.init_density_matrix(init_elec)
    # dyn_var.init_active_states(init_elec, rnd)

    # Setup the hierarchy of Hamiltonians
    ham = nHamiltonian(ndia, nadi, ndof)
    if is_nbra == 1:
        ham.add_new_children(ndia, nadi, ndof, 1)
    else:
        ham.add_new_children(ndia, nadi, ndof, ntraj)
    ham.init_all(2, 1)

    # Compute internals of the Hamiltonian objects
    model_params1 = dict(model_params)
    model_params1.update({"model": model_params["model0"], "timestep": 0})

    # We set up this temporary dict such that we we request diabatic Ham calculations
    # followed by the transformation to the adiabatic one - this is what we need to get
    # the transformation matrices to convert amplitudes between the representations
    dyn_params1 = dict(dyn_params)

    # sys.exit(0)
    if (dyn_params["ham_update_method"] == 2):
        pass
        # update_Hamiltonian_variables( dyn_params1, dyn_var, ham, ham, compute_model, model_params1, 0)
        # update_Hamiltonian_variables( dyn_params1, dyn_var, ham, ham, compute_model, model_params1, 1)
    else:
        dyn_params1.update({"ham_update_method": 1, "ham_transform_method": 1})
        # update_Hamiltonian_q( dyn_params1, dyn_var, ham, compute_model, model_params1)
        # print( dyn_params1 )
        # print(model_params1)
        # sys.exit()
        update_Hamiltonian_variables(dyn_params1, dyn_var, ham, ham, compute_model, model_params1, 0)
        # sys.exit(0)
        update_Hamiltonian_variables(dyn_params1, dyn_var, ham, ham, compute_model, model_params1, 1)

    # sys.exit(0)

    # Update internal dynamical variables using the computed properties of the Hamiltonian objects
    # Set up the "rep_tdse" variable here to the representation that coinsides with the initial representation
    # of electronic variables - this will convert the amplitudes to the proper representation
    dyn_var.update_basis_transform(ham)
    dyn_var.update_amplitudes({"rep_tdse": init_elec["rep"]}, ham)
    dyn_var.update_density_matrix(dyn_params, ham, 1)

    if dyn_params["rep_sh"] == 1:
        dyn_var.init_active_states(init_elec, rnd)
    elif dyn_params["rep_sh"] == 0:
        dyn_var.init_active_states_dia(init_elec, rnd)

    # print("Initial adiabatic amplitudes")
    # dyn_var.get_ampl_adi().show_matrix()

    # print("Initial diabatic amplitudes")
    # dyn_var.get_ampl_dia().show_matrix()

    # print("Initial adiabatic DM")
    # dyn_var.get_dm_adi(0).show_matrix()

    # print("Initial diabatic DM")
    # dyn_var.get_dm_dia(0).show_matrix()

    if dyn_params["rep_sh"] == 1:
        print("Active states (adiabatic)")
        print(Cpp2Py(dyn_var.act_states))

        print("Initial adiabatic populations")
        pops_sh1 = dyn_var.compute_average_sh_pop(1)
        print(Cpp2Py(pops_sh1))

    elif dyn_params["rep_sh"] == 0:
        print("Active states (diabatic)")
        print(Cpp2Py(dyn_var.act_states_dia))

        print("Initial diabatic populations")
        pops_sh0 = dyn_var.compute_average_sh_pop(0)
        print(Cpp2Py(pops_sh0))

    # Finally, start the dynamics calculations
    res = run_dynamics(dyn_var, dyn_params, ham, compute_model, model_params, rnd)
    return res


def run_multiple_sets(init_cond, _dyn_params, compute_model, _model_params, _init_nucl, _init_elec, rnd):
    """
    This function runs a series of NA-MD calculations for (potentially) different initial
    conditions and returns an object that contains results for each of the
    runs, and/or just stores the results into the corresponding folders.

    Args:

      init_cond ( list of 3-element lists ): mean coordinates and momenta as well as the masses of all DOFs

          The format is this:
          init_cond =
            init_cond[0] = [q0, p0, M0 ]  <- 0-th set of initial variables
            init_cond[1] = [q1, p1, M1 ]  <- 1-st set of initial variables
            ...                           etc.
            init_cond[N] = [qN, pN, MN ]

          Here, q0 = [q00, q01, ... q0<ndof>], same for p0 and M0

      _dyn_params ( list of dictionaries ): control parameters for running the dynamics.
          Each item of the list is a dictionary for a particular set of initial conditions.
          see the documentation of the :func:`libra_py.dynamics.run_dynamics` function

      compute_model ( list of PyObject objects ): the pointers to the Python functions
          that performs the Hamiltonian calculations, one per each initial conditions set

      _model_params ( list of dictionaries ): contains the selection of a model and the parameters
          for that model Hamiltonian. Each item corresponds to a particular simulation set.
          In addition to all the parameters, should contain the key:

          * **_model_params["model0"]** ( int ): the selection of the model to be used to compute diabatic-to-adiabatic
          transformation. The particular value chosen here depends on the way the `compute_model` functions is defined

              :Note: the function selected should be able to generate the transformation matrix.

      _init_nucl ( list of dictionary ): control parameters for initialization of nuclear variables
          One item per each initial conditions set.
          see the documentation of the :func:`libra_py.dynamics.init_nuclear_dyn_var` function

      _init_elec ( dictionary ): control parameters for initialization of electronic variables.
          One item per each initial conditions set.
          see the documentation of the :func:`libra_py.dynamics.init_electronic_dyn_var` function

      rnd ( Random ): random numbers generator object

    Returns:
        list : each element of the list is a tuple returned by the
            :func:`libra_py.dynamics.run_dynamics` function
            Each such element represents the results for a particular set of calculations

    """

    res = []
    for icond_indx, icond in enumerate(init_cond):

        output_prefix = _dyn_params[icond_indx]["prefix"]
        dyn_params = dict(_dyn_params[icond_indx])

        q0, p0, M0 = icond[0], icond[1], icond[2]  # each one is a ndof-elements list
        q, p, iM = init_nuclear_dyn_var(q0, p0, M0, _init_nucl[icond_indx], rnd)

        dyn_params.update({"prefix": F"{output_prefix}_{icond_indx}"})

        res_i = generic_recipe(
            q,
            p,
            iM,
            dyn_params,
            compute_model[icond_indx],
            _model_params[icond_indx],
            _init_elec[icond_indx],
            rnd)
        res.append(res_i)

    return res
