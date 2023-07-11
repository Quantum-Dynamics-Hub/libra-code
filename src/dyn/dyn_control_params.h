/*********************************************************************************
* Copyright (C) 2019 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file dyn_control_params.h
  \brief The file implements a class to store the control parameters
*/


#ifndef DYN_CONTROL_PARAMS_H
#define DYN_CONTROL_PARAMS_H

#if defined(USING_PCH)
#include "../pch.h"
#else
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#endif

#include "../math_linalg/liblinalg.h"
#include "../Units.h"


/// liblibra namespace
namespace liblibra{


using namespace liblinalg;

/// libdyn namespace
namespace libdyn{

namespace bp = boost::python;



class dyn_control_params{

  public:

  ///===============================================================================
  ///================= Computing Hamiltonian-related properties ====================
  ///===============================================================================
  /**
    Selects the representation in which nuclear/electronic (Ehrenfest core) dynamics 
    is executed. This is the representation chosen as the main one to represent the 
    time-dependent wavefunction and the one used to integrate the TD-SE
    
    Options:
     - 0: diabatic representation, wfc
     - 1: adiabatic representation, wfc [ default ]
     - 2: diabatic representation, density matrix (e.g. Liouville's picture)
     - 3: adiabatic representation, density matrix (e.g. Liouville's picture)
  */
  int rep_tdse;


  /** 
   The representation of the Hamiltonian update. This is the representation in which the 
   computed properties are assumed to be. For instance, we may have set it to 1, to read the 
   adiabatic energies and nonadiabatic couplings, to bypass the diabatic-to-adiabatic transformation,
   which may be useful in some atomistic calculations, or with the NBRA

   Options:
     - 0: diabatic representation [ default ]
     - 1: adiabatic representation 

  */
  //int rep_ham;  /// TO BE REPLACED BY the ham_update_method


  /** 
   How to update Hamiltonian and which type of Hamiltonian to update

   Options:
     - 0: don't update any Hamiltonians
     - 1: recompute only diabatic Hamiltonian [ default ]
     - 2: recompute only adiabatic Hamiltonian
  */
  int ham_update_method;


  /** 
   How to transform the Hamiltonians between representations

   Options:
     - 0: don't do any transforms
     - 1: diabatic->adiabatic according to internal diagonalization [ default ]
     - 2: diabatic->adiabatic according to internally stored basis transformation matrix
     - 3: adiabatic->diabatic according to internally stored basis transformation matrix
     - 4: adiabatic->diabatic according to local diabatization method

  */
  int ham_transform_method;


  /**
   The representation to run the SH.

   Options: 
     - 0: diabatic
     - 1: adiabatic [ default ]
  */
  int rep_sh;


  /** 
    The representation to compute LZ probabilitieis.
 
    Options:
      - 0: diabatic, Eq. 1 of the Belyaev-Lebedev paper [ default ]
      - 1: adiabatic, Eq. 3 of the Belyaev-Lebedev paper, crossing point is determined
           by the sign change of the diabatic gap
      - 2: adiabatic, Eq. 3 of the Belyaev-Lebedev paper, crossing point is determined
           by the sign change of the NAC
  */
  int rep_lz;


  /** 
    In which representation to compute forces. To clarify - the forces
    in both representations shall be equivalent, so this flag actually
    selects the type of the properties needed to compute such forces.
    For instance, if it is set to 0, we may be using the derivatives
    of the diabatic Hamiltonians, diabatic states overlaps, etc.

    Options:
      - 0: diabatic
      - 1: adiabatic [ default ]
  */
  int rep_force;


  /** 
    How to compute forces in the dynamics.
 
    Options:
      - 0: don't compute forces at all - e.g. we do not really need them
      - 1: state-specific  as in the TSH or adiabatic (including adiabatic excited states) [ default ]
      - 2: Ehrenfest
  */
  int force_method;


  /** 
    Wheather we want to enforce nuclear dynamics to be on a given state, regardlenss of the TSH transitions
 
    Options:
      - 0: no [ default ]
      - 1: yes

    Note: only matters is `force_method == 1`
  */
  int enforce_state_following; 

  /** 
    If we enforce the nuclear dynamics to be on a given state, what is the index of that state [any integer >- 0, default = 0 ]
 
    The default value of 0 enforces the nuclear dynamics to be on the ground state. This is a convenient way of doing NBRA calculations
    with model systems without the need for pre-computing the trajectories 

  */
  int enforced_state_index;  


  /**
    How do get the time-overlaps in the dynamics.

    Options:
      - 0: don't compute it (perhaps because it was already pre-computed or read)
      - 1: explicitly compute it from the wavefunctions (the Hamiltonian shall have the basis_transform variables updated)  [ default ]
  */
  int time_overlap_method;


  /** 
    How to update NACs 

    Options:
      - 0: don't update them (e.g. for simplest NAC)
      - 1: update according to changed momentum and existing derivative couplings [ default ]
      - 2: update according to time-overlaps (only time-derivative NACs)
  */
  int nac_update_method;


  /**
    How to compute time-derivative NACs
 
    Options:
      - -1: don't, e.g. we use NACs from somewhere else [ default ]
      -  0: use HST formula (if nac_update_method==2)
      -  1: use NPI of Meek and Levine (if nac_update_method==2)
  */
  int nac_algo;


  /** 
    How to update Hvib 

    Options:
      - 0: don't update them (e.g. if it is read externally)
      - 1: update according to regular formula: Hvib = Ham - i * hbar * NAC [ default ]
  */
  int hvib_update_method;


  /** 
    Whether to modify the Hamiltonian in the dynamics according the Shenvi-Subotnik-Yang (SSY)
    method, see my Chapter Eq. 3.27
    Note, that this is only applied in the adiabatic representation

    Options:
      - 0: don't [ default ]
      - 1: do
  */
  int do_ssy;


  /** 
    The algorithm to correct phases on adiabatic states

    Options: 
      - 0: no phase correction
      - 1: according to our phase correction algorithm [ default ]
  */
  int do_phase_correction;


  /** 
    The minimal magnutude of the matrix element for which we'll be computing the phase correction
    If the overlap is zero, then we don't really care about the phase, but if it is not, then this
    parameter sets out threshold for when we do.  [ default: 1e-3 ]
  */
  double phase_correction_tol;


  /** 
    State tracking algorithm:
      - 0: no state tracking
      - 1: method of Kosuke Sato (may fail by getting trapped into an infinite loop)
      - 2: Munkres-Kuhn (Hungarian) algorithm [ default ]
      - 3: experimental stochastic algorithm, the original version with elimination (known problems)
      - 32: experimental stochastic algorithms with all permutations (too expensive)
      - 33: the improved stochastic algorithm with good scaling and performance, on par with the mincost


  */
  int state_tracking_algo;


  /** 
    Munkres-Kuhn alpha (selects the range of orbitals included in reordering) [default: 0.0]
  */
  double MK_alpha;

  /**
    Munkres-Kuhn verbosity.
 
    Options:
      - 0: no extra output [ default ]
      - 1: prints extra details on what the algorithm is doing, for debugging
  */
  int MK_verbosity;

  /**
    A swtich for stochastic reordering algorithm 3 to choose what happens when an acceptable permutation isn't generated in the set number of attempts:
                0: returns the identity permutation (does not require convergence)
                1: exits and prints an error (requires convergence)
  */
  int convergence;
  
  /**
  The maximum number of hops that an be attempted before either choosing the identity or exiting in stochastic reordering algorithm 3. 
  */
  int max_number_attempts;

  /**
  The probability threshold for stochastic state reordering algorithm. 
  If a probability for a multi-state stransition is below this value, it will be disregarded and set to 0
  The rest of the probabilities will be renormalized
  Default: 0.0 
  */
  double min_probability_reordering; 


  ///===============================================================================
  ///================= Surface hopping: proposal, acceptance =======================
  ///===============================================================================

  /** 
    Surface hop proposal methodology.

    Options: 
      - [-1]: adiabatic dynamics, no hops [ default ]
      - 0: Fewest Switches Surface Hopping (FSSH)
      - 1: Global Flux Surface Hopping (GFSH)
      - 2: Markov-State Surface Hopping (MSSH)
      - 3: Landau-Zener (LZ) options
      - 4: Zhu-Nakamura (ZN) options
      - 5: DISH
      - 6: MASH
      - 7: FSSH2
  */
  int tsh_method;


  /**
    Options to control the acceptance of the proposed hops.

    Options:
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
  */
  int hop_acceptance_algo;


  /**
    Options to control nuclear momenta changes upon successful or frustrated hops.

    Options:

      - 0: don't rescale [ default ]

      - 100: based on adiabatic energy, don't reverse on frustrated hops
      - 101: based on adiabatic energy, reverse on frustrated hops
      - 110: based on diabatic energy, don't reverse on frustrated hops
      - 111: based on diabatic energy, reverse on frustrated hops

      - 200: along derivative coupling vectors, don't reverse on frustrated hops
      - 201: along derivative coupling vectors, reverse on frustrated hops
      - 210: along difference of state-specific forces, don't reverse on frustrated hops
      - 211: along difference of state-specific forces, reverse on frustrated hops
  */
  int momenta_rescaling_algo;

  /**
    A flag to scale the proposed hopping probabilities by the
    Boltzmann factor. This is needed for the libra_py/workflows/nbra calculations, where the hop proposal 
    probability also includes the factor to account for the hop acceptance probabilities

    Options:

      - 0: don't scale [ default ]
      - 1: do scale
  */
   int use_boltz_factor;


  
  ///===============================================================================
  ///================= Decoherence options =========================================
  ///===============================================================================

  /**
    Selector of the method to incorporate decoherence.

    Options:

      - [-1]: no decoherence [ default ]
      - 0: SDM and alike
      - 1: instantaneous decoherence options (ID-S, ID-A, ID-C)
      - 2: AFSSH
      - 3: BCSH of Linjun Wang
      - 4: MF-SD of Bedard-Hearn, Larsen, Schwartz
  */
  double decoherence_algo;


  /**
    Corresponds to the "tol" parameter in the sdm function. It controls 
    how much the norm of the old state can be larger than 1.0  before the 
    code stops with the error message [ default: 0.0 ]

    Note: only matters if decoherence_algo == 0
  **/
  double sdm_norm_tolerance;


  /**
    Selects the how to sample decoherence events in the DISH.
    Possible options:
      - 0: compare the coherence time counter with the decoherence time (simplified DISH) 
      - 1: compare the coherence time counter with the time drawn from the exponential distribution
           with the parameter lambda = 1/decoherence time - this distribution corresponds to 
           the statistics of wait times between the Poisson-distributed events (decoherence)
           This is what the original DISH meant to do [ default ]

    Note: only matters if tsh_method == 3
  **/
  int dish_decoherence_event_option;


  /**
    Type of dephasing times/rates calculation:

      - -1: set all dephasing rates to zero [ default ]
      - 0: use the rates read out from the input 
      - 1: use the energy-based decoherence method (EDC)    
      - 2: Schwartz - mean-field Force-based decoherence
      - 3: Schwartz - pair-wise-based decoherences
  */  
  int decoherence_times_type;


  /**
    MATRIX(ndof, 1) of 1/alpha - the parameters used in GWP in
    computing decoherence rates [ default: NULL ]
  */
  MATRIX* schwartz_decoherence_inv_alpha;


  /**
    An empirical parameter used in the EDC method: [ default = 1.0 Ha ]
  */ 
  double decoherence_C_param;


  /**
    An empirical parameter used in the EDC method: [ default = 0.1 Ha ]
  */ 
  double decoherence_eps_param;


  /**
    A flag to apply the dephasing-informed approach of Sifain et al 
    to correct dephasing times: 

      - 0: don't apply [ default ]
      - 1: use it 
  */
  int dephasing_informed; 


  /**
    Option to control the instantaneous decoherence methodology,
    only used with decoherence_algo == 1

      - 0: ID-S
      - 1: ID-A [default]
      - 2: ID-C - consistent ID - an experimental algorithm
  */
  int instantaneous_decoherence_variant;


  /**
    How to collapse wavefunction amplitudes in the decoherence schemes:
      - 0: by rescaling the magnitude of the amplitude vector elements, but preserving "phase" [ default ]
      - 1: by resetting the amplitudes to 1.0+0.0j. This option changes phase 

  */
  int collapse_option;


  /**
    Dephasing rates provided by user
    [ default : NULL ]
  */
  MATRIX* decoherence_rates;


  /**
    A matrix that contains the averaged moduli of the energy gaps:
    E_ij = <|E_i - E_j|>
    It is needed when dephasing_informed option is used
    [ default : NULL ]
  */
  MATRIX* ave_gaps;


  /**
  A flag for NBRA calculations. Since in NBRA, the Hamiltonian is the same for all the trajectories
  we can only compute the Hamiltonian related properties once for one trajectory and increase the speed of calculations.
  If we set the value to 1 it will consider the NBRA type calculations and other integers the Hamiltonian related properties
  are computed for all trajectories.

  Options:
      - 0: no NBRA - Hamiltonians for all trajectories are computed explicitly [ default ]
      - 1: the NBRA is involved - the calculations of the Hamiltonian are conducted for only 1 trajectory, 
           and re-used by all other SH trajectories.  
  */
  int isNBRA;


  ///===============================================================================
  ///================= Entanglement of trajectories ================================
  ///===============================================================================

  /**
    A selector of a method to couple the trajectories in this ensemble.

    Options:
      - 0: no coupling [ default ]
      - 1: ETHD
      - 2: ETHD3 (experimental)
      - 22: another flavor of ETHD3 (experimental)
  */
  int entanglement_opt;


  /**
    Gaussian exponents that dresses up the trajectories in the ETHD3 method
    in the coordinate space, that is   ~exp(-alpha*(R-R0)^2 ) [ default: 0.0 ]
  */
  double ETHD3_alpha;


  /**
    Gaussian exponents that dresses up the trajectories in the ETHD3 method
    in the momentum space, that is   ~exp(-beta*(P-P0)^2 ) [ default: 0.0 ]
  */
  double ETHD3_beta;


  ///===============================================================================
  ///================= QTAG parameters =============================================
  ///===============================================================================

  /**
    How to approximate the Hamiltonian matrix elements for trajectories that belong 
    to different (or same) surfaces

    Options:
      - 0 : BAT [ default ]
      - 1 : LHA
      - 2 : LHAe
      - 3 : BATe 
  */
  int qtag_pot_approx_method;

  ///===============================================================================
  ///================= Bath, Constraints, and Dynamical controls ===================
  ///===============================================================================

  /**
    Temperature of the system. This parameter could be used even in the NVE simulations
    e.g. as a parameters to compute hop acceptance probabilities based on Boltzmann factors [ default: 300 K]
  */ 
  double Temperature;


  /**
    Which ensemble to use in the dynamics. 

    Options:
      - 0: NVE [ default ]
      - 1: NVT
  */
  int ensemble;

  
  /**
    Thermostat parameters 
  */
  bp::dict thermostat_params;


  /**
    Thermostat DOFs

    This list contains the indices of nuclear DOFs which shall be coupled to a thermostat directly.
    [ default: [] ]
  */
  vector<int> thermostat_dofs;


  /**
    Quantum-classical partitioning

    This list of integers contains the indices of nuclear DOFs which chall be treated "quantum-mechanically", well
    including with TSH that is. These DOFs will determine the velocity-rescaling-based acceptance of the hops,
    and these DOFs will be rescaled when the transition is accepted 
    [ default: [0] ]
  */
  vector<int> quantum_dofs;


  /**
    Constrained DOFs

    This list of integers contains the indices of the nuclear DOFs to be constrained - their momenta will be constantly 
    reset to zero, so the corresponding coordinates will stay fixed
    [ default: [] ]
  */
  vector<int> constrained_dofs; 


  /** 
    the nuclear and electronic integration timesteps [ units: a.u. of time, default: 41.0 a.u. = 1 fs ]
  */
  double dt;


  /**
    the number of electronic integration substeps per a nuclear step, such that dt_el = dt_nucl / num_electronic_substeps
  */
  int num_electronic_substeps; 


  /**
    the method for electronic TD-SE integration:

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

  */
  int electronic_integrator;

   
  /**
    If set to True (1), we will force the reprojection matrix T_new to be the identity matrix. This effectively
    removes basis-reprojection (local diabatization) approach and turns on the "naive" approach where
    no trivial crossings exist.
    Default: No (0) - we do want to use the LD approaches by default.
    If Yes (1), one may need to turn on additional state tracking and phase correction methods
  */
  int assume_always_consistent;



  dyn_control_params();
  dyn_control_params(const dyn_control_params& x);
  ~dyn_control_params();

  void sanity_check();
  void set_parameters(bp::dict params);



  friend bool operator == (const dyn_control_params& n1, const dyn_control_params& n2){
    return &n1 == &n2;
  }
  friend bool operator != (const dyn_control_params& n1, const dyn_control_params& n2){
    return !(n1 == n2);  // only compare addresses
  }


};




} // libdyn
}// liblibra

#endif // DYN_CONTROL_PARAMS_H
