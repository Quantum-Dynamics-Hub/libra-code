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

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "../math_linalg/liblinalg.h"


/// liblibra namespace
namespace liblibra{


using namespace liblinalg;

/// libdyn namespace
namespace libdyn{

namespace bp = boost::python;



class dyn_control_params{

  public:



  /**
    The representation to run the Ehrenfest : 0 - diabatic, 1 - adiabatic
  */
  int rep_tdse;


  /**
   The representation of the Hamiltonian update: 0 - diabatic, 1 - adiabatic
   this is the representation in which the computed properties are assumed to be
   for instance, we may have set it to 1, to read the adiabatic energies and couplings,
   to bypass the diabatic-to-adiabatic transformation, which may be useful in some atomistic
   calculations, or with the NBRA
  */
  int rep_ham;


  /**
   The representation to run the SH : 0 - diabatic, 1 - adiabatic
  */
  int rep_sh;


  /**
    The representation to compute LZ probabilitieis: 0 - diabatic, 1- adiabatic
  */
  int rep_lz;


  /**
    Surface hopping mthodology:
    -1 - adiabatic (no need to compute)
     0 - FSSH
     1 - GFSH
     2 - MSSH
     3 - DISH
  */
  int tsh_method;


  /**
    How to compute forces in the dynamics:
      0 - don't compute forces at all - e.g. we do not really need them
      1 - state-specific  as in the TSH or adiabatic (including adiabatic excited states)
      2 - Ehrenfest
  */
  int force_method;

  /**
    How to update NACs and vibronic Hamiltonian before electronic TD-SE propagation
      0 - don't update them (e.g. for simplest NAC)
      1 - update according to changed momentum and existing derivative couplings
  */
  int nac_update_method;


  /**
    In which representation to compute forces:
      0 - diabatic
      1 - adiabatic
  */
  int rep_force;


  /**
    Options to control the acceptance of the proposed hops:
      0 - accept all

      10 - based on adiabatic energy
      11 - based on diabatic energy

      20 - derivative coupling vectors
      21 - difference of state-specific forces

      31 - quantum Boltzmann
      32 - Maxwell-Boltzmann
      33 - updated quantum Boltzmann
  */
  int hop_acceptance_algo;


  /**
    Options to control momenta changes upon successful or frustrated hops:

      0 - don't rescale

      100 - based on adiabatic energy, don't reverse on frustrated hops
      101 - based on adiabatic energy, reverse on frustrated hops
      110 - based on diabatic energy, don't reverse on frustrated hops
      111 - based on diabatic energy, reverse on frustrated hops

      200 - along derivative coupling vectors, don't reverse on frustrated hops
      201 - along derivative coupling vectors, reverse on frustrated hops
      210 - along difference of state-specific forces, don't reverse on frustrated hops
      211 - along difference of state-specific forces, reverse on frustrated hops
  */
  int momenta_rescaling_algo;


  /**
    Whether to scale the SH probabilities by the Boltzmann factor: 0 - do not scale, 1 - scale
  */
  int use_boltz_factor;


  /**
    Temperature of the system
  */
  double Temperature;

/**
   THESE VARIABLES ARE DEPRECATED

//    Do not revert momenta at the frustrated hops, 1 - do revert the momenta
  int do_reverse;


    How to rescale momenta if the hops are successful:

# How to rescale momenta if the hops are successful:
# -1: do not rescale, as in the NBRA [ default ]
#  0: rescale in the diabatic basis - don't care about the
#     velocity directions, just a uniform rescaling
#  1: rescale along the directions of derivative couplings

      0 - rescale along the directions of derivative couplings
      1 - rescale in the diabatic basis - don't care about the velocity directions, just a uniform rescaling,
      2 - do not rescale, as in the NBRA.

  int vel_rescale_opt;
*/

  /**
    integration timestep [units: a.u., default: 41 a.u. = 1 fs]
  */
  double dt;

  /**
    Option to perform the phase correction: 0 - no, 1 - yes (default)
  */
  int do_phase_correction;


  /**
    The minimal magnutude of the matrix element for which we'll be computing the phase correction
    If the overlap is zero, then we don't really care about the phase, but if it is not, then this
    parameter sets out threshold for when we do.  Default: 1e-3
  */
  double phase_correction_tol;


  /**
    State tracking algorithm:
      0 - no state tracking
      1 - method of Kosuke Sato (may fail by getting trapped into an infinite loop)
      2 - Munkres-Kuhn (Hungarian) algorithm (default)
      3 - stochastic reordering
  */
  int state_tracking_algo;

  /**
    Munkres-Kuhn alpha (selects the range of orbitals included in reordering) [default: 0.0]
  */
  double MK_alpha;

  /**
    Munkres-Kuhn verbosity: 0 - no extra output (default), 1 - details
  */
  int MK_verbosity;


  /**
    A selector of a method to couple the trajectories in this ensemble:
      0 - no coupling, 1 - ETHD, 2 - ETHD3 (experimental), 22 - another flavor of ETHD3 (experimental)
  */
  int entanglement_opt;


  /**
    Gaussian exponents that dresses up the trajectories in the ETHD3 method
    in the coordinate space, that is   ~exp(-alpha*(R-R0)^2 )
  */
  double ETHD3_alpha;


  /**
    Gaussian exponents that dresses up the trajectories in the ETHD3 method
    in the momentum space, that is   ~exp(-beta*(P-P0)^2 )
  */
  double ETHD3_beta;


  /**
    Selector of the method to incorporate decoherence:

     -1 - no decoherence [default]
      0 - SDM and alike
      1 - instantaneous decoherence options (ID-S, ID-A, ID-C)
  */
  double decoherence_algo;


  /**
     Corresponds to the "tol" parameter in the sdm function. It controls 
     how much the norm of the old state can be larger than 1.0  before the 
     code stops with the error message

     Default: 0.0
  **/
  double sdm_norm_tolerance;


  /**
    Dephasing rates provided by user
  */
  MATRIX* decoherence_rates;


  /**
    Type of dephasing times/rates calculation:

      0 - use the rates read out from the input  [default]
      1 - use the energy-based decoherence method (EDC)
  */
  int decoherence_times_type;


  /**
    An empirical parameter used in the EDC method: [default = 1.0 Ha]
  */
  double decoherence_C_param;


  /**
    An empirical parameter used in the EDC method: [default = 0.1 Ha]
  */
  double decoherence_eps_param;


  /**
    A flag to apply the dephasing-informed approach of Sifain et al
    to correct dephasing times:

      0 - don't apply [default]
      1 - use it
  */
  int dephasing_informed;


  /**
    A matrix that contains the averaged moduli of the energy gaps:
    E_ij = <|E_i - E_j|>
    It is needed when dephasing_informed option is used
  */
  MATRIX* ave_gaps;


  /**
    Option to control the instantaneous decoherence methodology,
    only used with decoherence_algo == 1

      0 - ID-S
      1 - ID-A [default]
      2 - ID-C - consistent ID - an experimental algorithm
  */
  int instantaneous_decoherence_variant;


  /**
    How to collapse wavefunction amplitudes in the decoherence schemes:
      0 - by rescaling the magnitude of the amplitude vector elements, but preserving "phase" [ default ]
      1 - by resetting the amplitudes to 1.0+0.0j. This option changes phase

  */
  int collapse_option;

  /**
    Ensemble: which ensemble to use: 0 - NVE, 1 - NVT
  */
  int ensemble;


  /**
    Thermostat parameters
  */
  bp::dict thermostat_params;



  dyn_control_params();
  dyn_control_params(const dyn_control_params& x){
    *this = x;
<<<<<<< HEAD
    decoherence_rates = new MATRIX( *x.decoherence_rates );
=======
    decoherence_rates = new MATRIX( *x.decoherence_rates );  
>>>>>>> c3601c27552479ab6c01776ebfbe8bbf69ad8fc8
  }
 ~dyn_control_params() { ;; }

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
