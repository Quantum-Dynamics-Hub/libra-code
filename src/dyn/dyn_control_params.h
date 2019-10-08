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

/// liblibra namespace
namespace liblibra{


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
    Formula for computing SH probabilities: -1 - adiabatic, 0 - FSSH, 1 - GFSH, 2 - MSSH
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
    Whether to scale the SH probabilities by the Boltzmann factor: 0 - do not scale, 1 - scale                            
  */
  int use_boltz_factor;

  
  /**
    Temperature of the system
  */ 
  double Temperature;


  /** 
    Do not revert momenta at the frustrated hops, 1 - do revert the momenta
  */
  int do_reverse;


  /** 
    How to rescale momenta if the hops are successful:
      0 - rescale along the directions of derivative couplings
      1 - rescale in the diabatic basis - don't care about the velocity directions, just a uniform rescaling,
      2 - do not rescale, as in the NBRA.  
  */
  int vel_rescale_opt;


  /** 
    integration timestep [units: a.u., default: 41 a.u. = 1 fs]
  */
  double dt;

  /** 
    Option to perform the phase correction: 0 - no, 1 - yes (default)
  */
  int do_phase_correction;


  /** 
    State tracking algorithm:
      0 - no state tracking
      1 - method of Kosuke Sato (may fail by getting trapped into an infinite loop)
      2 - Munkres-Kuhn (Hungarian) algorithm (default)
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




  dyn_control_params();
  dyn_control_params(const dyn_control_params& x){  *this = x; }
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
