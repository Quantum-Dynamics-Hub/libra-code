/*********************************************************************************
* Copyright (C) 2019-2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file dyn_control_params.cpp
  \brief The file implements the methods to setup control parameters for dynamics
*/

#include "dyn_control_params.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;

/// libdyn namespace
namespace libdyn{

namespace bp = boost::python;

dyn_control_params::dyn_control_params(){

/**

  This function initializes the default values of control parameters

*/


  ///================= Computing Hamiltonian-related properties ====================
  rep_tdse = 1;
//  rep_ham = 0;
  ham_update_method = 1;
  ham_transform_method = 1;
  rep_sh = 1;
  rep_lz = 0;
  rep_force = 1;
  force_method = 1;
  enforce_state_following = 0; 
  enforced_state_index = 0; 
  time_overlap_method = 0;
  nac_update_method = 1;
  nac_algo = -1;
  hvib_update_method = 1;
  do_ssy = 0;
  do_phase_correction = 1;
  phase_correction_tol = 1e-3; 
  state_tracking_algo = 2;
  MK_alpha = 0.0;
  MK_verbosity = 0;
  convergence = 0;
  max_number_attempts = 100;
  min_probability_reordering = 0.0;
  isNBRA = 0;

  ///================= Surface hopping: proposal, acceptance =======================
  tsh_method = 0;
  hop_acceptance_algo = 0;
  momenta_rescaling_algo = 0;
  use_boltz_factor = 0;

  ///================= Decoherence options =========================================
  decoherence_algo = -1; 
  sdm_norm_tolerance = 0.0;
  dish_decoherence_event_option = 1;
  decoherence_times_type = -1;
  schwartz_decoherence_inv_alpha = NULL;
  decoherence_C_param = 1.0;
  decoherence_eps_param = 0.1;
  dephasing_informed = 0;
  instantaneous_decoherence_variant = 1; 
  collapse_option = 0;
  decoherence_rates = NULL;
  ave_gaps = NULL;

  ///================= Entanglement of trajectories ================================
  entanglement_opt = 0;
  ETHD3_alpha = 1.0;
  ETHD3_beta = 1.0;

  ///============================ QTAG =============================================
  qtag_pot_approx_method = 0;
  ///================= Bath, Constraints, and Dynamical controls ===================

  Temperature = 300.0;
  ensemble = 0;
  thermostat_params = bp::dict();

  thermostat_dofs = vector<int>();
  quantum_dofs = vector<int>(1, 0);
  constrained_dofs = vector<int>();

  dt = 41.0;
  num_electronic_substeps = 1;
  electronic_integrator = 0;
  assume_always_consistent = 0;
}


dyn_control_params::dyn_control_params(const dyn_control_params& x){ 
  //cout<<"dyn_control_params cctor\n";

  rep_tdse = x.rep_tdse;
//  rep_ham = x.rep_ham;
  ham_update_method = x.ham_update_method;
  ham_transform_method = x.ham_transform_method;
  rep_sh = x.rep_sh;
  rep_lz = x.rep_lz;
  rep_force = x.rep_force;
  force_method = x.force_method;
  enforce_state_following = x.enforce_state_following;
  enforced_state_index = x.enforced_state_index; 
  time_overlap_method = x.time_overlap_method;
  nac_update_method = x.nac_update_method;
  nac_algo = x.nac_algo;
  hvib_update_method = x.hvib_update_method;
  do_ssy = x.do_ssy;
  do_phase_correction = x.do_phase_correction;
  phase_correction_tol = x.phase_correction_tol; 
  state_tracking_algo = x.state_tracking_algo;
  MK_alpha = x.MK_alpha;
  MK_verbosity = x.MK_verbosity;
  convergence = x.convergence;
  max_number_attempts = x.max_number_attempts;
  min_probability_reordering = x.min_probability_reordering;
  isNBRA = x.isNBRA;

  ///================= Surface hopping: proposal, acceptance =======================
  tsh_method = x.tsh_method;
  hop_acceptance_algo = x.hop_acceptance_algo;
  momenta_rescaling_algo = x.momenta_rescaling_algo;
  use_boltz_factor = x.use_boltz_factor;

  ///================= Decoherence options =========================================
  decoherence_algo = x.decoherence_algo; 
  sdm_norm_tolerance = x.sdm_norm_tolerance;
  dish_decoherence_event_option = x.dish_decoherence_event_option;
  decoherence_times_type = x.decoherence_times_type;
  decoherence_C_param = x.decoherence_C_param;
  decoherence_eps_param = x.decoherence_eps_param;
  dephasing_informed = x.dephasing_informed;
  instantaneous_decoherence_variant = x.instantaneous_decoherence_variant; 
  collapse_option = x.collapse_option;

  ///================= Entanglement of trajectories ================================
  entanglement_opt = x.entanglement_opt;
  ETHD3_alpha = x.ETHD3_alpha;
  ETHD3_beta = x.ETHD3_beta;

  ///============================ QTAG =============================================
  qtag_pot_approx_method = x.qtag_pot_approx_method;
  ///================= Bath, Constraints, and Dynamical controls ===================

  Temperature = x.Temperature;
  ensemble = x.ensemble;

  thermostat_params = bp::dict(x.thermostat_params);

  thermostat_dofs = x.thermostat_dofs;
  quantum_dofs = x.quantum_dofs;
  constrained_dofs = x.constrained_dofs;

  dt = x.dt;
  num_electronic_substeps = x.num_electronic_substeps;
  electronic_integrator = x.electronic_integrator;
  assume_always_consistent = x. assume_always_consistent;

  decoherence_rates = new MATRIX(x.decoherence_rates->n_rows, x.decoherence_rates->n_cols);  
  *decoherence_rates = *x.decoherence_rates;

  ave_gaps = new MATRIX( x.ave_gaps->n_rows, x.ave_gaps->n_cols );
  *ave_gaps = *x.ave_gaps;

  schwartz_decoherence_inv_alpha = new MATRIX( x.schwartz_decoherence_inv_alpha->n_rows, x.schwartz_decoherence_inv_alpha->n_cols );
  *schwartz_decoherence_inv_alpha = *x.schwartz_decoherence_inv_alpha;


}
dyn_control_params::~dyn_control_params() {  

  //cout<<"dyn_control_params destructor\n";

  delete decoherence_rates;  
  delete ave_gaps;
  delete schwartz_decoherence_inv_alpha;
}


void dyn_control_params::sanity_check(){

  ///=================== Options for state tracking ======================
  if(state_tracking_algo==0 || state_tracking_algo==1 ||
     state_tracking_algo==2 || state_tracking_algo==3 ||
     state_tracking_algo==32 || state_tracking_algo==33){ ; ; }
  else{
    std::cout<<"Error in dyn_control_params::sanity_check: state_tracking_algo = "
        <<state_tracking_algo<<" is not allowed\nExiting...\n";
    exit(0);
  }

  /// Shall not use DISH with the decoherence correction
  if(tsh_method==3 && decoherence_algo!=-1){
    cout<<"Error in dyn_control_params::sanity_check: Shall not use DISH (tsh_method == 3) with the \
           decoherence correction (any but decoherence_algo == -1)\nExiting...\n";
    exit(0);
  }

  /// Check that we don't constrian the thermostatted DOFs and also quantum DOFs
  int sz1 = constrained_dofs.size();
  if(sz1!=0){

    int sz2 = thermostat_dofs.size();

    for(int i1=0; i1<sz1; i1++){
      for(int i2=0; i2<sz2; i2++){
        if( constrained_dofs[i1] == thermostat_dofs[i2] ){
           cout<<"Error in dyn_control_params::sanity_check: can not constrain the thermostatted DOF "<<thermostat_dofs[i2]<<endl;
           cout<<"Exiting...\n";
           exit(0);
        }
      }// for i2
    }// for i1

    sz2 = quantum_dofs.size();
    for(int i1=0; i1<sz1; i1++){
      for(int i2=0; i2<sz2; i2++){
        if( constrained_dofs[i1] == quantum_dofs[i2] ){
           cout<<"Error in dyn_control_params::sanity_check: can not constrain the quantum DOF "<<quantum_dofs[i2]<<endl;
           cout<<"Exiting...\n";
           exit(0);
        }
      }// for i2
    }// for i1

  }// sz1!=0


  if(num_electronic_substeps<=0){
      cout<<"Error in dyn_control_params::sanity_check: num_electronic_substeps = "<<num_electronic_substeps
          <<" should be a positive integer"<<endl;
      cout<<"Exiting...\n";
  }

}



void dyn_control_params::set_parameters(bp::dict params){
/**
  Extract the parameters from the input dictionary
*/

  std::string key;
  for(int i=0;i<len(params.values());i++){
    key = bp::extract<std::string>(params.keys()[i]);

    ///================= Computing Hamiltonian-related properties ====================
    if(key=="rep_tdse") { rep_tdse = bp::extract<int>(params.values()[i]); }
//    else if(key=="rep_ham") { rep_ham = bp::extract<int>(params.values()[i]);   }
    else if(key=="ham_update_method") { ham_update_method = bp::extract<int>(params.values()[i]);   }
    else if(key=="ham_transform_method") { ham_transform_method = bp::extract<int>(params.values()[i]);   }
    else if(key=="rep_sh") { rep_sh = bp::extract<int>(params.values()[i]);  }
    else if(key=="rep_lz") { rep_lz = bp::extract<int>(params.values()[i]);  }
    else if(key=="rep_force") { rep_force = bp::extract<int>(params.values()[i]);  }
    else if(key=="force_method") { force_method = bp::extract<int>(params.values()[i]);  }
    else if(key=="enforce_state_following") { enforce_state_following = bp::extract<int>(params.values()[i]);  }
    else if(key=="enforced_state_index") { enforced_state_index = bp::extract<int>(params.values()[i]);  }
    else if(key=="time_overlap_method"){ time_overlap_method = bp::extract<double>(params.values()[i]); }
    else if(key=="nac_update_method") { nac_update_method = bp::extract<int>(params.values()[i]);  }
    else if(key=="nac_algo") { nac_algo = bp::extract<int>(params.values()[i]);  }
    else if(key=="hvib_update_method") { hvib_update_method = bp::extract<int>(params.values()[i]);   }
    else if(key=="do_ssy") { do_ssy = bp::extract<int>(params.values()[i]);   }
    else if(key=="do_phase_correction") { do_phase_correction = bp::extract<int>(params.values()[i]);  }
    else if(key=="phase_correction_tol") { phase_correction_tol = bp::extract<double>(params.values()[i]);  }
    else if(key=="state_tracking_algo"){  state_tracking_algo = bp::extract<int>(params.values()[i]);  }
    else if(key=="MK_alpha") { MK_alpha = bp::extract<double>(params.values()[i]);  }
    else if(key=="MK_verbosity") { MK_verbosity = bp::extract<int>(params.values()[i]);  }
    else if(key=="convergence") { convergence = bp::extract<int>(params.values()[i]);  }
    else if(key=="max_number_attempts") { max_number_attempts = bp::extract<int>(params.values()[i]);  }
    else if(key=="min_probability_reordering") { min_probability_reordering = bp::extract<double>(params.values()[i]);  }
    else if(key=="isNBRA") { isNBRA = bp::extract<int>(params.values()[i]);   }

    ///================= Surface hopping: proposal, acceptance =======================
    else if(key=="tsh_method") { tsh_method = bp::extract<int>(params.values()[i]);  }
    else if(key=="hop_acceptance_algo") { hop_acceptance_algo = bp::extract<int>(params.values()[i]);  }
    else if(key=="momenta_rescaling_algo"){ momenta_rescaling_algo = bp::extract<int>(params.values()[i]);  }
    else if(key=="use_boltz_factor"){ use_boltz_factor = bp::extract<int>(params.values()[i]);  }

    ///================= Decoherence options =========================================
    else if(key=="decoherence_algo"){ decoherence_algo = bp::extract<int>(params.values()[i]); }
    else if(key=="sdm_norm_tolerance"){ sdm_norm_tolerance = bp::extract<double>(params.values()[i]); }
    else if(key=="dish_decoherence_event_option"){ dish_decoherence_event_option = bp::extract<int>(params.values()[i]); }
    else if(key=="decoherence_times_type"){ decoherence_times_type = bp::extract<int>(params.values()[i]); }
    else if(key=="schwartz_decoherence_inv_alpha"){ 
      MATRIX x( bp::extract<MATRIX>(params.values()[i]) );
      schwartz_decoherence_inv_alpha = new MATRIX(x.n_rows, x.n_cols);      
      for(int a=0;a<x.n_rows;a++){
        for(int b=0;b<x.n_cols;b++){ schwartz_decoherence_inv_alpha->set(a, b, x.get(a,b));   }
      } 
    }
    else if(key=="decoherence_C_param"){ decoherence_C_param = bp::extract<double>(params.values()[i]); }
    else if(key=="decoherence_eps_param"){ decoherence_eps_param = bp::extract<double>(params.values()[i]); }
    else if(key=="dephasing_informed"){ dephasing_informed = bp::extract<int>(params.values()[i]); }
    else if(key=="instantaneous_decoherence_variant"){ instantaneous_decoherence_variant = bp::extract<int>(params.values()[i]); }
    else if(key=="collapse_option"){ collapse_option = bp::extract<int>(params.values()[i]); }
    else if(key=="decoherence_rates"){ 
      MATRIX x( bp::extract<MATRIX>(params.values()[i]) );
      decoherence_rates = new MATRIX(x.n_rows, x.n_cols);      
      for(int a=0;a<x.n_rows;a++){
        for(int b=0;b<x.n_cols;b++){ decoherence_rates->set(a, b, x.get(a,b));   }
      } 
    }
    else if(key=="ave_gaps"){ 
      MATRIX x( bp::extract<MATRIX>(params.values()[i]) );
      ave_gaps = new MATRIX(x.n_rows, x.n_cols);      
      for(int a=0;a<x.n_rows;a++){
        for(int b=0;b<x.n_cols;b++){ ave_gaps->set(a, b, x.get(a,b));   }
      } 
    }

    ///================= Entanglement of trajectories ================================
    else if(key=="entanglement_opt"){ entanglement_opt = bp::extract<int>(params.values()[i]); }
    else if(key=="ETHD3_alpha") { ETHD3_alpha = bp::extract<double>(params.values()[i]);   }
    else if(key=="ETHD3_beta") { ETHD3_beta = bp::extract<double>(params.values()[i]);   }

    ///================= Entanglement of trajectories ================================
    else if(key=="qtag_pot_approx_method"){ qtag_pot_approx_method = bp::extract<int>(params.values()[i]); }
    
    ///================= Bath, Constraints, and Dynamical controls ===================
    else if(key=="Temperature") { Temperature = bp::extract<double>(params.values()[i]);  }
    else if(key=="ensemble"){ ensemble = bp::extract<int>(params.values()[i]); }    
    else if(key=="thermostat_params"){ thermostat_params = bp::extract<bp::dict>(params.values()[i]); }    
    else if(key=="thermostat_dofs"){  
      thermostat_dofs.clear();
      boost::python::list tmp = extract<boost::python::list>(params.values()[i]);
      for(int j=0; j<len(tmp); j++){  thermostat_dofs.push_back( extract<double>(tmp[j]) );  }
    }
    else if(key=="quantum_dofs"){  
      quantum_dofs.clear();
      boost::python::list tmp = extract<boost::python::list>(params.values()[i]);
      for(int j=0; j<len(tmp); j++){  quantum_dofs.push_back( extract<double>(tmp[j]) );  }
    }
    else if(key=="constrained_dofs"){  
      constrained_dofs.clear();
      boost::python::list tmp = extract<boost::python::list>(params.values()[i]);
      for(int j=0; j<len(tmp); j++){  constrained_dofs.push_back( extract<double>(tmp[j]) );  }
    }
    else if(key=="dt") { dt = bp::extract<double>(params.values()[i]);  }
    else if(key=="num_electronic_substeps") { num_electronic_substeps = bp::extract<int>(params.values()[i]);  }
    else if(key=="electronic_integrator"){ electronic_integrator = bp::extract<int>(params.values()[i]); }
    else if(key=="assume_always_consistent"){  assume_always_consistent = bp::extract<int>(params.values()[i]); }

  }// for i

  sanity_check();

}


}// namespace libdyn
}// liblibra

