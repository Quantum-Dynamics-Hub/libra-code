/*********************************************************************************
* Copyright (C) 2015-2017 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file libdyn.cpp
  \brief The file implements Python export function
    
*/

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "libdyn.h"


/// liblibra namespace
namespace liblibra{

using namespace boost::python;

/// libdyn namespace
namespace libdyn{

using namespace libnuclear;
using namespace libelectronic;
//using namespace librigidbody;
using namespace libthermostat;
using namespace libbarostat;
using namespace libwfcgrid;
using namespace libwfcgrid2;
using namespace libensemble;
using namespace libgwp;
using namespace libheom;

using namespace libthermostat;




void export_dyn_control_params_objects(){

  // Arbitrary wavefunction
  void (dyn_control_params::*expt_sanity_check_v1)(boost::python::dict params) = &dyn_control_params::set_parameters;
  void (dyn_control_params::*expt_set_parameters_v1)(boost::python::dict params) = &dyn_control_params::set_parameters;


  class_<dyn_control_params>("dyn_control_params",init<>())
      .def("__copy__", &generic__copy__<dyn_control_params>)
      .def("__deepcopy__", &generic__deepcopy__<dyn_control_params>)

      .def_readwrite("rep_tdse", &dyn_control_params::rep_tdse)
      .def_readwrite("rep_ham", &dyn_control_params::rep_ham)
      .def_readwrite("rep_sh", &dyn_control_params::rep_sh)
      .def_readwrite("rep_lz", &dyn_control_params::rep_lz)
      .def_readwrite("tsh_method", &dyn_control_params::tsh_method)
      .def_readwrite("force_method", &dyn_control_params::force_method)
      .def_readwrite("nac_update_method", &dyn_control_params::nac_update_method)
      .def_readwrite("rep_force", &dyn_control_params::rep_force)
      .def_readwrite("hop_acceptance_algo", &dyn_control_params::hop_acceptance_algo)
      .def_readwrite("momenta_rescaling_algo", &dyn_control_params::momenta_rescaling_algo)
      .def_readwrite("use_boltz_factor", &dyn_control_params::use_boltz_factor)
      .def_readwrite("Temperature", &dyn_control_params::Temperature)
//      .def_readwrite("do_reverse", &dyn_control_params::do_reverse)
//      .def_readwrite("vel_rescale_opt", &dyn_control_params::vel_rescale_opt)
      .def_readwrite("dt", &dyn_control_params::dt)
      .def_readwrite("do_phase_correction", &dyn_control_params::do_phase_correction)
      .def_readwrite("phase_correction_tol", &dyn_control_params::phase_correction_tol)
      .def_readwrite("state_tracking_algo", &dyn_control_params::state_tracking_algo)
      .def_readwrite("MK_alpha", &dyn_control_params::MK_alpha)
      .def_readwrite("MK_verbosity", &dyn_control_params::MK_verbosity)
      .def_readwrite("entanglement_opt", &dyn_control_params::entanglement_opt)
      .def_readwrite("ETHD3_alpha", &dyn_control_params::ETHD3_alpha)
      .def_readwrite("ETHD3_beta", &dyn_control_params::ETHD3_beta)
      .def_readwrite("decoherence_algo", &dyn_control_params::decoherence_algo)
      .def_readwrite("decoherence_times_type", &dyn_control_params::decoherence_times_type)
      .def_readwrite("decoherence_C_param", &dyn_control_params::decoherence_C_param)
      .def_readwrite("decoherence_eps_param", &dyn_control_params::decoherence_eps_param)
      .def_readwrite("dephasing_informed", &dyn_control_params::dephasing_informed)
      .def_readwrite("ave_gaps", &dyn_control_params::ave_gaps)
      .def_readwrite("instantaneous_decoherence_variant", &dyn_control_params::instantaneous_decoherence_variant)
      .def_readwrite("collapse_option", &dyn_control_params::collapse_option)
      .def_readwrite("ensemble", &dyn_control_params::ensemble)
      .def_readwrite("thermostat_params", &dyn_control_params::thermostat_params)

      .def("sanity_check", expt_sanity_check_v1)
      .def("set_parameters", expt_set_parameters_v1)
  ;
}


void export_dyn_decoherence_objects(){
 
  //================== ID-A =======================
/*
  int (*expt_ida_v1)(CMATRIX& Coeff, int old_st, int new_st, double E_old, double E_new, double T, double ksi) = &ida;
  def("ida", expt_ida_v1);
*/

  ///=================== dyn_decoherence_methods.cpp =======================

  CMATRIX (*expt_sdm_v1)
  (CMATRIX& Coeff, double dt, int act_st, MATRIX& decoh_rates) = &sdm;
  def("sdm", expt_sdm_v1);

  CMATRIX (*expt_sdm_v2)
  (CMATRIX& Coeff, double dt, vector<int>& act_st, vector<MATRIX>& decoh_rates) = &sdm;
  def("sdm", expt_sdm_v2);

  void (*expt_project_out_v1)(CMATRIX& Coeff, int traj, int i) = &project_out;
  def("project_out", expt_project_out_v1);

  void (*expt_collapse_v1)(CMATRIX& Coeff, int traj, int i, int collapse_option) = &collapse;
  def("collapse", expt_collapse_v1);

  void (*expt_instantaneous_decoherence_v1)(CMATRIX& Coeff, 
  vector<int>& accepted_states, vector<int>& proposed_states, vector<int>& initial_states,
  int instantaneous_decoherence_variant, int collapse_option) = &instantaneous_decoherence;
  def("instantaneous_decoherence", expt_instantaneous_decoherence_v1);


  ///================  In dyn_decoherence_time.cpp  ===================================

  MATRIX (*expt_edc_rates_v1)
  (CMATRIX& Hvib, double Ekin, double C_param, double eps_param) = &edc_rates;
  def("edc_rates", expt_edc_rates_v1);

  vector<MATRIX> (*expt_edc_rates_v2)
  (vector<CMATRIX>& Hvib, vector<double>& Ekin, 
  double C_param, double eps_param) = &edc_rates;
  def("edc_rates", expt_edc_rates_v2);


  void (*expt_dephasing_informed_correction_v1)
  (MATRIX& decoh_rates, CMATRIX& Hvib, MATRIX& ave_gaps) = &dephasing_informed_correction;
  def("dephasing_informed_correction", expt_dephasing_informed_correction_v1);

  void (*expt_dephasing_informed_correction_v2)
  (vector<MATRIX>& decoh_rates, vector<CMATRIX>& Hvib, MATRIX& ave_gaps) = &dephasing_informed_correction;
  def("dephasing_informed_correction", expt_dephasing_informed_correction_v2);

  
  MATRIX (*expt_coherence_intervals_v1)(CMATRIX& Coeff, MATRIX& rates) = &coherence_intervals;
  def("coherence_intervals", expt_coherence_intervals_v1);

  MATRIX (*expt_coherence_intervals_v2)(CMATRIX& Coeff, vector<MATRIX>& rates) = &coherence_intervals;
  def("coherence_intervals", expt_coherence_intervals_v2);




  //================== DISH =======================

/*


  int (*expt_dish_v1)(Electronic& el, MATRIX& t_m, const MATRIX& tau_m, const CMATRIX& Hvib,
          int use_boltz_flag, double Ekin, double T, double ksi1, double ksi2) = &dish;

  int (*expt_dish_v2)(Electronic& el, Nuclear& mol, Hamiltonian& ham, 
          MATRIX& t_m, const MATRIX& tau_m, int use_boltz_flag, double T, double ksi1, double ksi2) = &dish;
  def("dish", expt_dish_v1);
  def("dish", expt_dish_v2);

*/
}


void export_dyn_hop_acceptance_objects(){

  //============= dyn_hop_proposal.cpp ======================

  int (*expt_can_rescale_along_vector_v1)
  (double E_old, double E_new, MATRIX& p, MATRIX& invM, MATRIX& t) = &can_rescale_along_vector;
  def("can_rescale_along_vector", expt_can_rescale_along_vector_v1);

  void (*expt_rescale_along_vector_v1)
  (double E_old, double E_new, MATRIX& p, MATRIX& invM, MATRIX& t, int do_reverse) = &rescale_along_vector;
  def("rescale_along_vector", expt_rescale_along_vector_v1);

  vector<double> (*expt_Boltz_quant_prob_v1)
  (vector<double>& E, double T) = &Boltz_quant_prob;
  def("Boltz_quant_prob", expt_Boltz_quant_prob_v1);

  double (*expt_Boltz_cl_prob_v1)(double E, double T) = &Boltz_cl_prob;
  def("Boltz_cl_prob", expt_Boltz_cl_prob_v1);

  double (*expt_Boltz_cl_prob_up_v1)(double E, double T) = Boltz_cl_prob_up;
  def("Boltz_cl_prob_up", expt_Boltz_cl_prob_up_v1);

  double (*expt_HO_prob_v1)
  (vector<double>& E, vector<int>& qn, double T, vector<double>& prob) = &HO_prob;
  def("HO_prob", expt_HO_prob_v1);

  double (*expt_HO_prob_up_v1)
  (vector<double>& E, vector<int>& qn, double T, vector<double>& prob) = &HO_prob_up;
  def("HO_prob_up", expt_HO_prob_up_v1);

  double (*expt_boltz_factor_v1)
  (double E_new, double E_old, double T, int boltz_opt) = &boltz_factor;
  def("boltz_factor", expt_boltz_factor_v1);

  vector<int> (*expt_accept_hops_v1)
  (dyn_control_params& prms,
   MATRIX& q, MATRIX& p, MATRIX& invM, CMATRIX& C, vector<CMATRIX>& projectors, 
   nHamiltonian& ham, vector<int>& proposed_states, vector<int>& initial_states, Random& rnd ) = &accept_hops;
  def("accept_hops", expt_accept_hops_v1);

}


void export_dyn_hop_proposal_objects(){

  //============= dyn_hop_proposal.cpp ======================

  MATRIX (*expt_hopping_probabilities_fssh_v1)
  (dyn_control_params& prms, CMATRIX& Coeff, CMATRIX& Hvib) = &hopping_probabilities_fssh;
  def("hopping_probabilities_fssh", expt_hopping_probabilities_fssh_v1);

  MATRIX (*expt_hopping_probabilities_gfsh_v1)
  (dyn_control_params& prms, CMATRIX& Coeff, CMATRIX& Hvib) = &hopping_probabilities_gfsh;
  def("hopping_probabilities_gfsh", expt_hopping_probabilities_gfsh_v1);

  MATRIX (*expt_hopping_probabilities_mssh_v1)
  (dyn_control_params& prms, CMATRIX& Coeff, CMATRIX& Hvib) = &hopping_probabilities_mssh;
  def("hopping_probabilities_mssh", expt_hopping_probabilities_mssh_v1);


  vector<MATRIX> (*expt_hop_proposal_probabilities_v1)
  (dyn_control_params& prms,
   MATRIX& q, MATRIX& p, MATRIX& invM, CMATRIX& C, vector<CMATRIX>& projectors,
   nHamiltonian& ham, vector<MATRIX>& prev_ham_dia) = &hop_proposal_probabilities;
  def("hop_proposal_probabilities", expt_hop_proposal_probabilities_v1);

  int (*expt_hop_v1)(int initstate, MATRIX& g, double ksi) = &hop;
  def("hop", expt_hop_v1);

  vector<int> (*expt_propose_hops_v1)
  (vector<MATRIX>& g, vector<int>& act_states, Random& rnd) = &propose_hops;
  def("propose_hops", expt_propose_hops_v1);

}


void export_dyn_methods_objects(){

  vector<int> (*expt_decoherence_event_v1)
  (MATRIX& coherence_time, MATRIX& coherence_interval, Random& rnd) = &decoherence_event;
  def("decoherence_event", expt_decoherence_event_v1);


  vector<int> (*expt_dish_hop_proposal_v1)
  (vector<int>& act_states, CMATRIX& Coeff, 
  MATRIX& coherence_time, vector<MATRIX>& decoherence_rates, Random& rnd) = &dish_hop_proposal;
  def("dish_hop_proposal", expt_dish_hop_proposal_v1);

  void (*expt_dish_project_out_collapse_v1)
  (vector<int>& old_states, vector<int>& proposed_states, vector<int>& new_states, 
  CMATRIX& Coeff, MATRIX& coherence_time, int collapse_option) = &dish_project_out_collapse;
  def("dish_project_out_collapse", expt_dish_project_out_collapse_v1);



}

void export_dyn_projectors_objects(){

  //============= dyn_projectors.cpp ======================

  CMATRIX (*expt_compute_phase_corrections_v1)(CMATRIX& S, double tol) = &compute_phase_corrections;
  def("compute_phase_corrections", expt_compute_phase_corrections_v1);

  CMATRIX (*expt_compute_phase_corrections_v2)(CMATRIX& S) = &compute_phase_corrections;
  def("compute_phase_corrections", expt_compute_phase_corrections_v2);

  vector<int> (*expt_get_reordering_v1)(CMATRIX& time_overlap) = &get_reordering;
  def("get_reordering", expt_get_reordering_v1);  

  MATRIX (*expt_make_cost_mat_v1)(CMATRIX& orb_mat_inp, CMATRIX& en_mat_inp, double alpha) = &make_cost_mat;
  def("make_cost_mat", expt_make_cost_mat_v1);  

  vector<int> (*expt_Munkres_Kuhn_v1)(CMATRIX& orb_mat_inp, CMATRIX& en_mat_inp, double alpha, int verbosity) = &Munkres_Kuhn;
  def("Munkres_Kuhn", expt_Munkres_Kuhn_v1);  

  CMATRIX (*expt_permutation2cmatrix_v1)(vector<int>& permutation) = &permutation2cmatrix;
  def("permutation2cmatrix", expt_permutation2cmatrix_v1);  

  void (*expt_update_projectors_v1)(dyn_control_params& prms, vector<CMATRIX>& projectors, 
  vector<CMATRIX>& Eadi, vector<CMATRIX>& St, Random& rnd) = &update_projectors;
  def("update_projectors", expt_update_projectors_v1);  


  CMATRIX (*expt_raw_to_dynconsyst_v1)
  (CMATRIX& amplitudes, vector<CMATRIX>& projectors) = &raw_to_dynconsyst;
  def("raw_to_dynconsyst", expt_raw_to_dynconsyst_v1);

  CMATRIX (*expt_dynconsyst_to_raw_v1)
  (CMATRIX& amplitudes, vector<CMATRIX>& projectors) = &dynconsyst_to_raw;
  def("dynconsyst_to_raw", expt_dynconsyst_to_raw_v1);




}

void export_LZ_hopping_probabilities_objects(){


  MATRIX (*expt_compute_hopping_probabilities_lz_v1)
  (nHamiltonian& ham, int rep, MATRIX& p, const MATRIX& invM, MATRIX& prev_ham_dia) = &compute_hopping_probabilities_lz;

  def("compute_hopping_probabilities_lz",expt_compute_hopping_probabilities_lz_v1);

}



void export_permutation_objects(){


  vector<int> (*expt_get_permutation_v1)(vector<vector<int> >& inp) = &get_permutation; 
  vector<int> (*expt_Munkres_Kuhn_minimize_v1)(MATRIX& _X, int verbosity) = &Munkres_Kuhn_minimize;
  vector<int> (*expt_Munkres_Kuhn_maximize_v1)(MATRIX& _X, int verbosity) = &Munkres_Kuhn_maximize;

  def("get_permutation", expt_get_permutation_v1);  
  def("Munkres_Kuhn_minimize", expt_Munkres_Kuhn_minimize_v1);  
  def("Munkres_Kuhn_maximize", expt_Munkres_Kuhn_maximize_v1);  


}

void export_Energy_Forces_objects(){


  double (*expt_compute_kinetic_energy_v1)(MATRIX& p, MATRIX& invM) = &compute_kinetic_energy;
  def("compute_kinetic_energy",expt_compute_kinetic_energy_v1);

  vector<double> (*expt_compute_kinetic_energies_v1)(MATRIX& p, MATRIX& invM) = &compute_kinetic_energies;
  def("compute_kinetic_energies",expt_compute_kinetic_energies_v1);

  CMATRIX (*expt_tsh_indx2ampl_v1)(vector<int>& res, int nstates) = &tsh_indx2ampl;
  def("tsh_indx2ampl", expt_tsh_indx2ampl_v1);


  MATRIX (*expt_aux_get_forces_v1)
  (dyn_control_params& prms, CMATRIX& amplitudes, vector<CMATRIX>& projectors, 
  vector<int>& act_states, nHamiltonian& ham) = &aux_get_forces;
  def("aux_get_forces", expt_aux_get_forces_v1);

  MATRIX (*expt_aux_get_forces_v2)
  (bp::dict prms, CMATRIX& amplitudes, vector<CMATRIX>& projectors, 
  vector<int>& act_states, nHamiltonian& ham) = &aux_get_forces;
  def("aux_get_forces", expt_aux_get_forces_v2);


  vector<CMATRIX> (*expt_get_Eadi_v1)(nHamiltonian& ham) = &get_Eadi;
  def("get_Eadi", expt_get_Eadi_v1);





  double (*expt_compute_kinetic_energy_v2)(Nuclear& mol) = &compute_kinetic_energy;
  def("compute_kinetic_energy",expt_compute_kinetic_energy_v2);

  double (*expt_compute_kinetic_energy_v3)(Ensemble& ens) = &compute_kinetic_energy;
  def("compute_kinetic_energy",expt_compute_kinetic_energy_v3);


  double (*expt_compute_potential_energy_v1)(Nuclear& mol, Electronic& el, Hamiltonian& ham, int opt) = &compute_potential_energy;
  def("compute_potential_energy",expt_compute_potential_energy_v1);

  double (*expt_compute_potential_energy_v2)(Ensemble& ens, int opt) = &compute_potential_energy;
  def("compute_potential_energy",expt_compute_potential_energy_v2);

  double (*expt_compute_forces_v1)(Nuclear& mol, Electronic& el, Hamiltonian& ham, int opt) = &compute_forces;
  def("compute_forces",expt_compute_forces_v1);

  double (*expt_compute_forces_v2)(Ensemble& ens, int opt) = &compute_forces;
  def("compute_forces",expt_compute_forces_v2);



}


void export_Dyn_objects(){
/** 
  \brief Exporter of libdyn classes and functions

*/


  export_Nuclear_objects();
  export_Electronic_objects();
  export_Thermostat_objects();
  export_Barostat_objects();
  export_Wfcgrid_objects();
  export_Wfcgrid2_objects();
  export_Ensemble_objects();
  export_gwp_objects();
  export_heom_objects();

  export_dyn_control_params_objects();
  export_dyn_decoherence_objects();
  export_dyn_hop_acceptance_objects();
  export_dyn_hop_proposal_objects();
  export_dyn_methods_objects();
  export_dyn_projectors_objects();

  
  export_LZ_hopping_probabilities_objects();
  export_permutation_objects();

  export_Energy_Forces_objects();



  //============= Dynamics.cpp ======================



  void (*expt_update_Hamiltonian_q_v1)
  (dyn_control_params& prms, MATRIX& q, vector<CMATRIX>& projectors, 
   nHamiltonian& ham, 
   bp::object py_funct, bp::object model_params) = &update_Hamiltonian_q;

  void (*expt_update_Hamiltonian_q_v2)
  (bp::dict prms, MATRIX& q, vector<CMATRIX>& projectors, 
   nHamiltonian& ham, 
   bp::object py_funct, bp::object model_params) = &update_Hamiltonian_q;

  def("update_Hamiltonian_q", expt_update_Hamiltonian_q_v1);
  def("update_Hamiltonian_q", expt_update_Hamiltonian_q_v2);


  void (*expt_update_Hamiltonian_q_ethd_v1)
  (dyn_control_params& prms, MATRIX& q, MATRIX& p, vector<CMATRIX>& projectors,
   nHamiltonian& ham, 
   bp::object py_funct, bp::object model_params, MATRIX& invM) = &update_Hamiltonian_q_ethd;

  void (*expt_update_Hamiltonian_q_ethd_v2)
  (bp::dict prms, MATRIX& q, MATRIX& p, vector<CMATRIX>& projectors,
   nHamiltonian& ham, 
   bp::object py_funct, bp::object model_params, MATRIX& invM) = &update_Hamiltonian_q_ethd;

  def("update_Hamiltonian_q_ethd", expt_update_Hamiltonian_q_ethd_v1);
  def("update_Hamiltonian_q_ethd", expt_update_Hamiltonian_q_ethd_v2);



  void (*expt_update_Hamiltonian_p_v1)
  (dyn_control_params& prms, nHamiltonian& ham, MATRIX& p, MATRIX& invM) = &update_Hamiltonian_p;

  void (*expt_update_Hamiltonian_p_v2)
  (bp::dict prms, nHamiltonian& ham, MATRIX& p, MATRIX& invM) = &update_Hamiltonian_p;

  def("update_Hamiltonian_p", expt_update_Hamiltonian_p_v1);
  def("update_Hamiltonian_p", expt_update_Hamiltonian_p_v2);


  CMATRIX (*expt_transform_amplitudes_v1)
  (int rep_in, int rep_out, CMATRIX& C, nHamiltonian& ham) = &transform_amplitudes;

  def("transform_amplitudes", expt_transform_amplitudes_v1);



  vector<CMATRIX> (*expt_compute_St_v1)(nHamiltonian& ham, vector<CMATRIX>& Uprev) = &compute_St;
  def("compute_St", expt_compute_St_v1);



  void (*expt_compute_dynamics_v1)
  (MATRIX& q, MATRIX& p, MATRIX& invM, CMATRIX& C, vector<CMATRIX>& projectors, vector<int>& act_states,
   nHamiltonian& ham, bp::object py_funct, bp::dict model_params, 
   bp::dict dyn_params, Random& rnd) = &compute_dynamics;
  def("compute_dynamics", expt_compute_dynamics_v1);

  void (*expt_compute_dynamics_v2)
  (MATRIX& q, MATRIX& p, MATRIX& invM, CMATRIX& C, vector<CMATRIX>& projectors, vector<int>& act_states, 
   nHamiltonian& ham, bp::object py_funct, bp::dict model_params, bp::dict dyn_params, Random& rnd, 
   vector<Thermostat>& therm) = &compute_dynamics;
  def("compute_dynamics", expt_compute_dynamics_v2);





}// export_Dyn_objects()




#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygdyn){
#else
BOOST_PYTHON_MODULE(libdyn){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

  export_Dyn_objects();

}


}// libdyn
}// liblibra

