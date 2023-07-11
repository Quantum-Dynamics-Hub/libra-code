/*********************************************************************************
* Copyright (C) 2015-2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file libdyn.cpp
  \brief The file implements Python export function
    
*/

#if defined(USING_PCH)
#include "../pch.h"
#else
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#endif 

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
//using namespace libensemble;
using namespace libgwp;
using namespace libheom;
using namespace libqtag;

using namespace libthermostat;




void export_dyn_control_params_objects(){

  // Arbitrary wavefunction
  void (dyn_control_params::*expt_sanity_check_v1)(boost::python::dict params) = &dyn_control_params::set_parameters;
  void (dyn_control_params::*expt_set_parameters_v1)(boost::python::dict params) = &dyn_control_params::set_parameters;


  class_<dyn_control_params>("dyn_control_params",init<>())
      .def("__copy__", &generic__copy__<dyn_control_params>)
      .def("__deepcopy__", &generic__deepcopy__<dyn_control_params>)

      ///================= Computing Hamiltonian-related properties ====================
      .def_readwrite("rep_tdse", &dyn_control_params::rep_tdse)
//      .def_readwrite("rep_ham", &dyn_control_params::rep_ham)
      .def_readwrite("ham_update_method", &dyn_control_params::ham_update_method)    
      .def_readwrite("ham_transform_method", &dyn_control_params::ham_transform_method)    
      .def_readwrite("rep_sh", &dyn_control_params::rep_sh)
      .def_readwrite("rep_lz", &dyn_control_params::rep_lz)
      .def_readwrite("rep_force", &dyn_control_params::rep_force)
      .def_readwrite("force_method", &dyn_control_params::force_method)
      .def_readwrite("time_overlap_method", &dyn_control_params::time_overlap_method)
      .def_readwrite("nac_update_method", &dyn_control_params::nac_update_method)
      .def_readwrite("nac_algo", &dyn_control_params::nac_algo)
      .def_readwrite("hvib_update_method", &dyn_control_params::hvib_update_method)    
      .def_readwrite("do_phase_correction", &dyn_control_params::do_phase_correction)
      .def_readwrite("phase_correction_tol", &dyn_control_params::phase_correction_tol)
      .def_readwrite("state_tracking_algo", &dyn_control_params::state_tracking_algo)
      .def_readwrite("MK_alpha", &dyn_control_params::MK_alpha)
      .def_readwrite("MK_verbosity", &dyn_control_params::MK_verbosity)
      .def_readwrite("convergence", &dyn_control_params::convergence)
      .def_readwrite("max_number_attempts", &dyn_control_params::max_number_attempts)
      .def_readwrite("isNBRA", &dyn_control_params::isNBRA)

      ///================= Surface hopping: proposal, acceptance =======================
      .def_readwrite("tsh_method", &dyn_control_params::tsh_method)
      .def_readwrite("hop_acceptance_algo", &dyn_control_params::hop_acceptance_algo)
      .def_readwrite("momenta_rescaling_algo", &dyn_control_params::momenta_rescaling_algo)
      .def_readwrite("use_boltz_factor", &dyn_control_params::use_boltz_factor)

      ///================= Decoherence options =========================================
      .def_readwrite("decoherence_algo", &dyn_control_params::decoherence_algo)
      .def_readwrite("sdm_norm_tolerance", &dyn_control_params::sdm_norm_tolerance)
      .def_readwrite("dish_decoherence_event_option", &dyn_control_params::dish_decoherence_event_option)
      .def_readwrite("decoherence_times_type", &dyn_control_params::decoherence_times_type)
      .def_readwrite("schwartz_decoherence_inv_alpha", &dyn_control_params::schwartz_decoherence_inv_alpha)
      .def_readwrite("decoherence_C_param", &dyn_control_params::decoherence_C_param)
      .def_readwrite("decoherence_eps_param", &dyn_control_params::decoherence_eps_param)
      .def_readwrite("dephasing_informed", &dyn_control_params::dephasing_informed)
      .def_readwrite("instantaneous_decoherence_variant", &dyn_control_params::instantaneous_decoherence_variant)
      .def_readwrite("collapse_option", &dyn_control_params::collapse_option)
      .def_readwrite("decoherence_rates", &dyn_control_params::decoherence_rates)
      .def_readwrite("ave_gaps", &dyn_control_params::ave_gaps)

      ///================= Entanglement of trajectories ================================
      .def_readwrite("entanglement_opt", &dyn_control_params::entanglement_opt)
      .def_readwrite("ETHD3_alpha", &dyn_control_params::ETHD3_alpha)
      .def_readwrite("ETHD3_beta", &dyn_control_params::ETHD3_beta)


      ///================= Bath, Constraints, and Dynamical controls ===================
      .def_readwrite("Temperature", &dyn_control_params::Temperature)
      .def_readwrite("ensemble", &dyn_control_params::ensemble)
      .def_readwrite("thermostat_params", &dyn_control_params::thermostat_params)
      .def_readwrite("thermostat_dofs", &dyn_control_params::thermostat_dofs)
      .def_readwrite("quantum_dofs", &dyn_control_params::quantum_dofs)
      .def_readwrite("constrained_dofs", &dyn_control_params::constrained_dofs)
      .def_readwrite("dt", &dyn_control_params::dt)
      .def_readwrite("num_electronic_substeps", &dyn_control_params::num_electronic_substeps)
      .def_readwrite("electronic_integrator", &dyn_control_params::electronic_integrator)
      .def_readwrite("assume_always_consistent", &dyn_control_params::assume_always_consistent)

      .def("sanity_check", expt_sanity_check_v1)
      .def("set_parameters", expt_set_parameters_v1)
  ;
}




void export_dyn_variables_objects(){


  CMATRIX (*expt_transform_amplitudes_v1)
  (int rep_in, int rep_out, CMATRIX& C, nHamiltonian& ham) = &transform_amplitudes;

  def("transform_amplitudes", expt_transform_amplitudes_v1);

  CMATRIX (*expt_orthogonalized_T_v1)(CMATRIX& T) = &orthogonalized_T;
  def("orthogonalized_T", expt_orthogonalized_T_v1);


  CMATRIX (dyn_variables::*expt_get_dm_adi_v1)(int i, int prev_steps) = &dyn_variables::get_dm_adi;
  CMATRIX (dyn_variables::*expt_get_dm_adi_v2)(int i) = &dyn_variables::get_dm_adi;

  CMATRIX (dyn_variables::*expt_get_dm_dia_v1)(int i, int prev_steps) = &dyn_variables::get_dm_dia;
  CMATRIX (dyn_variables::*expt_get_dm_dia_v2)(int i) = &dyn_variables::get_dm_dia;


  // Arbitrary wavefunction
  void (dyn_variables::*expt_set_parameters_v1)(boost::python::dict params) = &dyn_variables::set_parameters;

  void (dyn_variables::*expt_update_amplitudes_v1)
  (dyn_control_params& dyn_params, nHamiltonian& ham) = &dyn_variables::update_amplitudes;
  void (dyn_variables::*expt_update_amplitudes_v2)
  (bp::dict dyn_params, nHamiltonian& ham) = &dyn_variables::update_amplitudes;
  void (dyn_variables::*expt_update_amplitudes_v3)
  (dyn_control_params& dyn_params, bp::object compute_model, bp::dict model_params) = &dyn_variables::update_amplitudes;
  void (dyn_variables::*expt_update_amplitudes_v4)
  (bp::dict dyn_params, bp::object compute_model, bp::dict model_params) = &dyn_variables::update_amplitudes;


  void (dyn_variables::*expt_update_density_matrix_v1)
  (dyn_control_params& dyn_params, nHamiltonian& ham, int lvl) = &dyn_variables::update_density_matrix;
  void (dyn_variables::*expt_update_density_matrix_v2)
  (bp::dict dyn_params, nHamiltonian& ham, int lvl) = &dyn_variables::update_density_matrix;
  void (dyn_variables::*expt_update_density_matrix_v3)
  (dyn_control_params& dyn_params, bp::object compute_model, bp::dict model_params, int lvl) = &dyn_variables::update_density_matrix;
  void (dyn_variables::*expt_update_density_matrix_v4)
  (bp::dict dyn_params, bp::object compute_model, bp::dict model_params, int lvl) = &dyn_variables::update_density_matrix;


  double (dyn_variables::*expt_compute_average_kinetic_energy_v1)() = &dyn_variables::compute_average_kinetic_energy;
  double (dyn_variables::*expt_compute_average_kinetic_energy_v2)(vector<int>& which_dofs) = &dyn_variables::compute_average_kinetic_energy;
  double (dyn_variables::*expt_compute_kinetic_energy_v1)(int itraj) = &dyn_variables::compute_kinetic_energy;
  double (dyn_variables::*expt_compute_kinetic_energy_v2)(int itraj, vector<int>& which_dofs) = &dyn_variables::compute_kinetic_energy;
  vector<double> (dyn_variables::*expt_compute_kinetic_energies_v1)() = &dyn_variables::compute_kinetic_energies;
  vector<double> (dyn_variables::*expt_compute_kinetic_energies_v2)(vector<int>& which_dofs) = &dyn_variables::compute_kinetic_energies;



  class_<dyn_variables>("dyn_variables",init<int, int, int, int>())
      .def("__copy__", &generic__copy__<dyn_variables>)
      .def("__deepcopy__", &generic__deepcopy__<dyn_variables>)

      ///================= Dimension numbers ===================
      .def_readwrite("ndia", &dyn_variables::ndia)
      .def_readwrite("nadi", &dyn_variables::nadi)
      .def_readwrite("ndof", &dyn_variables::ndof)
      .def_readwrite("ntraj", &dyn_variables::ntraj)

      .def_readwrite("electronic_vars_status", &dyn_variables::electronic_vars_status)
      .def_readwrite("act_states", &dyn_variables::act_states)
      .def_readwrite("nuclear_vars_status", &dyn_variables::nuclear_vars_status)
      .def_readwrite("afssh_vars_status", &dyn_variables::afssh_vars_status)
      .def_readwrite("bcsh_vars_status", &dyn_variables::bcsh_vars_status)
      .def_readwrite("dish_vars_status", &dyn_variables::dish_vars_status)
      .def_readwrite("fssh2_vars_status", &dyn_variables::fssh2_vars_status)


      .def("set_parameters", expt_set_parameters_v1)

      .def("allocate_electronic_vars", &dyn_variables::allocate_electronic_vars)
      .def("allocate_nuclear_vars", &dyn_variables::allocate_nuclear_vars)
      .def("allocate_afssh", &dyn_variables::allocate_afssh)
      .def("allocate_bcsh", &dyn_variables::allocate_bcsh)
      .def("allocate_dish", &dyn_variables::allocate_dish)

      .def("set_q", &dyn_variables::set_q)
      .def("set_p", &dyn_variables::set_p)
      .def("set_f", &dyn_variables::set_f)
      .def("get_ampl_adi", &dyn_variables::get_ampl_adi)
      .def("get_ampl_dia", &dyn_variables::get_ampl_dia)
      .def("get_dm_adi", expt_get_dm_adi_v1)
      .def("get_dm_adi", expt_get_dm_adi_v2)
      .def("get_dm_dia", expt_get_dm_dia_v1)
      .def("get_dm_dia", expt_get_dm_dia_v2)
      .def("get_imass", &dyn_variables::get_imass)
      .def("get_coords", &dyn_variables::get_coords)
      .def("get_momenta", &dyn_variables::get_momenta)
      .def("get_forces", &dyn_variables::get_forces)

      .def("init_nuclear_dyn_var", &dyn_variables::init_nuclear_dyn_var)
      .def("compute_average_kinetic_energy", expt_compute_average_kinetic_energy_v1)
      .def("compute_average_kinetic_energy", expt_compute_average_kinetic_energy_v2)
      .def("compute_kinetic_energy", expt_compute_kinetic_energy_v1)
      .def("compute_kinetic_energy", expt_compute_kinetic_energy_v2)
      .def("compute_kinetic_energies", expt_compute_kinetic_energies_v1)
      .def("compute_kinetic_energies", expt_compute_kinetic_energies_v2)


      .def("update_amplitudes", expt_update_amplitudes_v1)
      .def("update_amplitudes", expt_update_amplitudes_v2)
      .def("update_amplitudes", expt_update_amplitudes_v3)
      .def("update_amplitudes", expt_update_amplitudes_v4)

      .def("update_density_matrix", expt_update_density_matrix_v1)
      .def("update_density_matrix", expt_update_density_matrix_v2)
      .def("update_density_matrix", expt_update_density_matrix_v3)
      .def("update_density_matrix", expt_update_density_matrix_v4)

      .def("update_active_states", &dyn_variables::update_active_states)

      .def("init_amplitudes", &dyn_variables::init_amplitudes)
      .def("init_density_matrix", &dyn_variables::init_density_matrix)
      .def("init_active_states", &dyn_variables::init_active_states)
      .def("init_electronic_dyn_var", &dyn_variables::init_electronic_dyn_var)

      .def("compute_average_dm", &dyn_variables::compute_average_dm)
      .def("compute_average_se_pop", &dyn_variables::compute_average_se_pop)
      .def("compute_average_sh_pop", &dyn_variables::compute_average_sh_pop)

      .def("save_curr_dm_into_prev", &dyn_variables::save_curr_dm_into_prev)

  ;
}




void export_dyn_decoherence_objects(){
 
  //================== ID-A =======================

  ///=================== dyn_decoherence_methods.cpp =======================

  CMATRIX (*expt_sdm_v1)
  (CMATRIX& Coeff, double dt, int act_st, MATRIX& decoh_rates, double tol) = &sdm;
  def("sdm", expt_sdm_v1);

  CMATRIX (*expt_sdm_v2)
  (CMATRIX& Coeff, double dt, int act_st, MATRIX& decoh_rates) = &sdm;
  def("sdm", expt_sdm_v2);


  CMATRIX (*expt_sdm_v3)
  (CMATRIX& Coeff, double dt, vector<int>& act_st, vector<MATRIX>& decoh_rates, double tol, int isNBRA) = &sdm;
  def("sdm", expt_sdm_v3);

  CMATRIX (*expt_sdm_v4)
  (CMATRIX& Coeff, double dt, vector<int>& act_st, vector<MATRIX>& decoh_rates, double tol) = &sdm;
  def("sdm", expt_sdm_v4);

  CMATRIX (*expt_sdm_v5)
  (CMATRIX& Coeff, double dt, vector<int>& act_st, vector<MATRIX>& decoh_rates) = &sdm;
  def("sdm", expt_sdm_v5);




  void (*expt_project_out_v1)(CMATRIX& Coeff, int traj, int i) = &project_out;
  def("project_out", expt_project_out_v1);

  void (*expt_collapse_v1)(CMATRIX& Coeff, int traj, int i, int collapse_option) = &collapse;
  def("collapse", expt_collapse_v1);

  void (*expt_instantaneous_decoherence_v1)(CMATRIX& Coeff, 
  vector<int>& accepted_states, vector<int>& proposed_states, vector<int>& initial_states,
  int instantaneous_decoherence_variant, int collapse_option) = &instantaneous_decoherence;
  def("instantaneous_decoherence", expt_instantaneous_decoherence_v1);



  void (*expt_wp_reversal_events_v1)
  (dyn_variables& dyn_var, nHamiltonian& ham, double dt) = &wp_reversal_events;
  def("wp_reversal_events", expt_wp_reversal_events_v1);

  CMATRIX (*expt_bcsh_v1)
  (CMATRIX& Coeff, double dt, vector<int>& act_states, MATRIX& reversal_events) = &bcsh;
  def("bcsh", expt_bcsh_v1);


  CMATRIX (*expt_mfsd_v1)
  (MATRIX& p, CMATRIX& Coeff, MATRIX& invM, double dt, vector<MATRIX>& decoherence_rates, 
   nHamiltonian& ham, Random& rnd, int isNBRA) = &mfsd;
  def("mfsd", expt_mfsd_v1);

  CMATRIX (*expt_mfsd_v2)
  (MATRIX& p, CMATRIX& Coeff, MATRIX& invM, double dt, vector<MATRIX>& decoherence_rates, 
   nHamiltonian& ham, Random& rnd) = &mfsd;
  def("mfsd", expt_mfsd_v2);


  ///================  In dyn_decoherence_time.cpp  ===================================

  MATRIX (*expt_edc_rates_v1)
  (CMATRIX& Hvib, double Ekin, double C_param, double eps_param, int isNBRA) = &edc_rates;
  def("edc_rates", expt_edc_rates_v1);

  MATRIX (*expt_edc_rates_v2)
  (CMATRIX& Hvib, double Ekin, double C_param, double eps_param) = &edc_rates;
  def("edc_rates", expt_edc_rates_v2);


  vector<MATRIX> (*expt_edc_rates_v3)
  (vector<CMATRIX>& Hvib, vector<double>& Ekin, 
  double C_param, double eps_param, int isNBRA) = &edc_rates;
  def("edc_rates", expt_edc_rates_v3);

  vector<MATRIX> (*expt_edc_rates_v4)
  (vector<CMATRIX>& Hvib, vector<double>& Ekin, 
  double C_param, double eps_param, int isNBRA) = &edc_rates;
  def("edc_rates", expt_edc_rates_v4);



  void (*expt_dephasing_informed_correction_v1)
  (MATRIX& decoh_rates, CMATRIX& Hvib, MATRIX& ave_gaps, int isNBRA) = &dephasing_informed_correction;
  def("dephasing_informed_correction", expt_dephasing_informed_correction_v1);

  void (*expt_dephasing_informed_correction_v2)
  (MATRIX& decoh_rates, CMATRIX& Hvib, MATRIX& ave_gaps) = &dephasing_informed_correction;
  def("dephasing_informed_correction", expt_dephasing_informed_correction_v2);


  void (*expt_dephasing_informed_correction_v3)
  (vector<MATRIX>& decoh_rates, vector<CMATRIX>& Hvib, MATRIX& ave_gaps, int isNBRA) = &dephasing_informed_correction;
  def("dephasing_informed_correction", expt_dephasing_informed_correction_v3);

  void (*expt_dephasing_informed_correction_v4)
  (vector<MATRIX>& decoh_rates, vector<CMATRIX>& Hvib, MATRIX& ave_gaps) = &dephasing_informed_correction;
  def("dephasing_informed_correction", expt_dephasing_informed_correction_v4);

  

  MATRIX (*expt_coherence_intervals_v1)(CMATRIX& Coeff, MATRIX& rates) = &coherence_intervals;
  def("coherence_intervals", expt_coherence_intervals_v1);

  MATRIX (*expt_coherence_intervals_v2)(CMATRIX& Coeff, vector<MATRIX>& rates) = &coherence_intervals;
  def("coherence_intervals", expt_coherence_intervals_v2);


  vector<MATRIX> (*expt_schwartz_1_v1)
  (dyn_control_params& prms, CMATRIX& amplitudes, nHamiltonian& ham, MATRIX& inv_alp) = &schwartz_1;

  vector<MATRIX> (*expt_schwartz_2_v1)
  (dyn_control_params& prms, CMATRIX& amplitudes, nHamiltonian& ham, MATRIX& inv_alp) = &schwartz_2;



  ///================== In dyn_methods_dish.cpp  =======================

  vector<int> (*expt_dish_v1)
  (dyn_control_params& prms, MATRIX& q, MATRIX& p,  MATRIX& invM, CMATRIX& Coeff, 
  /*vector<CMATRIX>& projectors,*/ nHamiltonian& ham, vector<int>& act_states, 
  MATRIX& coherence_time, vector<MATRIX>& decoherence_rates, Random& rnd) = &dish;
  def("dish", expt_dish_v1);


  CMATRIX (*expt_afssh_dzdt_v1)
  (CMATRIX& dz, CMATRIX& Hvib, CMATRIX& F, CMATRIX& C, double mass, int act_state) = &afssh_dzdt;
  def("afssh_dzdt", expt_afssh_dzdt_v1);

  void (*expt_integrate_afssh_moments_v1)
  (CMATRIX& dR, CMATRIX& dP, CMATRIX& Hvib, CMATRIX& F, CMATRIX& C, 
  double mass, int act_state, double dt, int nsteps) = &integrate_afssh_moments;
  def("integrate_afssh_moments", expt_integrate_afssh_moments_v1);

  ///================== In dyn_methods_qtsh.cpp ===========================
  MATRIX (*expt_compute_dkinemat_v1)
  (dyn_variables& dyn_var, nHamiltonian& ham) = &compute_dkinemat; 
  def("compute_dkinemat", expt_compute_dkinemat_v1);


}


void export_dyn_hop_acceptance_objects(){

  //============= dyn_hop_proposal.cpp ======================

  int (*expt_can_rescale_along_vector_v1)
  (double E_old, double E_new, MATRIX& p, MATRIX& invM, MATRIX& t, vector<int>& which_dofs) = &can_rescale_along_vector;
  def("can_rescale_along_vector", expt_can_rescale_along_vector_v1);

  int (*expt_can_rescale_along_vector_v2)
  (double E_old, double E_new, MATRIX& p, MATRIX& invM, MATRIX& t) = &can_rescale_along_vector;
  def("can_rescale_along_vector", expt_can_rescale_along_vector_v2);



  void (*expt_rescale_along_vector_v1)
  (double E_old, double E_new, MATRIX& p, MATRIX& invM, MATRIX& t, int do_reverse, vector<int>& which_dofs) = &rescale_along_vector;
  def("rescale_along_vector", expt_rescale_along_vector_v1);

  void (*expt_rescale_along_vector_v2)
  (double E_old, double E_new, MATRIX& p, MATRIX& invM, MATRIX& t, int do_reverse) = &rescale_along_vector;
  def("rescale_along_vector", expt_rescale_along_vector_v2);


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
   MATRIX& q, MATRIX& p, MATRIX& invM, CMATRIX& C, /*vector<CMATRIX>& projectors, */
   nHamiltonian& ham, vector<int>& proposed_states, vector<int>& initial_states, Random& rnd, 
   vector<int>& which_trajectories) = &accept_hops;
  def("accept_hops", expt_accept_hops_v1);

  vector<int> (*expt_accept_hops_v2)
  (dyn_control_params& prms,
   MATRIX& q, MATRIX& p, MATRIX& invM, CMATRIX& C, /*vector<CMATRIX>& projectors, */
   nHamiltonian& ham, vector<int>& proposed_states, vector<int>& initial_states, Random& rnd ) = &accept_hops;
  def("accept_hops", expt_accept_hops_v2);



  vector<int> (*expt_where_can_we_hop_v1)
  (int traj, dyn_control_params& prms,
   MATRIX& q, MATRIX& p,  MATRIX& invM, CMATRIX& Coeff, /*vector<CMATRIX>& projectors, */
   nHamiltonian& ham, vector<int>& act_states, Random& rnd) = &where_can_we_hop;
  def("where_can_we_hop", expt_where_can_we_hop_v1);


  void (*expt_handle_hops_nuclear_v1)
  (dyn_control_params& prms,
   MATRIX& q, MATRIX& p, MATRIX& invM, CMATRIX& C, /*vector<CMATRIX>& projectors,*/
   nHamiltonian& ham, vector<int>& new_states, vector<int>& old_states) = &handle_hops_nuclear;
  def("handle_hops_nuclear", expt_handle_hops_nuclear_v1);

}


void export_dyn_hop_proposal_objects(){

  //============= dyn_hop_proposal.cpp ======================

  MATRIX (*expt_hopping_probabilities_fssh_v1)
  (dyn_control_params& prms, CMATRIX& Coeff, CMATRIX& Hvib) = &hopping_probabilities_fssh;
  def("hopping_probabilities_fssh", expt_hopping_probabilities_fssh_v1);
  vector<double> (*expt_hopping_probabilities_fssh_v2)
  (dyn_control_params& prms, CMATRIX& denmat, CMATRIX& Hvib, int act_state_indx) = &hopping_probabilities_fssh;
  def("hopping_probabilities_fssh", expt_hopping_probabilities_fssh_v2);



  MATRIX (*expt_hopping_probabilities_gfsh_v1)
  (dyn_control_params& prms, CMATRIX& Coeff, CMATRIX& Hvib) = &hopping_probabilities_gfsh;
  def("hopping_probabilities_gfsh", expt_hopping_probabilities_gfsh_v1);
  vector<double> (*expt_hopping_probabilities_gfsh_v2)
  (dyn_control_params& prms, CMATRIX& denmat, CMATRIX& Hvib, int act_state_indx) = &hopping_probabilities_gfsh;
  def("hopping_probabilities_gfsh", expt_hopping_probabilities_gfsh_v2);


  vector<double> (*expt_hopping_probabilities_fssh2_v1)
  (dyn_control_params& prms, CMATRIX& denmat, CMATRIX& denmat_old, int act_state_indx) = &hopping_probabilities_fssh2;
  def("hopping_probabilities_fssh2", expt_hopping_probabilities_fssh2_v1);


  MATRIX (*expt_hopping_probabilities_mssh_v1)
  (dyn_control_params& prms, CMATRIX& Coeff, CMATRIX& Hvib) = &hopping_probabilities_mssh;
  def("hopping_probabilities_mssh", expt_hopping_probabilities_mssh_v1);
  vector<double> (*expt_hopping_probabilities_mssh_v2)
  (dyn_control_params& prms, CMATRIX& denmat, CMATRIX& Hvib, int act_state_indx) = &hopping_probabilities_mssh;
  def("hopping_probabilities_mssh", expt_hopping_probabilities_mssh_v2);


  vector<double> (*expt_hopping_probabilities_lz_v1)
  (nHamiltonian& ham, nHamiltonian& ham_prev, int act_state_indx, int rep, 
  MATRIX& p, const MATRIX& invM) = &hopping_probabilities_lz;
  def("hopping_probabilities_lz", expt_hopping_probabilities_lz_v1);

  vector<double> (*expt_hopping_probabilities_zn_v1)
  (nHamiltonian& ham, nHamiltonian& ham_prev, int act_state_indx, int rep,
  MATRIX& p, const MATRIX& invM) = &hopping_probabilities_lz;
  def("hopping_probabilities_zn", expt_hopping_probabilities_zn_v1);



  vector<double> (*expt_hopping_probabilities_mash_v1)
  (dyn_control_params& prms, CMATRIX& denmat) = &hopping_probabilities_mash;
  def("hopping_probabilities_mash", expt_hopping_probabilities_mash_v1);



  vector<MATRIX> (*expt_hop_proposal_probabilities_v1)
  (dyn_control_params& prms,
   MATRIX& q, MATRIX& p, MATRIX& invM, CMATRIX& C,
   nHamiltonian& ham, vector<MATRIX>& prev_ham_dia) = &hop_proposal_probabilities;
  def("hop_proposal_probabilities", expt_hop_proposal_probabilities_v1);

  vector< vector<double> > (*expt_hop_proposal_probabilities_v2)
  (dyn_control_params& prms, dyn_variables& dyn_var,
   nHamiltonian& ham, nHamiltonian& ham_prev) = &hop_proposal_probabilities;
  def("hop_proposal_probabilities", expt_hop_proposal_probabilities_v2);


  int (*expt_hop_v1)(vector<double>& prob, double ksi) = &hop;
  def("hop", expt_hop_v1);
  int (*expt_hop_v2)(int initstate, MATRIX& g, double ksi) = &hop;
  def("hop", expt_hop_v2);
  int (*expt_hop_v3)(int initstate, vector<double>& g, double ksi) = &hop;
  def("hop", expt_hop_v3);

  vector<int> (*expt_propose_hops_v1)
  (vector<MATRIX>& g, vector<int>& act_states, Random& rnd) = &propose_hops;
  def("propose_hops", expt_propose_hops_v1);
  vector<int> (*expt_propose_hops_v2)
  (vector< vector<double> >& g, vector<int>& act_states, Random& rnd) = &propose_hops;
  def("propose_hops", expt_propose_hops_v2);


}


void export_dyn_methods_objects(){

  vector<int> (*expt_decoherence_event_v1)
  (MATRIX& coherence_time, MATRIX& coherence_interval, int decoherence_event_option, Random& rnd) = &decoherence_event;
  def("decoherence_event", expt_decoherence_event_v1);

  vector<int> (*expt_decoherence_event_v2)
  (MATRIX& coherence_time, MATRIX& coherence_interval, Random& rnd) = &decoherence_event;
  def("decoherence_event", expt_decoherence_event_v2);

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

  vector<int> (*expt_permute_states_v1)
  (vector<vector<int> >& perms, vector<int>& act_states) = &permute_states;
  def("permute_states", expt_permute_states_v1);



  void (*expt_update_projectors_v1)(dyn_control_params& prms, vector<CMATRIX>& projectors, 
  vector<CMATRIX>& Eadi, vector<CMATRIX>& St, Random& rnd) = &update_projectors;
  def("update_projectors", expt_update_projectors_v1);  


  vector< vector<int> > (*expt_compute_permutations_v1)
  (dyn_control_params& prms, vector<CMATRIX>& Eadi, vector<CMATRIX>& St, Random& rnd) = &compute_permutations;


  vector<CMATRIX> (*expt_compute_projectors_v1)
  (dyn_control_params& prms, vector<CMATRIX>& Eadi, vector<CMATRIX>& St, Random& rnd) = &compute_projectors;
  def("compute_projectors", expt_compute_projectors_v1);

  vector<CMATRIX> (*expt_compute_projectors_v2)
  (dyn_control_params& prms, vector<CMATRIX>& St, vector<vector<int> >& perms) = &compute_projectors;
  def("compute_projectors", expt_compute_projectors_v2);


  CMATRIX (*expt_raw_to_dynconsyst_v1)
  (CMATRIX& amplitudes, vector<CMATRIX>& projectors) = &raw_to_dynconsyst;
  def("raw_to_dynconsyst", expt_raw_to_dynconsyst_v1);

  CMATRIX (*expt_dynconsyst_to_raw_v1)
  (CMATRIX& amplitudes, vector<CMATRIX>& projectors) = &dynconsyst_to_raw;
  def("dynconsyst_to_raw", expt_dynconsyst_to_raw_v1);

  vector<int> (*expt_get_stochastic_reordering_v1)
  (CMATRIX& time_overlap, Random& rnd) = &get_stochastic_reordering;
  def("get_stochastic_reordering", expt_get_stochastic_reordering_v1);

  vector<int> (*expt_get_stochastic_reordering2_v1)
  (CMATRIX& time_overlap, Random& rnd) = &get_stochastic_reordering2;
  def("get_stochastic_reordering2", expt_get_stochastic_reordering2_v1);

  vector<int> (*expt_get_stochastic_reordering3_v1)
  (CMATRIX& time_overlap, Random& rnd, int convergence, int max_number_attempts) = &get_stochastic_reordering3;
  def("get_stochastic_reordering3", expt_get_stochastic_reordering3_v1);

  vector<int> (*expt_get_stochastic_reordering3_v2)
  (CMATRIX& time_overlap, Random& rnd, int convergence, int max_number_attempts,
  double filter_tol, int verbosity_level) = &get_stochastic_reordering3;
  def("get_stochastic_reordering3", expt_get_stochastic_reordering3_v2);



}

/*
void export_LZ_hopping_probabilities_objects(){


  MATRIX (*expt_compute_hopping_probabilities_lz_v1)
  (nHamiltonian& ham, int rep, MATRIX& p, const MATRIX& invM, MATRIX& prev_ham_dia) = &compute_hopping_probabilities_lz;

  def("compute_hopping_probabilities_lz",expt_compute_hopping_probabilities_lz_v1);

}
*/


void export_permutation_objects(){


  vector<int> (*expt_get_permutation_v1)(vector<vector<int> >& inp) = &get_permutation; 
  vector<int> (*expt_Munkres_Kuhn_minimize_v1)(MATRIX& _X, int verbosity) = &Munkres_Kuhn_minimize;
  vector<int> (*expt_Munkres_Kuhn_maximize_v1)(MATRIX& _X, int verbosity) = &Munkres_Kuhn_maximize;

  def("get_permutation", expt_get_permutation_v1);  
  def("Munkres_Kuhn_minimize", expt_Munkres_Kuhn_minimize_v1);  
  def("Munkres_Kuhn_maximize", expt_Munkres_Kuhn_maximize_v1);  


}

void export_Energy_Forces_objects(){


  double (*expt_compute_kinetic_energy_v11)(MATRIX& p, MATRIX& invM, vector<int>& which_dofs) = &compute_kinetic_energy;
  def("compute_kinetic_energy",expt_compute_kinetic_energy_v11);
  double (*expt_compute_kinetic_energy_v12)(MATRIX& p, MATRIX& invM) = &compute_kinetic_energy;
  def("compute_kinetic_energy",expt_compute_kinetic_energy_v12);

  vector<double> (*expt_compute_kinetic_energies_v1)(MATRIX& p, MATRIX& invM, vector<int>& which_dofs) = &compute_kinetic_energies;
  def("compute_kinetic_energies",expt_compute_kinetic_energies_v1);

  vector<double> (*expt_compute_kinetic_energies_v2)(MATRIX& p, MATRIX& invM) = &compute_kinetic_energies;
  def("compute_kinetic_energies",expt_compute_kinetic_energies_v2);


  CMATRIX (*expt_tsh_indx2ampl_v1)(vector<int>& res, int nstates) = &tsh_indx2ampl;
  def("tsh_indx2ampl", expt_tsh_indx2ampl_v1);


  double (*expt_average_potential_energy_v1)
  (dyn_control_params& prms, dyn_variables& dyn_vars, nHamiltonian& ham) = &average_potential_energy;
  def("average_potential_energy", expt_average_potential_energy_v1);

  double (*expt_average_potential_energy_v2)
  (bp::dict prms, dyn_variables& dyn_vars, nHamiltonian& ham) = &average_potential_energy;
  def("average_potential_energy", expt_average_potential_energy_v2);

  vector<double> (*expt_potential_energies_v1)
  (dyn_control_params& prms, dyn_variables& dyn_vars, nHamiltonian& ham) = &potential_energies;
  def("potential_energies", expt_potential_energies_v1);

  vector<double> (*expt_potential_energies_v2)
  (bp::dict prms, dyn_variables& dyn_vars, nHamiltonian& ham) = &potential_energies;
  def("potential_energies", expt_potential_energies_v2);


  void (*expt_update_forces_v1)
  (dyn_control_params& prms, dyn_variables& dynvars, nHamiltonian& ham) = &update_forces;
  def("update_forces", expt_update_forces_v1);

  void (*expt_update_forces_v2)
  (bp::dict params, dyn_variables& dynvars, nHamiltonian& ham) = &update_forces;
  def("update_forces", expt_update_forces_v2);



  vector<CMATRIX> (*expt_get_Eadi_v1)(nHamiltonian& ham) = &get_Eadi;
  def("get_Eadi", expt_get_Eadi_v1);


}


void export_dyn_ham(){


  void (*expt_update_Hamiltonian_variables_v1)
  (dyn_control_params& prms, dyn_variables& dyn_var, nHamiltonian& ham, nHamiltonian& ham_prev,
   bp::object py_funct, bp::object model_params, int update_type) = &update_Hamiltonian_variables;

  void (*expt_update_Hamiltonian_variables_v2)
  (bp::dict prms, dyn_variables& dyn_var, nHamiltonian& ham, nHamiltonian& ham_prev,
   bp::object py_funct, bp::object model_params, int update_type) = &update_Hamiltonian_variables;
  
  def("update_Hamiltonian_variables", expt_update_Hamiltonian_variables_v1);
  def("update_Hamiltonian_variables", expt_update_Hamiltonian_variables_v2);

/*
  void (*expt_update_Hamiltonian_q_v1)
  (dyn_control_params& prms, MATRIX& q, nHamiltonian& ham, 
   bp::object py_funct, bp::object model_params) = &update_Hamiltonian_q;

  void (*expt_update_Hamiltonian_q_v2)
  (dyn_control_params& prms, dyn_variables& dyn_var, nHamiltonian& ham, 
   bp::object py_funct, bp::object model_params) = &update_Hamiltonian_q;

  void (*expt_update_Hamiltonian_q_v3)
  (bp::dict prms, MATRIX& q, nHamiltonian& ham, 
   bp::object py_funct, bp::object model_params) = &update_Hamiltonian_q;

  void (*expt_update_Hamiltonian_q_v4)
  (bp::dict prms, dyn_variables& dyn_var, nHamiltonian& ham, 
   bp::object py_funct, bp::object model_params) = &update_Hamiltonian_q;

  def("update_Hamiltonian_q", expt_update_Hamiltonian_q_v1);
  def("update_Hamiltonian_q", expt_update_Hamiltonian_q_v2);
  def("update_Hamiltonian_q", expt_update_Hamiltonian_q_v3);
  def("update_Hamiltonian_q", expt_update_Hamiltonian_q_v4);


  void (*expt_update_Hamiltonian_q_ethd_v1)
  (dyn_control_params& prms, MATRIX& q, MATRIX& p, nHamiltonian& ham, 
   bp::object py_funct, bp::object model_params, MATRIX& invM) = &update_Hamiltonian_q_ethd;

  void (*expt_update_Hamiltonian_q_ethd_v2)
  (dyn_control_params& prms, dyn_variables& dyn_var, nHamiltonian& ham, bp::object py_funct, bp::object model_params) = &update_Hamiltonian_q_ethd;

  void (*expt_update_Hamiltonian_q_ethd_v3)
  (bp::dict prms, MATRIX& q, MATRIX& p, nHamiltonian& ham, 
   bp::object py_funct, bp::object model_params, MATRIX& invM) = &update_Hamiltonian_q_ethd;

  void (*expt_update_Hamiltonian_q_ethd_v4)
  (bp::dict prms, dyn_variables& dyn_var, nHamiltonian& ham, bp::object py_funct, bp::object model_params) = &update_Hamiltonian_q_ethd;

  def("update_Hamiltonian_q_ethd", expt_update_Hamiltonian_q_ethd_v1);
  def("update_Hamiltonian_q_ethd", expt_update_Hamiltonian_q_ethd_v2);
  def("update_Hamiltonian_q_ethd", expt_update_Hamiltonian_q_ethd_v3);
  def("update_Hamiltonian_q_ethd", expt_update_Hamiltonian_q_ethd_v4);


  void (*expt_update_Hamiltonian_p_v1)
  (dyn_control_params& prms, nHamiltonian& ham, MATRIX& p, MATRIX& invM) = &update_Hamiltonian_p;

  void (*expt_update_Hamiltonian_p_v2)
  (dyn_control_params& prms, dyn_variables& dyn_var, nHamiltonian& ham) = &update_Hamiltonian_p;

  void (*expt_update_Hamiltonian_p_v3)
  (bp::dict prms, nHamiltonian& ham, MATRIX& p, MATRIX& invM) = &update_Hamiltonian_p;

  void (*expt_update_Hamiltonian_p_v4)
  (bp::dict prms, dyn_variables& dyn_var, nHamiltonian& ham) = &update_Hamiltonian_p;

  def("update_Hamiltonian_p", expt_update_Hamiltonian_p_v1);
  def("update_Hamiltonian_p", expt_update_Hamiltonian_p_v2);
  def("update_Hamiltonian_p", expt_update_Hamiltonian_p_v3);
  def("update_Hamiltonian_p", expt_update_Hamiltonian_p_v4);


  void (*expt_update_nacs_v1)(dyn_control_params& prms, nHamiltonian& ham) = &update_nacs;
  def("update_nacs", expt_update_nacs_v1);
*/

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
  //export_Ensemble_objects();
  export_gwp_objects();
  export_heom_objects();
  export_qtag_objects();

  export_dyn_control_params_objects();
  export_dyn_variables_objects();
  export_dyn_decoherence_objects();
  export_dyn_hop_acceptance_objects();
  export_dyn_hop_proposal_objects();
  export_dyn_methods_objects();
  export_dyn_projectors_objects();

  
//  export_LZ_hopping_probabilities_objects();
  export_permutation_objects();

  export_Energy_Forces_objects();

  export_dyn_ham();



  //============= Dynamics.cpp ======================


  vector<CMATRIX> (*expt_compute_St_v1)(nHamiltonian& ham, int isNBRA) = &compute_St;
  def("compute_St", expt_compute_St_v1);
  vector<CMATRIX> (*expt_compute_St_v2)(nHamiltonian& ham) = &compute_St;
  def("compute_St", expt_compute_St_v2);  
  vector<CMATRIX> (*expt_compute_St_v3)(nHamiltonian& ham, nHamiltonian& ham_prev, int isNBRA) = &compute_St;
  def("compute_St", expt_compute_St_v3);
  vector<CMATRIX> (*expt_compute_St_v4)(nHamiltonian& ham, nHamiltonian& ham_prev) = &compute_St;
  def("compute_St", expt_compute_St_v4);


  MATRIX (*expt_momenta_on_excited_states_v1)
  (dyn_variables& dyn_var, nHamiltonian& ham, int itraj) = &momenta_on_excited_states;
  def("momenta_on_excited_states", expt_momenta_on_excited_states_v1);

  void (*expt_SSY_correction_v1)
  (CMATRIX& Ham, dyn_variables& dyn_var, nHamiltonian& ham, int itraj) = &SSY_correction;
  def("SSY_correction", expt_SSY_correction_v1);

  CMATRIX (*expt_Zhu_Liouvillian_v1)(double Etot, CMATRIX& Ham, CMATRIX& rho) = Zhu_Liouvillian;
  def("Zhu_Liouvillian", expt_Zhu_Liouvillian_v1);

  void (*expt_propagate_electronic_v1)
  (dyn_variables& dyn_var, nHamiltonian& ham, nHamiltonian& ham_prev, dyn_control_params& prms) = &propagate_electronic;
  def("propagate_electronic", expt_propagate_electronic_v1);

/*
  void (*expt_compute_dynamics_v1)
  (MATRIX& q, MATRIX& p, MATRIX& invM, CMATRIX& C, vector<CMATRIX>& projectors, vector<int>& act_states,
   nHamiltonian& ham, bp::object py_funct, bp::dict model_params, 
   bp::dict dyn_params, Random& rnd) = &compute_dynamics;
  def("compute_dynamics", expt_compute_dynamics_v1);

  void (*expt_compute_dynamics_v2)
  (MATRIX& q, MATRIX& p, MATRIX& invM, CMATRIX& C, vector<CMATRIX>& projectors, vector<int>& act_states, 
   nHamiltonian& ham, bp::object py_funct, bp::dict& model_params, bp::dict& dyn_params, Random& rnd, 
   vector<Thermostat>& therm) = &compute_dynamics;
  def("compute_dynamics", expt_compute_dynamics_v2);

  void (*expt_compute_dynamics_v3)
  (MATRIX& q, MATRIX& p, MATRIX& invM, CMATRIX& C, vector<CMATRIX>& projectors, vector<int>& act_states, 
   nHamiltonian& ham, bp::object py_funct, bp::dict& model_params, bp::dict& dyn_params, Random& rnd, 
   vector<Thermostat>& therm, dyn_variables&) = &compute_dynamics;
  def("compute_dynamics", expt_compute_dynamics_v3);
*/
  void (*expt_compute_dynamics_v4)
  (dyn_variables& dyn_var, bp::dict dyn_params, nHamiltonian& ham, nHamiltonian& ham_aux, 
   bp::object py_funct, bp::dict model_params, Random& rnd, vector<Thermostat>& therm) = &compute_dynamics;
  def("compute_dynamics", expt_compute_dynamics_v4);






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

