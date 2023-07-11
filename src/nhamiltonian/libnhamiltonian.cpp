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
  \file libhamiltonian.cpp
  \brief The file implements Python export function
    
*/

#if defined(USING_PCH)
#include "../pch.h"
#else
#include <memory> // for std::auto_ptr<>
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#endif 

#include "libnhamiltonian.h"

/// liblibra namespace
namespace liblibra{


using namespace boost::python;

/// libnhamiltonian namespace
namespace libnhamiltonian{


void export_nhamiltonian_objects(){
/** 
  \brief Exporter of the libnhamiltonian_generic classes and functions

*/


//  void (Hamiltonian::*set_params)(boost::python::list) = &Hamiltonian::set_params;

  CMATRIX (nHamiltonian::*expt_get_ovlp_dia_v1)() = &nHamiltonian::get_ovlp_dia;
  CMATRIX (nHamiltonian::*expt_get_ovlp_dia_v2)(vector<int>&) = &nHamiltonian::get_ovlp_dia;

  CMATRIX (nHamiltonian::*expt_get_dc1_dia_v1)(int i) = &nHamiltonian::get_dc1_dia;
  CMATRIX (nHamiltonian::*expt_get_dc1_dia_v2)(int i, vector<int>&) = &nHamiltonian::get_dc1_dia;

  CMATRIX (nHamiltonian::*expt_get_ham_dia_v1)() = &nHamiltonian::get_ham_dia;
  CMATRIX (nHamiltonian::*expt_get_ham_dia_v2)(vector<int>&) = &nHamiltonian::get_ham_dia;

  CMATRIX (nHamiltonian::*expt_get_nac_dia_v1)() = &nHamiltonian::get_nac_dia;
  CMATRIX (nHamiltonian::*expt_get_nac_dia_v2)(vector<int>&) = &nHamiltonian::get_nac_dia;

  CMATRIX (nHamiltonian::*expt_get_hvib_dia_v1)() = &nHamiltonian::get_hvib_dia;
  CMATRIX (nHamiltonian::*expt_get_hvib_dia_v2)(vector<int>&) = &nHamiltonian::get_hvib_dia;

  CMATRIX (nHamiltonian::*expt_get_d1ham_dia_v1)(int i) = &nHamiltonian::get_d1ham_dia;
  CMATRIX (nHamiltonian::*expt_get_d1ham_dia_v2)(int i, vector<int>&) = &nHamiltonian::get_d1ham_dia;

  CMATRIX (nHamiltonian::*expt_get_d2ham_dia_v11)(int i) = &nHamiltonian::get_d2ham_dia;
  CMATRIX (nHamiltonian::*expt_get_d2ham_dia_v12)(int i, vector<int>&) = &nHamiltonian::get_d2ham_dia;

  CMATRIX (nHamiltonian::*expt_get_d2ham_dia_v21)(int i, int j) = &nHamiltonian::get_d2ham_dia;
  CMATRIX (nHamiltonian::*expt_get_d2ham_dia_v22)(int i, int j, vector<int>&) = &nHamiltonian::get_d2ham_dia;



  CMATRIX (nHamiltonian::*expt_get_dc1_adi_v1)(int i) = &nHamiltonian::get_dc1_adi;
  CMATRIX (nHamiltonian::*expt_get_dc1_adi_v2)(int i, vector<int>&) = &nHamiltonian::get_dc1_adi;

  CMATRIX (nHamiltonian::*expt_get_ham_adi_v1)() = &nHamiltonian::get_ham_adi;
  CMATRIX (nHamiltonian::*expt_get_ham_adi_v2)(vector<int>&) = &nHamiltonian::get_ham_adi;

  CMATRIX (nHamiltonian::*expt_get_nac_adi_v1)() = &nHamiltonian::get_nac_adi;
  CMATRIX (nHamiltonian::*expt_get_nac_adi_v2)(vector<int>&) = &nHamiltonian::get_nac_adi;

  CMATRIX (nHamiltonian::*expt_get_hvib_adi_v1)() = &nHamiltonian::get_hvib_adi;
  CMATRIX (nHamiltonian::*expt_get_hvib_adi_v2)(vector<int>&) = &nHamiltonian::get_hvib_adi;

  CMATRIX (nHamiltonian::*expt_get_d1ham_adi_v1)(int i) = &nHamiltonian::get_d1ham_adi;
  CMATRIX (nHamiltonian::*expt_get_d1ham_adi_v2)(int i, vector<int>&) = &nHamiltonian::get_d1ham_adi;

  CMATRIX (nHamiltonian::*expt_get_d2ham_adi_v11)(int i) = &nHamiltonian::get_d2ham_adi;
  CMATRIX (nHamiltonian::*expt_get_d2ham_adi_v12)(int i, vector<int>&) = &nHamiltonian::get_d2ham_adi;

  CMATRIX (nHamiltonian::*expt_get_d2ham_adi_v21)(int i, int j) = &nHamiltonian::get_d2ham_adi;
  CMATRIX (nHamiltonian::*expt_get_d2ham_adi_v22)(int i, int j, vector<int>&) = &nHamiltonian::get_d2ham_adi;



  CMATRIX (nHamiltonian::*expt_get_basis_transform_v1)() = &nHamiltonian::get_basis_transform;
  CMATRIX (nHamiltonian::*expt_get_basis_transform_v2)(vector<int>&) = &nHamiltonian::get_basis_transform;

  CMATRIX (nHamiltonian::*expt_get_time_overlap_adi_v1)() = &nHamiltonian::get_time_overlap_adi;
  CMATRIX (nHamiltonian::*expt_get_time_overlap_adi_v2)(vector<int>&) = &nHamiltonian::get_time_overlap_adi;

  CMATRIX (nHamiltonian::*expt_get_time_overlap_dia_v1)() = &nHamiltonian::get_time_overlap_dia;
  CMATRIX (nHamiltonian::*expt_get_time_overlap_dia_v2)(vector<int>&) = &nHamiltonian::get_time_overlap_dia;



  vector<int> (nHamiltonian::*expt_get_ordering_adi_v1)() = &nHamiltonian::get_ordering_adi; 
  vector<int> (nHamiltonian::*expt_get_ordering_adi_v2)(vector<int>&) = &nHamiltonian::get_ordering_adi; 

  CMATRIX (nHamiltonian::*expt_get_cum_phase_corr_v1)() = &nHamiltonian::get_cum_phase_corr;
  CMATRIX (nHamiltonian::*expt_get_cum_phase_corr_v2)(vector<int>&) = &nHamiltonian::get_cum_phase_corr;



  // for internal model types
  void (nHamiltonian::*expt_compute_diabatic_v1)(int model, vector<double>& q, vector<double>& params, int lvl)
  = &nHamiltonian::compute_diabatic; 

  void (nHamiltonian::*expt_compute_diabatic_v2)(int model, vector<double>& q, vector<double>& params)
  = &nHamiltonian::compute_diabatic; 


  // for models defined in Python
/*
  void (nHamiltonian::*expt_compute_diabatic_v3)(bp::object py_funct, bp::object q, bp::object params, int lvl)
  = &nHamiltonian::compute_diabatic;

  void (nHamiltonian::*expt_compute_diabatic_v4)(bp::object py_funct, bp::object q, bp::object params)
  = &nHamiltonian::compute_diabatic;
*/
  void (nHamiltonian::*expt_compute_diabatic_v3)(bp::object py_funct, MATRIX& q, bp::object params, int lvl)
  = &nHamiltonian::compute_diabatic;

  void (nHamiltonian::*expt_compute_diabatic_v4)(bp::object py_funct, MATRIX& q, bp::object params)
  = &nHamiltonian::compute_diabatic;



  void (nHamiltonian::*expt_update_ordering_v1)(vector<int>& perm_t, int lvl) = &nHamiltonian::update_ordering;
  void (nHamiltonian::*expt_update_ordering_v2)(vector<int>& perm_t) = &nHamiltonian::update_ordering;
/*
  void (nHamiltonian::*expt_apply_phase_corrections_v1)(CMATRIX& phase_corr, int lvl) = &nHamiltonian::apply_phase_corrections;
  void (nHamiltonian::*expt_apply_phase_corrections_v2)(CMATRIX& phase_corr) = &nHamiltonian::apply_phase_corrections;

  CMATRIX (nHamiltonian::*expt_update_phases_v1)(CMATRIX& phases, int lvl) = &nHamiltonian::update_phases;
  CMATRIX (nHamiltonian::*expt_update_phases_v2)(CMATRIX& phases) = &nHamiltonian::update_phases;
*/




  void (nHamiltonian::*expt_compute_adiabatic_v1)(int der_lvl, int lvl)
  = &nHamiltonian::compute_adiabatic;
  void (nHamiltonian::*expt_compute_adiabatic_v2)(int der_lvl)
  = &nHamiltonian::compute_adiabatic;
  void (nHamiltonian::*expt_compute_adiabatic_v3)(bp::object py_funct, MATRIX& q, bp::object params, int lvl)
  = &nHamiltonian::compute_adiabatic;
  void (nHamiltonian::*expt_compute_adiabatic_v4)(bp::object py_funct, MATRIX& q, bp::object params)
  = &nHamiltonian::compute_adiabatic;

/*
  void (nHamiltonian::*expt_compute_adiabatic_v3)(bp::object py_funct, bp::object q, bp::object params, int lvl)
  = &nHamiltonian::compute_adiabatic;
  void (nHamiltonian::*expt_compute_adiabatic_v4)(bp::object py_funct, bp::object q, bp::object params)
  = &nHamiltonian::compute_adiabatic;
*/


  void (nHamiltonian::*expt_ampl_dia2adi_v1)(CMATRIX& ampl_dia, CMATRIX& ampl_adi) 
  = &nHamiltonian::ampl_dia2adi;
  void (nHamiltonian::*expt_ampl_dia2adi_v2)(CMATRIX& ampl_dia, CMATRIX& ampl_adi, vector<int>& id_)
  = &nHamiltonian::ampl_dia2adi;
  void (nHamiltonian::*expt_ampl_dia2adi_v3)(CMATRIX& ampl_dia, CMATRIX& ampl_adi, int lvl, int split)
  = &nHamiltonian::ampl_dia2adi;
  void (nHamiltonian::*expt_ampl_adi2dia_v1)(CMATRIX& ampl_dia, CMATRIX& ampl_adi)
  = &nHamiltonian::ampl_adi2dia;
  void (nHamiltonian::*expt_ampl_adi2dia_v2)(CMATRIX& ampl_dia, CMATRIX& ampl_adi, vector<int>& id_)
  = &nHamiltonian::ampl_adi2dia;
  void (nHamiltonian::*expt_ampl_adi2dia_v3)(CMATRIX& ampl_dia, CMATRIX& ampl_adi, int lvl, int split)
  = &nHamiltonian::ampl_adi2dia;




  void (nHamiltonian::*expt_compute_nac_dia_v1)(MATRIX& p, const MATRIX& invM)
  = &nHamiltonian::compute_nac_dia;
  void (nHamiltonian::*expt_compute_nac_dia_v2)(MATRIX& p, const MATRIX& invM, vector<int>& id_)
  = &nHamiltonian::compute_nac_dia;
  void (nHamiltonian::*expt_compute_nac_dia_v3)(MATRIX& p, const MATRIX& invM, int lvl, int split)
  = &nHamiltonian::compute_nac_dia;
  void (nHamiltonian::*expt_compute_nac_adi_v1)(MATRIX& p, const MATRIX& invM)
  = &nHamiltonian::compute_nac_adi;
  void (nHamiltonian::*expt_compute_nac_adi_v2)(MATRIX& p, const MATRIX& invM, vector<int>& id_)
  = &nHamiltonian::compute_nac_adi;
  void (nHamiltonian::*expt_compute_nac_adi_v3)(MATRIX& p, const MATRIX& invM, int lvl, int split)
  = &nHamiltonian::compute_nac_adi;


  void (nHamiltonian::*expt_compute_hvib_dia_v1)()
  = &nHamiltonian::compute_hvib_dia;
  void (nHamiltonian::*expt_compute_hvib_dia_v2)(vector<int>& id_)
  = &nHamiltonian::compute_hvib_dia;
  void (nHamiltonian::*expt_compute_hvib_dia_v3)(int lvl)
  = &nHamiltonian::compute_hvib_dia;
  void (nHamiltonian::*expt_compute_hvib_adi_v1)()
  = &nHamiltonian::compute_hvib_adi;
  void (nHamiltonian::*expt_compute_hvib_adi_v2)(vector<int>& id_)
  = &nHamiltonian::compute_hvib_adi;
  void (nHamiltonian::*expt_compute_hvib_adi_v3)(int lvl)
  = &nHamiltonian::compute_hvib_adi;



  CMATRIX (nHamiltonian::*expt_forces_adi_v1)(vector<int>& act_states)
  = &nHamiltonian::forces_adi;

  CMATRIX (nHamiltonian::*expt_forces_dia_v1)(vector<int>& act_states)
  = &nHamiltonian::forces_dia;

  CMATRIX (nHamiltonian::*expt_all_forces_adi_v1)(vector<int>& id_) 
  = &nHamiltonian::all_forces_adi;


/*
  CMATRIX (nHamiltonian::*expt_forces_adi_v1)(CMATRIX& ampl_adi) 
  = &nHamiltonian::forces_adi;
  CMATRIX (nHamiltonian::*expt_forces_adi_v2)(CMATRIX& ampl_adi, vector<int>& id_)
  = &nHamiltonian::forces_adi;  
  CMATRIX (nHamiltonian::*expt_forces_adi_v3)(vector<int>& act_states)
  = &nHamiltonian::forces_adi;
  CMATRIX (nHamiltonian::*expt_forces_adi_v4)(int act_state)
  = &nHamiltonian::forces_adi;


  CMATRIX (nHamiltonian::*expt_forces_dia_v1)(CMATRIX& ampl_dia)
  = &nHamiltonian::forces_dia;  
  CMATRIX (nHamiltonian::*expt_forces_dia_v2)(CMATRIX& ampl_dia, vector<int>& id_)
  = &nHamiltonian::forces_dia;
  CMATRIX (nHamiltonian::*expt_forces_dia_v3)(vector<int>& act_states)
  = &nHamiltonian::forces_dia;
  CMATRIX (nHamiltonian::*expt_forces_dia_v4)(int act_state)
  = &nHamiltonian::forces_dia;


  vector<CMATRIX> (nHamiltonian::*expt_forces_tens_adi_v1)(CMATRIX& ampl_adi)
  = &nHamiltonian::forces_tens_adi; 
  vector<CMATRIX> (nHamiltonian::*expt_forces_tens_adi_v2)(CMATRIX& ampl_adi, vector<int>& id_)
  = &nHamiltonian::forces_tens_adi; 
  vector<CMATRIX> (nHamiltonian::*expt_forces_tens_dia_v1)(CMATRIX& ampl_dia)
  = &nHamiltonian::forces_tens_dia; 
  vector<CMATRIX> (nHamiltonian::*expt_forces_tens_dia_v2)(CMATRIX& ampl_dia, vector<int>& id_) 
  = &nHamiltonian::forces_tens_dia; 

*/


  complex<double> (nHamiltonian::*expt_Ehrenfest_energy_adi_v1)(CMATRIX& ampl_adi)
  = &nHamiltonian::Ehrenfest_energy_adi;
  complex<double> (nHamiltonian::*expt_Ehrenfest_energy_adi_v2)(CMATRIX& ampl_adi, vector<int>& id_)
  = &nHamiltonian::Ehrenfest_energy_adi;
  complex<double> (nHamiltonian::*expt_Ehrenfest_energy_dia_v1)(CMATRIX& ampl_dia) 
  = &nHamiltonian::Ehrenfest_energy_dia;
  complex<double> (nHamiltonian::*expt_Ehrenfest_energy_dia_v2)(CMATRIX& ampl_dia, vector<int>& id_)
  = &nHamiltonian::Ehrenfest_energy_dia;


  CMATRIX (nHamiltonian::*expt_Ehrenfest_forces_adi_v1)(CMATRIX& ampl_adi, int lvl, int option)
  = &nHamiltonian::Ehrenfest_forces_adi;
  CMATRIX (nHamiltonian::*expt_Ehrenfest_forces_adi_v2)(CMATRIX& ampl_adi, int lvl)
  = &nHamiltonian::Ehrenfest_forces_adi;
//  CMATRIX (nHamiltonian::*expt_Ehrenfest_forces_adi_v2)(CMATRIX& ampl_adi, vector<int>& id_)
//  = &nHamiltonian::Ehrenfest_forces_adi;
  CMATRIX (nHamiltonian::*expt_Ehrenfest_forces_dia_v1)(CMATRIX& ampl_dia, int lvl, int option)
  = &nHamiltonian::Ehrenfest_forces_dia;
  CMATRIX (nHamiltonian::*expt_Ehrenfest_forces_dia_v2)(CMATRIX& ampl_dia, int lvl)
  = &nHamiltonian::Ehrenfest_forces_dia;
//  CMATRIX (nHamiltonian::*expt_Ehrenfest_forces_dia_v2)(CMATRIX& ampl_dia, vector<int>& id_)
//  = &nHamiltonian::Ehrenfest_forces_dia;

  vector<CMATRIX> (nHamiltonian::*expt_Ehrenfest_forces_tens_adi_v1)(CMATRIX& ampl_adi)
  = &nHamiltonian::Ehrenfest_forces_tens_adi; 
  vector<CMATRIX> (nHamiltonian::*expt_Ehrenfest_forces_tens_adi_v2)(CMATRIX& ampl_adi, vector<int>& id_)
  = &nHamiltonian::Ehrenfest_forces_tens_adi;
  vector<CMATRIX> (nHamiltonian::*expt_Ehrenfest_forces_tens_dia_v1)(CMATRIX& ampl_dia)
  = &nHamiltonian::Ehrenfest_forces_tens_dia;
  vector<CMATRIX> (nHamiltonian::*expt_Ehrenfest_forces_tens_dia_v2)(CMATRIX& ampl_dia, vector<int>& id_)
  = &nHamiltonian::Ehrenfest_forces_tens_dia;


  void (nHamiltonian::*expt_add_ethd_dia_v1)(const MATRIX& q, const MATRIX& invM, int der_lvl) = &nHamiltonian::add_ethd_dia;
  void (nHamiltonian::*expt_add_ethd_adi_v1)(const MATRIX& q, const MATRIX& invM, int der_lvl) = &nHamiltonian::add_ethd_adi;

  void (nHamiltonian::*expt_add_ethd3_dia_v1)(const MATRIX& q, const MATRIX& invM, double alp, int der_lvl) = &nHamiltonian::add_ethd3_dia;
  void (nHamiltonian::*expt_add_ethd3_adi_v1)(const MATRIX& q, const MATRIX& invM, double alp, int der_lvl) = &nHamiltonian::add_ethd3_adi;
  void (nHamiltonian::*expt_add_ethd3_adi_v2)(const MATRIX& q, const MATRIX& p, const MATRIX& invM, double alp, double bet, int der_lvl) 
  = &nHamiltonian::add_ethd3_adi;


  void (nHamiltonian::*expt_init_all_v1)(int der_lvl) = &nHamiltonian::init_all;
  void (nHamiltonian::*expt_init_all_v2)(int der_lvl, int lvl) = &nHamiltonian::init_all;


  void (nHamiltonian::*expt_copy_content_v1)(const nHamiltonian& src) = &nHamiltonian::copy_content;


  void (nHamiltonian::*expt_show_memory_status_v1)(vector<int>& id_) = &nHamiltonian::show_memory_status;


  class_<nHamiltonian>("nHamiltonian",init<int,int,int>())
      .def(init<const nHamiltonian&>())
//      .def("__copy__", &generic__copy__<Hamiltonian>)
//      .def("__deepcopy__", &generic__deepcopy__<Hamiltonian>)
      .def_readwrite("id", &nHamiltonian::id)
      .def_readwrite("level", &nHamiltonian::level)
      .def_readwrite("ndia", &nHamiltonian::ndia)
      .def_readwrite("nadi", &nHamiltonian::nadi)
      .def_readwrite("nnucl", &nHamiltonian::nnucl)

      .def_readwrite("eigen_algo", &nHamiltonian::eigen_algo)
      .def_readwrite("phase_corr_ovlp_tol", &nHamiltonian::phase_corr_ovlp_tol)
  
      .def("set_levels", &nHamiltonian::set_levels)
      .def("add_child", &nHamiltonian::add_child)
      .def("add_new_children", &nHamiltonian::add_new_children)
      .def("get_full_id", &nHamiltonian::get_full_id)

      .def("copy_content", expt_copy_content_v1)
      .def("init_all", expt_init_all_v1)
      .def("init_all", expt_init_all_v2)
      .def("show_memory_status", expt_show_memory_status_v1)

      .def("init_ovlp_dia", &nHamiltonian::init_ovlp_dia)
      .def("set_ovlp_dia_by_ref", &nHamiltonian::set_ovlp_dia_by_ref)
      .def("set_ovlp_dia_by_val", &nHamiltonian::set_ovlp_dia_by_val)

      .def("init_dc1_dia", &nHamiltonian::init_dc1_dia)
      .def("set_dc1_dia_by_ref", &nHamiltonian::set_dc1_dia_by_ref)
      .def("set_dc1_dia_by_val", &nHamiltonian::set_dc1_dia_by_val)

      .def("init_ham_dia", &nHamiltonian::init_ham_dia)
      .def("set_ham_dia_by_ref", &nHamiltonian::set_ham_dia_by_ref)
      .def("set_ham_dia_by_val", &nHamiltonian::set_ham_dia_by_val)

      .def("init_nac_dia", &nHamiltonian::init_nac_dia)
      .def("set_nac_dia_by_ref", &nHamiltonian::set_nac_dia_by_ref)
      .def("set_nac_dia_by_val", &nHamiltonian::set_nac_dia_by_val)

      .def("init_hvib_dia", &nHamiltonian::init_hvib_dia)
      .def("set_hvib_dia_by_ref", &nHamiltonian::set_hvib_dia_by_ref)
      .def("set_hvib_dia_by_val", &nHamiltonian::set_hvib_dia_by_val)

      .def("init_d1ham_dia", &nHamiltonian::init_d1ham_dia)
      .def("set_d1ham_dia_by_ref", &nHamiltonian::set_d1ham_dia_by_ref)
      .def("set_d1ham_dia_by_val", &nHamiltonian::set_d1ham_dia_by_val)

      .def("init_d2ham_dia", &nHamiltonian::init_d2ham_dia)
      .def("set_d2ham_dia_by_ref", &nHamiltonian::set_d2ham_dia_by_ref)
      .def("set_d2ham_dia_by_val", &nHamiltonian::set_d2ham_dia_by_val)


      .def("init_dc1_adi", &nHamiltonian::init_dc1_adi)
      .def("set_dc1_adi_by_ref", &nHamiltonian::set_dc1_adi_by_ref)
      .def("set_dc1_adi_by_val", &nHamiltonian::set_dc1_adi_by_val)

      .def("init_ham_adi", &nHamiltonian::init_ham_adi)
      .def("set_ham_adi_by_ref", &nHamiltonian::set_ham_adi_by_ref)
      .def("set_ham_adi_by_val", &nHamiltonian::set_ham_adi_by_val)

      .def("init_nac_adi", &nHamiltonian::init_nac_adi)
      .def("set_nac_adi_by_ref", &nHamiltonian::set_nac_adi_by_ref)
      .def("set_nac_adi_by_val", &nHamiltonian::set_nac_adi_by_val)

      .def("init_hvib_adi", &nHamiltonian::init_hvib_adi)
      .def("set_hvib_adi_by_ref", &nHamiltonian::set_hvib_adi_by_ref)
      .def("set_hvib_adi_by_val", &nHamiltonian::set_hvib_adi_by_val)

      .def("init_d1ham_adi", &nHamiltonian::init_d1ham_adi)
      .def("set_d1ham_adi_by_ref", &nHamiltonian::set_d1ham_adi_by_ref)
      .def("set_d1ham_adi_by_val", &nHamiltonian::set_d1ham_adi_by_val)

      .def("init_d2ham_adi", &nHamiltonian::init_d2ham_adi)
      .def("set_d2ham_adi_by_ref", &nHamiltonian::set_d2ham_adi_by_ref)
      .def("set_d2ham_adi_by_val", &nHamiltonian::set_d2ham_adi_by_val)

      .def("init_basis_transform", &nHamiltonian::init_basis_transform)
      .def("set_basis_transform_by_ref", &nHamiltonian::set_basis_transform_by_ref)
      .def("set_basis_transform_by_val", &nHamiltonian::set_basis_transform_by_val)

      .def("init_time_overlap_adi", &nHamiltonian::init_time_overlap_adi)
      .def("set_time_overlap_adi_by_ref", &nHamiltonian::set_time_overlap_adi_by_ref)
      .def("set_time_overlap_adi_by_val", &nHamiltonian::set_time_overlap_adi_by_val)

      .def("init_time_overlap_dia", &nHamiltonian::init_time_overlap_dia)
      .def("set_time_overlap_dia_by_ref", &nHamiltonian::set_time_overlap_dia_by_ref)
      .def("set_time_overlap_dia_by_val", &nHamiltonian::set_time_overlap_dia_by_val)


      .def("set_ordering_adi_by_ref", &nHamiltonian::set_ordering_adi_by_ref)
      .def("set_ordering_adi_by_val", &nHamiltonian::set_ordering_adi_by_val)

      .def("init_cum_phase_corr", &nHamiltonian::init_cum_phase_corr)
      .def("set_cum_phase_corr_by_ref", &nHamiltonian::set_cum_phase_corr_by_ref)
      .def("set_cum_phase_corr_by_val", &nHamiltonian::set_cum_phase_corr_by_val)



      .def("get_ovlp_dia", expt_get_ovlp_dia_v1)
      .def("get_ovlp_dia", expt_get_ovlp_dia_v2)
      .def("get_dc1_dia", expt_get_dc1_dia_v1)
      .def("get_dc1_dia", expt_get_dc1_dia_v2)
      .def("get_ham_dia", expt_get_ham_dia_v1)
      .def("get_ham_dia", expt_get_ham_dia_v2)
      .def("get_nac_dia", expt_get_nac_dia_v1)
      .def("get_ham_dia", expt_get_ham_dia_v2)
      .def("get_hvib_dia", expt_get_hvib_dia_v1)
      .def("get_hvib_dia", expt_get_hvib_dia_v2)
      .def("get_d1ham_dia", expt_get_d1ham_dia_v1)
      .def("get_d1ham_dia", expt_get_d1ham_dia_v2)
      .def("get_d2ham_dia", expt_get_d2ham_dia_v11)
      .def("get_d2ham_dia", expt_get_d2ham_dia_v12)
      .def("get_d2ham_dia", expt_get_d2ham_dia_v21)
      .def("get_d2ham_dia", expt_get_d2ham_dia_v22)

      .def("get_dc1_adi", expt_get_dc1_adi_v1)
      .def("get_dc1_adi", expt_get_dc1_adi_v2)
      .def("get_ham_adi", expt_get_ham_adi_v1)
      .def("get_ham_adi", expt_get_ham_adi_v2)
      .def("get_nac_adi", expt_get_nac_adi_v1)
      .def("get_ham_adi", expt_get_ham_adi_v2)
      .def("get_hvib_adi", expt_get_hvib_adi_v1)
      .def("get_hvib_adi", expt_get_hvib_adi_v2)
      .def("get_d1ham_adi", expt_get_d1ham_adi_v1)
      .def("get_d1ham_adi", expt_get_d1ham_adi_v2)
      .def("get_d2ham_adi", expt_get_d2ham_adi_v11)
      .def("get_d2ham_adi", expt_get_d2ham_adi_v12)
      .def("get_d2ham_adi", expt_get_d2ham_adi_v21)
      .def("get_d2ham_adi", expt_get_d2ham_adi_v22)

      .def("get_basis_transform", expt_get_basis_transform_v1)
      .def("get_basis_transform", expt_get_basis_transform_v2)
      .def("get_time_overlap_adi", expt_get_time_overlap_adi_v1)
      .def("get_time_overlap_adi", expt_get_time_overlap_adi_v2)
      .def("get_time_overlap_dia", expt_get_time_overlap_dia_v1)
      .def("get_time_overlap_dia", expt_get_time_overlap_dia_v2)
      .def("get_ordering_adi", expt_get_ordering_adi_v1)
      .def("get_ordering_adi", expt_get_ordering_adi_v2)
      .def("get_cum_phase_corr", expt_get_cum_phase_corr_v1)
      .def("get_cum_phase_corr", expt_get_cum_phase_corr_v2)



      .def("compute_diabatic", expt_compute_diabatic_v1)
      .def("compute_diabatic", expt_compute_diabatic_v2)
      .def("compute_diabatic", expt_compute_diabatic_v3)
      .def("compute_diabatic", expt_compute_diabatic_v4)


      .def("update_ordering", expt_update_ordering_v1)
      .def("update_ordering", expt_update_ordering_v2)
//      .def("apply_phase_corrections", expt_apply_phase_corrections_v1)
//      .def("apply_phase_corrections", expt_apply_phase_corrections_v2)
//      .def("update_phases", expt_update_phases_v1)
//      .def("update_phases", expt_update_phases_v2)


      .def("compute_adiabatic", expt_compute_adiabatic_v1)
      .def("compute_adiabatic", expt_compute_adiabatic_v2)
      .def("compute_adiabatic", expt_compute_adiabatic_v3)
      .def("compute_adiabatic", expt_compute_adiabatic_v4)


      .def("ampl_adi2dia", expt_ampl_adi2dia_v1)
      .def("ampl_adi2dia", expt_ampl_adi2dia_v2)
      .def("ampl_adi2dia", expt_ampl_adi2dia_v3)
      .def("ampl_dia2adi", expt_ampl_dia2adi_v1)
      .def("ampl_dia2adi", expt_ampl_dia2adi_v2)
      .def("ampl_dia2adi", expt_ampl_dia2adi_v3)


//      .def("forces_tens_adi", expt_forces_tens_adi_v1)
//      .def("forces_tens_adi", expt_forces_tens_adi_v2)
//      .def("forces_tens_dia", expt_forces_tens_dia_v1)
//      .def("forces_tens_dia", expt_forces_tens_dia_v2)

      .def("forces_adi", expt_forces_adi_v1)
//      .def("forces_adi", expt_forces_adi_v2)
//      .def("forces_adi", expt_forces_adi_v3)
//      .def("forces_adi", expt_forces_adi_v4)

      .def("forces_dia", expt_forces_dia_v1)
//      .def("forces_dia", expt_forces_dia_v2)
//      .def("forces_dia", expt_forces_dia_v3)
//      .def("forces_dia", expt_forces_dia_v4)

      .def("all_forces_adi", expt_all_forces_adi_v1)


      .def("compute_nac_dia", expt_compute_nac_dia_v1)
      .def("compute_nac_dia", expt_compute_nac_dia_v2)
      .def("compute_nac_dia", expt_compute_nac_dia_v3)
      .def("compute_nac_adi", expt_compute_nac_adi_v1)
      .def("compute_nac_adi", expt_compute_nac_adi_v2)
      .def("compute_nac_adi", expt_compute_nac_adi_v3)

      .def("compute_hvib_dia", expt_compute_hvib_dia_v1)
      .def("compute_hvib_dia", expt_compute_hvib_dia_v2)
      .def("compute_hvib_dia", expt_compute_hvib_dia_v3)
      .def("compute_hvib_adi", expt_compute_hvib_adi_v1)
      .def("compute_hvib_adi", expt_compute_hvib_adi_v2)
      .def("compute_hvib_adi", expt_compute_hvib_adi_v3)

      .def("add_ethd_dia", expt_add_ethd_dia_v1)
      .def("add_ethd_adi", expt_add_ethd_adi_v1)

      .def("add_ethd3_dia", expt_add_ethd3_dia_v1)
      .def("add_ethd3_adi", expt_add_ethd3_adi_v1)
      .def("add_ethd3_adi", expt_add_ethd3_adi_v2)




      .def("Ehrenfest_energy_adi", expt_Ehrenfest_energy_adi_v1)
      .def("Ehrenfest_energy_adi", expt_Ehrenfest_energy_adi_v2)
      .def("Ehrenfest_energy_dia", expt_Ehrenfest_energy_dia_v1)
      .def("Ehrenfest_energy_dia", expt_Ehrenfest_energy_dia_v2)

      .def("Ehrenfest_forces_tens_adi", expt_Ehrenfest_forces_tens_adi_v1)
      .def("Ehrenfest_forces_tens_adi", expt_Ehrenfest_forces_tens_adi_v2)
      .def("Ehrenfest_forces_tens_dia", expt_Ehrenfest_forces_tens_dia_v1)
      .def("Ehrenfest_forces_tens_dia", expt_Ehrenfest_forces_tens_dia_v2)

      .def("Ehrenfest_forces_adi", expt_Ehrenfest_forces_adi_v1)
      .def("Ehrenfest_forces_adi", expt_Ehrenfest_forces_adi_v2)
//      .def("Ehrenfest_forces_adi", expt_Ehrenfest_forces_adi_v2)
      .def("Ehrenfest_forces_dia", expt_Ehrenfest_forces_dia_v1)
      .def("Ehrenfest_forces_dia", expt_Ehrenfest_forces_dia_v2)
//      .def("Ehrenfest_forces_dia", expt_Ehrenfest_forces_dia_v2)


  ;

  class_< nHamiltonianList >("nHamiltonianList")
      .def(vector_indexing_suite< nHamiltonianList >())
  ;



/*
  CMATRIX (*expt_compute_phase_corrections1_v1)(CMATRIX& S, double tol) = &compute_phase_corrections1;
  CMATRIX (*expt_compute_phase_corrections1_v2)(CMATRIX& U, CMATRIX& U_prev, double tol) = &compute_phase_corrections1;

  CMATRIX (*expt_compute_phase_corrections_v1)(CMATRIX& S) = &compute_phase_corrections;
  CMATRIX (*expt_compute_phase_corrections_v2)(CMATRIX& U, CMATRIX& U_prev) = &compute_phase_corrections;

  def("compute_phase_corrections1", expt_compute_phase_corrections1_v1);
  def("compute_phase_corrections1", expt_compute_phase_corrections1_v2);
  def("compute_phase_corrections", expt_compute_phase_corrections_v1);
  def("compute_phase_corrections", expt_compute_phase_corrections_v2);
*/


  double (*expt_ETHD_energy_v1)(const MATRIX& q, const MATRIX& invM) = &ETHD_energy;
  MATRIX (*expt_ETHD_forces_v1)(const MATRIX& q, const MATRIX& invM) = &ETHD_forces;

  def("ETHD_energy", expt_ETHD_energy_v1);
  def("ETHD_forces", expt_ETHD_forces_v1);


  double (*expt_ETHD3_energy_v1)(const MATRIX& q, const MATRIX& invM, double alp) = &ETHD3_energy;
  double (*expt_ETHD3_energy_v2)(const MATRIX& q, const MATRIX& p, const MATRIX& invM, double alp, double bet) = &ETHD3_energy;
  MATRIX (*expt_ETHD3_forces_v1)(const MATRIX& q, const MATRIX& invM, double alp) = &ETHD3_forces;
  MATRIX (*expt_ETHD3_forces_v2)(const MATRIX& q, const MATRIX& p, const MATRIX& invM, double alp, double bet) = &ETHD3_forces;
  MATRIX (*expt_ETHD3_friction_v1)(const MATRIX& q, const MATRIX& p, const MATRIX& invM, double alp, double bet) = &ETHD3_friction;

  def("ETHD3_energy", expt_ETHD3_energy_v1);
  def("ETHD3_energy", expt_ETHD3_energy_v2);
  def("ETHD3_forces", expt_ETHD3_forces_v1);
  def("ETHD3_forces", expt_ETHD3_forces_v2);
  def("ETHD3_friction", expt_ETHD3_friction_v1);


}


#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygnhamiltonian){
#else
BOOST_PYTHON_MODULE(libnhamiltonian){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

//  export_Mathematics_objects();
  export_nhamiltonian_objects();

}



}// namespace libnhamiltonian
}// liblibra

