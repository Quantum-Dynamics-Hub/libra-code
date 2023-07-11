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
  \file libensemble.cpp
  \brief The file implements Python export function
    
*/

#if defined(USING_PCH)
#include "../../pch.h"
#else
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#endif 

#include "libensemble.h"

/// liblibra namespace
namespace liblibra{


using namespace boost::python;

/// libdyn namespace 
namespace libdyn{

/// libensemble namespace 
namespace libensemble{


void export_Ensemble_objects(){
/** 
  \brief Exporter of the libensemble classes and functions

*/


//  def("SAC_Ham", expt_SAC_Ham1);
//  void (Electronic::*expt_propagate_electronic)(double,Hamiltonian&) = &Electronic::propagate_electronic;

  void (Ensemble::*expt_ham_set_v_v1)(int i, vector<double>& v) = &Ensemble::ham_set_v;
  void (Ensemble::*expt_ham_set_v_v2)(int i, boost::python::list v) = &Ensemble::ham_set_v;
  void (Ensemble::*expt_ham_set_v_v3)() = &Ensemble::ham_set_v;

  void (Ensemble::*expt_ham_set_q_v1)(int i, vector<double>& q) = &Ensemble::ham_set_q;
  void (Ensemble::*expt_ham_set_q_v2)(int i, boost::python::list q) = &Ensemble::ham_set_q;

  void (Ensemble::*expt_ham_set_ham_v1)(int i, std::string opt, int mopt) = &Ensemble::ham_set_ham;
  void (Ensemble::*expt_ham_set_ham_v2)(std::string opt, int mopt) = &Ensemble::ham_set_ham;
  void (Ensemble::*expt_ham_set_ham_v3)(int i, Hamiltonian& _ham) = &Ensemble::ham_set_ham;

  void (Ensemble::*expt_ham_set_rep_v1)(int i, int _rep) = &Ensemble::ham_set_rep;
  void (Ensemble::*expt_ham_set_rep_v2)(int _rep) = &Ensemble::ham_set_rep;

  void (Ensemble::*expt_ham_set_params_v1)(int i, vector<double>& v) = &Ensemble::ham_set_params;
  void (Ensemble::*expt_ham_set_params_v2)(int i, boost::python::list v) = &Ensemble::ham_set_params;
  void (Ensemble::*expt_ham_set_params_v3)(vector<double>& v) = &Ensemble::ham_set_params;
  void (Ensemble::*expt_ham_set_params_v4)(boost::python::list v) = &Ensemble::ham_set_params;


  void (Ensemble::*expt_ham_compute_v1)(int i) = &Ensemble::ham_compute;
  void (Ensemble::*expt_ham_compute_diabatic_v1)(int i) = &Ensemble::ham_compute_diabatic;
  void (Ensemble::*expt_ham_compute_adiabatic_v1)(int i) = &Ensemble::ham_compute_adiabatic;

  std::complex<double> (Ensemble::*expt_ham_H_v1)(int traj, int i,int j) = &Ensemble::ham_H;
  std::complex<double> (Ensemble::*expt_ham_dHdq_v1)(int traj, int i,int j,int n) = &Ensemble::ham_dHdq;
  std::complex<double> (Ensemble::*expt_ham_D_v1)(int traj, int i,int j,int n) = &Ensemble::ham_D;
  std::complex<double> (Ensemble::*expt_ham_nac_v1)(int traj,int i,int j) = &Ensemble::ham_nac; 
  std::complex<double> (Ensemble::*expt_ham_Hvib_v1)(int traj, int i,int j) = &Ensemble::ham_Hvib;




  void (Ensemble::*expt_el_propagate_electronic_v1)(int i, double dt) = &Ensemble::el_propagate_electronic;
  void (Ensemble::*expt_el_propagate_electronic_v2)(double dt) = &Ensemble::el_propagate_electronic;
  void (Ensemble::*expt_mol_propagate_q_v1)(int i, double dt) = &Ensemble::mol_propagate_q;
  void (Ensemble::*expt_mol_propagate_q_v2)(double dt) = &Ensemble::mol_propagate_q;
  void (Ensemble::*expt_mol_propagate_p_v1)(int i, double dt) = &Ensemble::mol_propagate_p;
  void (Ensemble::*expt_mol_propagate_p_v2)(double dt) = &Ensemble::mol_propagate_p;

  void (Ensemble::*expt_se_pop_v1)(vector<double>& v,double xmax,double xmin) = &Ensemble::se_pop;
  void (Ensemble::*expt_se_pop_v2)(vector<double>& v) = &Ensemble::se_pop;
  boost::python::list (Ensemble::*expt_se_pop_v3)(double xmax,double xmin) = &Ensemble::se_pop;
  boost::python::list (Ensemble::*expt_se_pop_v4)() = &Ensemble::se_pop;


  void (Ensemble::*expt_sh_pop_v1)(vector<double>& v,double xmax,double xmin) = &Ensemble::sh_pop;
  void (Ensemble::*expt_sh_pop_v2)(vector<double>& v) = &Ensemble::sh_pop;
  boost::python::list (Ensemble::*expt_sh_pop_v3)(double xmax,double xmin) = &Ensemble::sh_pop;
  boost::python::list (Ensemble::*expt_sh_pop_v4)() = &Ensemble::sh_pop;


  void (Ensemble::*expt_sh_pop1_v1)(vector<double>& v,double xmax,double xmin) = &Ensemble::sh_pop1;
  void (Ensemble::*expt_sh_pop1_v2)(vector<double>& v) = &Ensemble::sh_pop1;
  boost::python::list (Ensemble::*expt_sh_pop1_v3)(double xmax,double xmin) = &Ensemble::sh_pop1;
  boost::python::list (Ensemble::*expt_sh_pop1_v4)() = &Ensemble::sh_pop1;







  class_<Ensemble>("Ensemble",init<>())
      .def(init<int,int,int>())
      .def("__copy__", &generic__copy__<Ensemble>)
      .def("__deepcopy__", &generic__deepcopy__<Ensemble>)

      .def_readwrite("ntraj", &Ensemble::ntraj)
      .def_readwrite("nnucl", &Ensemble::nnucl)
      .def_readwrite("nelec", &Ensemble::nelec)
      .def_readwrite("is_active", &Ensemble::is_active)

      .def_readwrite("mol", &Ensemble::mol)
      .def_readwrite("el", &Ensemble::el)
//      .def_readwrite("ham", &Ensemble::ham)

//      .def("get_mol", &Ensemble::get_mol)
//      .def("get_el", &Ensemble::get_el)
//      .def("get_ham", &Ensemble::get_ham, return_value_policy<manage_new_object>())
//      .def("get_ham", &Ensemble::get_ham, return_value_policy<reference_existing_object>())
//      .def("set_mol", &Ensemble::set_mol)
//      .def("set_el", &Ensemble::set_el)
      .def("ham_set_ham", expt_ham_set_ham_v1)
      .def("ham_set_ham", expt_ham_set_ham_v2)
      .def("ham_set_ham", expt_ham_set_ham_v3)

      .def("ham_set_rep", expt_ham_set_rep_v1)
      .def("ham_set_rep", expt_ham_set_rep_v2)

      .def("ham_set_params", expt_ham_set_params_v1)
      .def("ham_set_params", expt_ham_set_params_v2)
      .def("ham_set_params", expt_ham_set_params_v3)
      .def("ham_set_params", expt_ham_set_params_v4)

      .def("ham_set_q", expt_ham_set_q_v1)
      .def("ham_set_q", expt_ham_set_q_v2)

      .def("ham_set_v", expt_ham_set_v_v1)
      .def("ham_set_v", expt_ham_set_v_v2)
      .def("ham_set_v", expt_ham_set_v_v3)

      .def("ham_compute", expt_ham_compute_v1)
      .def("ham_compute_diabatic", expt_ham_compute_diabatic_v1)
      .def("ham_compute_adiabatic", expt_ham_compute_adiabatic_v1)

      .def("ham_H", expt_ham_H_v1)
      .def("ham_dHdq", expt_ham_dHdq_v1)
      .def("ham_D", expt_ham_D_v1)
      .def("ham_nac", expt_ham_nac_v1)
      .def("ham_Hvib", expt_ham_Hvib_v1)


      .def("el_propagate_electronic", expt_el_propagate_electronic_v1)
      .def("el_propagate_electronic", expt_el_propagate_electronic_v2)

      .def("mol_propagate_q", expt_mol_propagate_q_v1)
      .def("mol_propagate_q", expt_mol_propagate_q_v2)
      .def("mol_propagate_p", expt_mol_propagate_p_v1)
      .def("mol_propagate_p", expt_mol_propagate_p_v2)


      .def("se_pop", expt_se_pop_v1)
      .def("se_pop", expt_se_pop_v2)
      .def("se_pop", expt_se_pop_v3)
      .def("se_pop", expt_se_pop_v4)
      .def("sh_pop", expt_sh_pop_v1)
      .def("sh_pop", expt_sh_pop_v2)
      .def("sh_pop", expt_sh_pop_v3)
      .def("sh_pop", expt_sh_pop_v4)

      .def("sh_pop1", expt_sh_pop1_v1)
      .def("sh_pop1", expt_sh_pop1_v2)
      .def("sh_pop1", expt_sh_pop1_v3)
      .def("sh_pop1", expt_sh_pop1_v4)


 
  ;






}


#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygensemble){
#else
BOOST_PYTHON_MODULE(libensemble){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

//  export_Mathematics_objects();
  export_Ensemble_objects();

}


}// namespace libensemble
}// namespace libdyn
}// liblibra

