/*********************************************************************************
* Copyright (C) 2021-2022 Mohammad Shakiba and Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file liblibint2_wrappers.cpp
  \brief The file implements Python export function
    
*/

#define BOOST_PYTHON_MAX_ARITY 30

#if defined(USING_PCH)
#include "../pch.h"
#else
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#endif 

#include "libint2_wrappers.h"


/// liblibra namespace
namespace liblibra{

using namespace boost::python;
using namespace liblinalg;
using namespace libspecialfunctions;



/// libqobjects namespace
namespace liblibint2_wrappers{




void export_libint2_wrappers_objects(){
/** 
  \brief Exporter of liblibint2_wrappers classes and functions

*/
/*
  double (*expt_electron_repulsion_integral_AO_v2)
  ( AO& AOa, AO& AOb, AO& AOc, AO& AOd, int is_normalize ) = &electron_repulsion_integral;
  boost::python::list (*expt_electron_repulsion_integral_AO_v3)
  ( AO& AOa, AO& AOb, AO& AOc, AO& AOd, int is_normalize, int is_derivs) = &electron_repulsion_integral;

  // ============ Now export functions =============
  // Overlaps
  def("gaussian_overlap", expt_gaussian_overlap_G_v1);


  void (PrimitiveG::*expt_set_position_PrimitiveG_v1)(VECTOR) = &PrimitiveG::set_position;
  void (AO::*expt_set_position_AO_v1)(VECTOR) = &AO::set_position;


*/

  size_t (*expt_nbasis_v1)(const std::vector<libint2::Shell>& shells) = &nbasis;
  def("nbasis", expt_nbasis_v1);


  std::vector<libint2::Shell> (*expt_initialize_shell_v1)(
  int l_val, bool is_spherical , const std::vector<double>& exponents, 
  const std::vector<double>& coeff, VECTOR& coords) = &initialize_shell;

  def("initialize_shell", expt_initialize_shell_v1);


  void (*expt_add_to_shell_v1)(std::vector<libint2::Shell>& shells, 
  int l_val, bool is_spherical , const std::vector<double>& exponents, 
  const std::vector<double>& coeff, VECTOR& coords) = &add_to_shell;

  def("add_to_shell", expt_add_to_shell_v1);


  void (*expt_print_shells_v1)(std::vector<libint2::Shell>& shells) = &print_shells;  
  def("print_shell", expt_print_shells_v1);


  class_<libint2::Shell>("libint2_Shell",init<>())
      .def(init<const libint2::Shell&>())
      .def("__copy__", &generic__copy__<libint2::Shell>) 
      .def("__deepcopy__", &generic__deepcopy__<libint2::Shell>)
      .def_readwrite("alpha",&libint2::Shell::alpha)
      .def_readwrite("contr",&libint2::Shell::contr)
      //.def("init",&PrimitiveG::init)
  ;

  class_< libint2_ShellList >("libint2_ShellList")
      .def(vector_indexing_suite< libint2_ShellList >())
  ;

  class_<Matrix>("Matrix",init<>())
      .def("__copy__", &generic__copy__<Matrix>) 
      .def("__deepcopy__", &generic__deepcopy__<Matrix>)
  ;

  MATRIX (*expt_compute_overlaps_v1)
  (const std::vector<libint2::Shell>& shells_1, 
   const std::vector<libint2::Shell>& shells_2, int number_of_threads) = &compute_overlaps; 
 
  def("compute_overlaps", expt_compute_overlaps_v1);

/*
  MATRIX (*expt_compute_overlaps_serial_v1)
  (const std::vector<libint2::Shell>& shells_1, const std::vector<libint2::Shell>& shells_2) = &compute_overlaps_serial;
  def("compute_overlaps_serial", expt_compute_overlaps_serial_v1);
*/


/**
  class_<libint2::Operator>("libint2_Operator",init<>())
      .def(init<const libint2::Operator&>())
      .def("__copy__", &generic__copy__<libint2::Operator>) 
      .def("__deepcopy__", &generic__deepcopy__<libint2::Operator>)
      //.def_readwrite("x_exp",&PrimitiveG::x_exp)
      .def("nbasis",&libint2::Operator::nbasis)
  ;

*/
/*
    py::class_<libint2::Operator>(m, "Operator");
    m.def("nbasis", &nbasis, "A function for number of basis sets", py::return_value_policy::reference_internal);
    m.def("compute_1body_ints_parallel", &compute_1body_ints_parallel, "A function for computing AO overlaps parallel", py::return_value_policy::reference_internal);
    m.def("compute_1body_ints", &compute_1body_ints, "A function for computing integrals serial", py::return_value_policy::reference_internal);
    m.def("initialize_shell", &initialize_shell, "A function for creating shells", py::return_value_policy::reference_internal);
    m.def("initialize_atom", &initialize_atom, "A function for initializing Atom", py::return_value_policy::reference_internal);
    m.def("add_to_shell", &add_to_shell, "A function to add the exponents, contraction coefficients, and Atom structure to a predefined shell", py::return_value_policy::reference_internal);
    m.def("print_shell", &print_shell, "A function for printing a shell", py::return_value_policy::reference_internal);
    m.def("compute_overlaps", &compute_overlaps, "A function that computes overlaps between pairs of shells", py::return_value_policy::reference_internal);

*/

}


#ifdef CYGWIN
BOOST_PYTHON_MODULE(cyglibint2_wrappers){
#else
BOOST_PYTHON_MODULE(liblibint2_wrappers){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

  export_libint2_wrappers_objects();

}


}// namespace libqobject
}// namespace liblibra

