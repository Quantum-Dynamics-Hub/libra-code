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
  \file libelectronic.cpp
  \brief The file implements Python export function
    
*/

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "libelectronic.h"

/// liblibra namespace
namespace liblibra{


using namespace boost::python;

/// libdyn namespace
namespace libdyn{

/// libelectronic namespace
namespace libelectronic{


void export_Electronic_objects(){
/** 
  \brief Exporter of libelectronic classes and functions

*/


  void (Electronic::*expt_propagate_electronic_v_cl1)(double,Hamiltonian&) = &Electronic::propagate_electronic;
  void (Electronic::*expt_propagate_electronic_v_cl2)(double,Hamiltonian&, CMATRIX& ) = &Electronic::propagate_electronic;

  void (Electronic::*expt_project_out_v1)(int, int) = &Electronic::project_out;
  void (Electronic::*expt_project_out_v2)(int) = &Electronic::project_out;

  void (Electronic::*expt_collapse_v1)(int, int) = &Electronic::collapse;
  void (Electronic::*expt_collapse_v2)(int) = &Electronic::collapse;



  class_<Electronic>("Electronic",init<>())
      .def(init<int>())
      .def(init<int,int>())
      .def(init<int,int,double>())
      .def(init<const Electronic&>())
      .def("__copy__", &generic__copy__<Electronic>)
      .def("__deepcopy__", &generic__deepcopy__<Electronic>)

      .def("c", &Electronic::c, "electronic amplitudes")
      .def("rho", &Electronic::rho, "density matrix")
      .def("C", &Electronic::C, "electronic amplitudes")
      .def("RHO", &Electronic::RHO, "density matrix")


      .def_readwrite("nstates", &Electronic::nstates)
      .def_readwrite("istate", &Electronic::istate)
      .def_readwrite("q", &Electronic::q)
      .def_readwrite("p", &Electronic::p)

      .def("project_out", expt_project_out_v1)
      .def("project_out", expt_project_out_v2)
      .def("collapse", expt_collapse_v1)
      .def("collapse", expt_collapse_v2)
      .def("propagate_electronic", expt_propagate_electronic_v_cl1)
      .def("propagate_electronic", expt_propagate_electronic_v_cl2)
 
  ;


  class_< ElectronicList >("ElectronicList")
      .def(vector_indexing_suite< ElectronicList >())
  ;


  
  void (*expt_propagate_electronic_v1)(double dt,Electronic& el, CMATRIX& Hvib) = &propagate_electronic;
  void (*expt_propagate_electronic_v2)(double dt,CMATRIX& Coeff, CMATRIX& Hvib) = &propagate_electronic;
  void (*expt_propagate_electronic_v3)(double dt,Electronic& el, CMATRIX& Hvib, MATRIX& S) = &propagate_electronic;
  void (*expt_propagate_electronic_v4)(double dt,Electronic& el, CMATRIX& Hvib, CMATRIX& S) = &propagate_electronic;
  void (*expt_propagate_electronic_v5)(double dt,CMATRIX& Coeff, CMATRIX& Hvib, CMATRIX& S) = &propagate_electronic;
  void (*expt_propagate_electronic_v6)(double dt, CMATRIX& C, nHamiltonian& ham, int rep) = &propagate_electronic;
  void (*expt_propagate_electronic_v7)(double dt, nHamiltonian& ham, int rep) = &propagate_electronic;

  def("propagate_electronic", expt_propagate_electronic_v1);
  def("propagate_electronic", expt_propagate_electronic_v2);
  def("propagate_electronic", expt_propagate_electronic_v3);
  def("propagate_electronic", expt_propagate_electronic_v4);
  def("propagate_electronic", expt_propagate_electronic_v5);
  def("propagate_electronic", expt_propagate_electronic_v6);
  def("propagate_electronic", expt_propagate_electronic_v7);

  void (*expt_propagate_electronic_nonHermitian_v1)(double dt,CMATRIX& Coeff, CMATRIX& Hvib) = &propagate_electronic_nonHermitian;
  def("propagate_electronic_nonHermitian", expt_propagate_electronic_nonHermitian_v1);


  void (*expt_grid_propagator_v1)(double dt, CMATRIX& Hvib, CMATRIX& S, CMATRIX& U) = &grid_propagator;
  def("grid_propagator", expt_grid_propagator_v1);


}


#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygelectronic){
#else
BOOST_PYTHON_MODULE(libelectronic){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

//  export_Mathematics_objects();
  export_Electronic_objects();

}


}// namespace libdyn
}// namespace libelectronic
}// liblibra
