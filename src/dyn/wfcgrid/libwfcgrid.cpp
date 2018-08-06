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
  \file libwfcgrid.cpp
  \brief The file implements Python export function
    
*/

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "libwfcgrid.h"

/// liblibra namespace
namespace liblibra{

using namespace boost::python;

/// libdyn namespace
namespace libdyn{

/// libwfcgrid namespace
namespace libwfcgrid{


void export_Wfcgrid_objects(){
/** 
  \brief Exporter of libwfcgrid classes and functions

*/


  boost::python::list (Wfcgrid::*expt_absorb_1D)(double dL) = &Wfcgrid::absorb_1D;

  void (Wfcgrid::*expt_update_potential_1D_v1)(Hamiltonian& ham) = &Wfcgrid::update_potential_1D;
  void (Wfcgrid::*expt_update_potential_1D_v2)(bp::object py_funct, bp::object params) = &Wfcgrid::update_potential_1D;
  void (Wfcgrid::*expt_update_potential_2D_v1)(Hamiltonian& ham) = &Wfcgrid::update_potential_2D;
  void (Wfcgrid::*expt_update_potential_2D_v2)(bp::object py_funct, bp::object params) = &Wfcgrid::update_potential_2D;




  class_<Wfcgrid>("Wfcgrid",init<double,double,double,int>())
      .def(init<double, double, double, double, double, double, int>())
      .def(init<const Wfcgrid&>())
      .def("__copy__", &generic__copy__<Wfcgrid>)
      .def("__deepcopy__", &generic__deepcopy__<Wfcgrid>)

      .def("init_wfc_1D", &Wfcgrid::init_wfc_1D)
      .def("init_wfc_2D", &Wfcgrid::init_wfc_2D)

      .def("print_wfc_1D", &Wfcgrid::print_wfc_1D)
      .def("print_wfc_2D", &Wfcgrid::print_wfc_2D)
      .def("print_reci_wfc_1D", &Wfcgrid::print_reci_wfc_1D)
      .def("print_ham_1D", &Wfcgrid::print_ham_1D)
      .def("print_expH_1D", &Wfcgrid::print_expH_1D)
      .def("print_expK_1D", &Wfcgrid::print_expK_1D)

      .def("print_populations_1D", &Wfcgrid::print_populations_1D)
      .def("print_populations_2D", &Wfcgrid::print_populations_2D)

      .def("update_potential_1D", expt_update_potential_1D_v1)
      .def("update_potential_1D", expt_update_potential_1D_v2)
      .def("update_potential_2D", expt_update_potential_2D_v1)
      .def("update_potential_2D", expt_update_potential_2D_v2)

      .def("update_propagator_1D", &Wfcgrid::update_propagator_1D)
      .def("update_propagator_2D", &Wfcgrid::update_propagator_2D)

      .def("update_propagator_K_1D", &Wfcgrid::update_propagator_K_1D)
      .def("update_propagator_K_2D", &Wfcgrid::update_propagator_K_2D)

      .def("propagate_exact_1D", &Wfcgrid::propagate_exact_1D)
      .def("propagate_exact_2D", &Wfcgrid::propagate_exact_2D)

      .def("absorb_1D",expt_absorb_1D)

      .def("e_kin_1D", &Wfcgrid::e_kin_1D)
      .def("e_pot_1D", &Wfcgrid::e_pot_1D)
      .def("e_tot_1D", &Wfcgrid::e_tot_1D)
 
  ;



}


#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygwfcgrid){
#else
BOOST_PYTHON_MODULE(libwfcgrid){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

//  export_Mathematics_objects();
  export_Wfcgrid_objects();

}


}// namespace libwfcgrid
}// namespace libdyn
}// liblibra
