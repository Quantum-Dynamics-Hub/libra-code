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
  \file libhamiltonian_generic.cpp
  \brief The file implements Python export function
    
*/

#if defined(USING_PCH)
#include "../../pch.h"
#else
#include <memory> // for std::auto_ptr<>
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#endif 

#include "libhamiltonian_generic.h"

/// liblibra namespace
namespace liblibra{


using namespace boost::python;

/// libhamiltonian namespace
namespace libatomistic{

/// libhamiltonian_generic namespace
namespace libhamiltonian_generic{


void export_hamiltonian_generic_objects(){
/** 
  \brief Exporter of the libhamiltonian_generic classes and functions

*/




  void (Hamiltonian::*set_params)(boost::python::list) = &Hamiltonian::set_params;
  void (Hamiltonian::*set_q)(boost::python::list) = &Hamiltonian::set_q;
  void (Hamiltonian::*set_v)(boost::python::list) = &Hamiltonian::set_v;


  class_<Hamiltonian>("Hamiltonian",init<>())
      .def("__copy__", &generic__copy__<Hamiltonian>)
      .def("__deepcopy__", &generic__deepcopy__<Hamiltonian>)

      .def("set_rep", &Hamiltonian::set_rep)
      .def("set_params", set_params)
      .def("set_q", set_q)
      .def("set_v", set_v)

      .def("get_ham_dia", &Hamiltonian::get_ham_dia)
      .def("get_ham_adi", &Hamiltonian::get_ham_adi)
      .def("get_d1ham_dia", &Hamiltonian::get_d1ham_dia)
      .def("get_d1ham_adi", &Hamiltonian::get_d1ham_adi)

      .def("compute",          &Hamiltonian::compute)
//      .def("compute_diabatic", &Hamiltonian_Model::compute_diabatic)
//      .def("compute_adiabatic",&Hamiltonian_Model::compute_adiabatic)

      .def_readwrite("status_adi", &Hamiltonian::status_adi)
      .def_readwrite("status_dia", &Hamiltonian::status_dia)

      .def("H", &Hamiltonian::H)
      .def("dHdq", &Hamiltonian::dHdq)
      .def("Hvib", &Hamiltonian::Hvib)
      .def("D", &Hamiltonian::D)
      .def("nac", &Hamiltonian::nac)
  ;

  class_< HamiltonianList >("HamiltonianList")
      .def(vector_indexing_suite< HamiltonianList >())
  ;


}


#ifdef CYGWIN
BOOST_PYTHON_MODULE(cyghamiltonian_generic){
#else
BOOST_PYTHON_MODULE(libhamiltonian_generic){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

//  export_Mathematics_objects();
  export_hamiltonian_generic_objects();

}



}// namespace libhamiltonian_generic
}// namespace libatomistic
}// liblibra

