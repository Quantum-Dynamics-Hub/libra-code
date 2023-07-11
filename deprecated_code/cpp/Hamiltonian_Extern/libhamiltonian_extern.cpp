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
  \file libhamiltonian_extern.cpp
  \brief The file implements Python export function
    
*/

#if defined(USING_PCH)
#include "../../pch.h"
#else
#include <memory> // for std::auto_ptr<>
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#endif 

#include "libhamiltonian_extern.h"

/// liblibra namespace
namespace liblibra{


using namespace boost::python;


/// libhamiltonian namespace
namespace libhamiltonian{

/// libhamiltonian_extern namespace
namespace libhamiltonian_extern{


void export_hamiltonian_extern_objects(){
/** 
  \brief Exporter of the libhamiltonian_extern classes and functions

*/


  class_<Hamiltonian_Extern, bases<Hamiltonian> >("Hamiltonian_Extern",init<int,int>())
      .def("__copy__", &generic__copy__<Hamiltonian_Extern>)
      .def("__deepcopy__", &generic__deepcopy__<Hamiltonian_Extern>)
      .def("set_adiabatic_opt", &Hamiltonian_Extern::set_adiabatic_opt)
      .def("set_vibronic_opt", &Hamiltonian_Extern::set_vibronic_opt)
      .def("bind_ham_dia", &Hamiltonian_Extern::bind_ham_dia)
      .def("bind_d1ham_dia", &Hamiltonian_Extern::bind_d1ham_dia)
      .def("bind_d2ham_dia", &Hamiltonian_Extern::bind_d2ham_dia)
      .def("bind_ham_adi", &Hamiltonian_Extern::bind_ham_adi)
      .def("bind_d1ham_adi", &Hamiltonian_Extern::bind_d1ham_adi)
      .def("bind_ham_vib", &Hamiltonian_Extern::bind_ham_vib)
      .def("compute_diabatic", &Hamiltonian_Extern::compute_diabatic)
      .def("compute_adiabatic",&Hamiltonian_Extern::compute_adiabatic)
      .def("Hvib",&Hamiltonian_Extern::Hvib)

  ;

}


#ifdef CYGWIN
BOOST_PYTHON_MODULE(cyghamiltonian_extern){
#else
BOOST_PYTHON_MODULE(libhamiltonian_extern){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

  export_hamiltonian_extern_objects();

}


}// namespace libhamiltonian_extern
}// namespace libhamiltonian
}// liblibra

