/*********************************************************************************
* Copyright (C) 2015 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file libhamiltonian_mm.cpp
  \brief The file implements Python export function
    
*/

#include <memory> // for std::auto_ptr<>
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "libhamiltonian_mm.h"

using namespace boost::python;

/// libhamiltonian namespace
namespace libhamiltonian{

/// libhamiltonian_atomistic namespace
namespace libhamiltonian_atomistic{

/// libhamiltonian_mm namespace
namespace libhamiltonian_mm{

using namespace libforcefield;


void export_Hamiltonian_MM_objects(){
/** 
  \brief Exporter of the libhamiltonian_mm classes and functions

*/


  export_forcefield_objects();

  // Also export self classes:
  bool (listHamiltonian_MM::*expt_is_active_v1)(Atom&,Atom&) = &listHamiltonian_MM::is_active;
  bool (listHamiltonian_MM::*expt_is_active_v2)(Atom&,Atom&,Atom&) = &listHamiltonian_MM::is_active;
  bool (listHamiltonian_MM::*expt_is_active_v3)(Atom&,Atom&,Atom&,Atom&) = &listHamiltonian_MM::is_active;



  class_<Hamiltonian_MM>("Hamiltonian_MM",init<>())
      .def("__copy__", &generic__copy__<Hamiltonian_MM>)
      .def("__deepcopy__", &generic__deepcopy__<Hamiltonian_MM>)

      .def("activate", &Hamiltonian_MM::activate)
      .def("deactivate", &Hamiltonian_MM::deactivate)
      .def("is_origin", &Hamiltonian_MM::is_origin)
      .def("set_respa_type", &Hamiltonian_MM::set_respa_type)
      .def("get_respa_type", &Hamiltonian_MM::get_respa_type)
      .def("set_interaction_type_and_functional", &Hamiltonian_MM::set_interaction_type_and_functional)
      .def("activate", &Hamiltonian_MM::activate)

  ;

  class_<listHamiltonian_MM>("listHamiltonian_MM",init<>())
      .def("__copy__", &generic__copy__<listHamiltonian_MM>)
      .def("__deepcopy__", &generic__deepcopy__<listHamiltonian_MM>)

      .def("is_new_interaction", &listHamiltonian_MM::is_new_interaction)
      .def("show_interactions_statistics", &listHamiltonian_MM::show_interactions_statistics)

      .def("set_atom_types", &listHamiltonian_MM::set_atom_types)
      .def("set_fragment_types", &listHamiltonian_MM::set_fragment_types)

      .def("set_atom_interactions_for_atoms", &listHamiltonian_MM::set_atom_interactions_for_atoms)
      .def("set_group_interactions_for_atoms", &listHamiltonian_MM::set_group_interactions_for_atoms)

      .def("set_interactions_for_atoms", &listHamiltonian_MM::set_interactions_for_atoms)
      .def("set_interactions_for_fragments", &listHamiltonian_MM::set_interactions_for_fragments)

      .def("apply_pbc_to_interactions", &listHamiltonian_MM::apply_pbc_to_interactions)
      .def("set_respa_types", &listHamiltonian_MM::set_respa_types)

      .def("is_active", expt_is_active_v1)
      .def("is_active", expt_is_active_v2)
      .def("is_active", expt_is_active_v3)


  ;




}


#ifdef CYGWIN
BOOST_PYTHON_MODULE(cyghamiltonian_mm){
#else
BOOST_PYTHON_MODULE(libhamiltonian_mm){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

//  export_Mathematics_objects();
  export_Hamiltonian_MM_objects();

}




}// namespace libhamiltonian_mm
}// namespace libhamiltonian_atomistic
}// namespace libhamiltonian


