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
  \file libhamiltonian_atomistic.cpp
  \brief The file implements Python export function
    
*/

#if defined(USING_PCH)
#include "../pch.h"
#else
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#endif

#include "libatomistic.h"


/// liblibra namespace
namespace liblibra{

using namespace boost::python;


/// libhamiltonian namespace
namespace libatomistic{

using namespace libhamiltonian_mm;
using namespace libhamiltonian_qm;
using namespace libhamiltonian_generic;



void export_atomistic_objects(){
/** 
  \brief Exporter of the libhamiltonian_atomistic classes and functions

*/

  export_hamiltonian_generic_objects();
  export_hamiltonian_mm_objects();
  export_hamiltonian_qm_objects();



//  void (Hamiltonian_Atomistic::*set_params)(boost::python::list) = &Hamiltonian_Atomistic::set_params;
  void (Hamiltonian_Atomistic::*expt_set_q_v1)(vector<double>&) = &Hamiltonian_Atomistic::set_q;
  void (Hamiltonian_Atomistic::*expt_set_q_v2)(boost::python::list) = &Hamiltonian_Atomistic::set_q;
  void (Hamiltonian_Atomistic::*expt_set_v_v1)(vector<double>&) = &Hamiltonian_Atomistic::set_v;
  void (Hamiltonian_Atomistic::*expt_set_v_v2)(boost::python::list) = &Hamiltonian_Atomistic::set_v;

  bool (Hamiltonian_Atomistic::*expt_is_active_v1a)(Atom&,Atom&) = &Hamiltonian_Atomistic::is_active;
  bool (Hamiltonian_Atomistic::*expt_is_active_v2a)(Atom&,Atom&,Atom&) = &Hamiltonian_Atomistic::is_active;
  bool (Hamiltonian_Atomistic::*expt_is_active_v3a)(Atom&,Atom&,Atom&,Atom&) = &Hamiltonian_Atomistic::is_active;


  class_<Hamiltonian_Atomistic, bases<Hamiltonian> >("Hamiltonian_Atomistic",init<int,int>())
      .def("__copy__", &generic__copy__<Hamiltonian_Atomistic>)
      .def("__deepcopy__", &generic__deepcopy__<Hamiltonian_Atomistic>)

/*
      .def("set_params", set_params)
      .def("set_rep", &Hamiltonian_Atomistic::set_rep)

      .def("compute",          &Hamiltonian_Atomistic::compute)
*/
      .def("set_q", expt_set_q_v1)
      .def("set_q", expt_set_q_v2)
      .def("set_v", expt_set_v_v1)
      .def("set_v", expt_set_v_v2)


      .def("set_Hamiltonian_type", &Hamiltonian_Atomistic::set_Hamiltonian_type)


//      .def("is_new_interaction", &Hamiltonian_Atomistic::is_new_interaction)
      .def("show_interactions_statistics", &Hamiltonian_Atomistic::show_interactions_statistics)
      .def("show_interactions", &Hamiltonian_Atomistic::show_interactions)

      .def("set_atom_types", &Hamiltonian_Atomistic::set_atom_types)
      .def("set_fragment_types", &Hamiltonian_Atomistic::set_fragment_types)

      .def("set_atom_interactions_for_atoms", &Hamiltonian_Atomistic::set_atom_interactions_for_atoms)
      .def("set_group_interactions_for_atoms", &Hamiltonian_Atomistic::set_group_interactions_for_atoms)

      .def("set_interactions_for_atoms", &Hamiltonian_Atomistic::set_interactions_for_atoms)
      .def("set_interactions_for_fragments", &Hamiltonian_Atomistic::set_interactions_for_fragments)

      .def("apply_pbc_to_interactions", &Hamiltonian_Atomistic::apply_pbc_to_interactions)
      .def("set_respa_types", &Hamiltonian_Atomistic::set_respa_types)

      .def("init_qm_Hamiltonian",&Hamiltonian_Atomistic::init_qm_Hamiltonian)
      .def("add_excitation",&Hamiltonian_Atomistic::add_excitation)
      .def("excite_alp",&Hamiltonian_Atomistic::excite_alp)
      .def("excite_bet",&Hamiltonian_Atomistic::excite_bet)

      .def("is_active", expt_is_active_v1a)
      .def("is_active", expt_is_active_v2a)
      .def("is_active", expt_is_active_v3a)


      .def("set_system", &Hamiltonian_Atomistic::set_system)
      .def("compute_diabatic", &Hamiltonian_Atomistic::compute_diabatic)
      .def("compute_adiabatic",&Hamiltonian_Atomistic::compute_adiabatic)

      .def("get_stress", &Hamiltonian_Atomistic::get_stress)

/*
      .def("H", &Hamiltonian_Atomistic::H)
      .def("dHdq", &Hamiltonian_Atomistic::dHdq)
      .def("Hvib", &Hamiltonian_Atomistic::Hvib)
      .def("D", &Hamiltonian_Atomistic::D)
      .def("nac", &Hamiltonian_Atomistic::nac)
*/
 
  ;



}


#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygatomistic){
#else
BOOST_PYTHON_MODULE(libatomistic){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

  export_atomistic_objects();

}


}// namespace libatomistic
}// liblibra

