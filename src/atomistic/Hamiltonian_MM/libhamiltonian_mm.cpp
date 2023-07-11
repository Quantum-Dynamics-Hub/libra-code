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
  \file libhamiltonian_mm.cpp
  \brief The file implements Python export function
    
*/

#if defined(USING_PCH)
#include "../../pch.h"
#else
#include <memory> // for std::auto_ptr<>
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#endif 

#include "libhamiltonian_mm.h"


/// liblibra namespace
namespace liblibra{


using namespace boost::python;

namespace libatomistic{

/// libhamiltonian_mm namespace
namespace libhamiltonian_mm{

using namespace libforcefield;

void export_hamiltonian_mm_objects(){
/** 
  \brief Exporter of the libhamiltonian_mm classes and functions
*/

  void (Interaction_N_Body::*expt_set_coords_v1)(VECTOR& r_, int indx) = &Interaction_N_Body::set_coords;
  void (Interaction_N_Body::*expt_set_coords_v2)(vector<VECTOR>& r_, vector<int>& indxs_) = &Interaction_N_Body::set_coords;

  void (Interaction_N_Body::*expt_set_transl_v1)(VECTOR& r_, int indx) = &Interaction_N_Body::set_transl;
  void (Interaction_N_Body::*expt_set_transl_v2)(vector<VECTOR>& t_, vector<int>& indxs_) = &Interaction_N_Body::set_transl;

  void (Interaction_N_Body::*expt_set_forces_v1)(VECTOR& r_, int indx) = &Interaction_N_Body::set_forces;
  void (Interaction_N_Body::*expt_set_forces_v2)(vector<VECTOR>& r_, vector<int>& indxs_) = &Interaction_N_Body::set_forces;

  void (Interaction_N_Body::*expt_set_charges_v1)(double& q_, int indx) = &Interaction_N_Body::set_charges;
  void (Interaction_N_Body::*expt_set_charges_v2)(vector<double>& q_, vector<int>& indxs_) = &Interaction_N_Body::set_charges;

  void (Interaction_N_Body::*expt_set_hessian_v1)(MATRIX& hess_, vector<int>& hess_stenc_) = &Interaction_N_Body::set_hessian;
  void (Interaction_N_Body::*expt_set_functional_v1)(std::string f) = &Interaction_N_Body::set_functional;



  class_<Interaction_N_Body>("Interaction_N_Body",init<int>())
      .def_readwrite("energy", &Interaction_N_Body::energy)
      .def_readwrite("stress_at", &Interaction_N_Body::stress_at)

      .def_readwrite("Nbody", &Interaction_N_Body::Nbody)
      .def_readwrite("is_active", &Interaction_N_Body::is_active)
      .def_readwrite("int_type", &Interaction_N_Body::int_type)
      .def_readwrite("functional", &Interaction_N_Body::functional)
      .def("set_coords", expt_set_coords_v1)
      .def("set_coords", expt_set_coords_v2)
      .def("set_transl", expt_set_transl_v1)
      .def("set_transl", expt_set_transl_v2)
      .def("set_forces", expt_set_forces_v1)
      .def("set_forces", expt_set_forces_v2)
      .def("set_charges", expt_set_charges_v1)
      .def("set_charges", expt_set_charges_v2)
      .def("set_hessian", expt_set_hessian_v1)
      .def("set_functional", expt_set_functional_v1)

  ;



  void (Interaction_2_Body::*expt_set_coords_v2_2b)(VECTOR& r1_, VECTOR& r2_) = &Interaction_2_Body::set_coords;
  void (Interaction_2_Body::*expt_set_transl_v2_2b)(VECTOR& t1_, VECTOR& t2_) = &Interaction_2_Body::set_transl;
  void (Interaction_2_Body::*expt_set_forces_v2_2b)(VECTOR& f1_, VECTOR& f2_) = &Interaction_2_Body::set_forces;
  void (Interaction_2_Body::*expt_set_charges_v2_2b)(double& q1_, double& q2_) = &Interaction_2_Body::set_charges;

  class_<Interaction_2_Body, bases<Interaction_N_Body> >("Interaction_2_Body",init<>())
      .def_readwrite("energy", &Interaction_N_Body::energy)
      .def_readwrite("stress_at", &Interaction_N_Body::stress_at)

      .def("set_coords", expt_set_coords_v1)
      .def("set_coords", expt_set_coords_v2)
      .def("set_transl", expt_set_transl_v1)
      .def("set_transl", expt_set_transl_v2)
      .def("set_forces", expt_set_forces_v1)
      .def("set_forces", expt_set_forces_v2)
      .def("set_charges", expt_set_charges_v1)
      .def("set_charges", expt_set_charges_v2)
      .def("set_hessian", expt_set_hessian_v1)
      .def("set_functional", expt_set_functional_v1)

      .def("set_coords", expt_set_coords_v2_2b)
      .def("set_transl", expt_set_transl_v2_2b)
      .def("set_forces", expt_set_forces_v2_2b)
      .def("set_charges", expt_set_charges_v2_2b)

  ;

  void (Bond_Interaction::*expt_set_params_v1_2b)(boost::python::dict) = &Bond_Interaction::set_params;
  void (VdW_Interaction::*expt_set_params_v1_2v)(boost::python::dict) = &VdW_Interaction::set_params;
  void (Elec_Interaction::*expt_set_params_v1_2e)(boost::python::dict) = &Elec_Interaction::set_params;

  void (Bond_Interaction::*expt_compute_v1_2b)() = &Bond_Interaction::compute;
  void (VdW_Interaction::*expt_compute_v1_2v)(double Ron, double Roff) = &VdW_Interaction::compute;
  void (VdW_Interaction::*expt_compute_v2_2v)() = &VdW_Interaction::compute;
  void (Elec_Interaction::*expt_compute_v1_2e)(double Ron, double Roff) = &Elec_Interaction::compute;
  void (Elec_Interaction::*expt_compute_v2_2e)() = &Elec_Interaction::compute;



  class_<Bond_Interaction, bases<Interaction_2_Body> >("Bond_Interaction",init<>())
      .def_readwrite("energy", &Interaction_N_Body::energy)
      .def_readwrite("stress_at", &Interaction_N_Body::stress_at)

      .def_readwrite("K", &Bond_Interaction::K)
      .def_readwrite("D", &Bond_Interaction::D)
      .def_readwrite("r0", &Bond_Interaction::r0)
      .def_readwrite("alpha", &Bond_Interaction::alpha)
      .def("set_functional", &Bond_Interaction::set_functional)
      .def("set_params", expt_set_params_v1_2b)
      .def("compute", expt_compute_v1_2b)
  ;

  class_<VdW_Interaction, bases<Interaction_2_Body> >("VdW_Interaction",init<>())
      .def_readwrite("energy", &Interaction_N_Body::energy)
      .def_readwrite("stress_at", &Interaction_N_Body::stress_at)

      .def_readwrite("sigma", &VdW_Interaction::sigma)
      .def_readwrite("epsilon", &VdW_Interaction::epsilon)
      .def_readwrite("D", &VdW_Interaction::D)
      .def_readwrite("r0", &VdW_Interaction::r0)
      .def_readwrite("alpha", &VdW_Interaction::alpha)
      .def_readwrite("scale", &VdW_Interaction::scale)
      .def("set_functional", &VdW_Interaction::set_functional)
      .def("set_params", expt_set_params_v1_2v)
      .def("compute", expt_compute_v1_2v)
      .def("compute", expt_compute_v2_2v)
  ;

  class_<Elec_Interaction, bases<Interaction_2_Body> >("Elec_Interaction",init<>())
      .def_readwrite("energy", &Interaction_N_Body::energy)
      .def_readwrite("stress_at", &Interaction_N_Body::stress_at)

      .def_readwrite("J", &Elec_Interaction::J)
      .def_readwrite("xi1", &Elec_Interaction::xi1)
      .def_readwrite("xi2", &Elec_Interaction::xi2)
      .def_readwrite("eps", &Elec_Interaction::eps)
      .def_readwrite("delta", &Elec_Interaction::delta)
      .def_readwrite("scale", &Elec_Interaction::scale)
      .def("set_functional", &Elec_Interaction::set_functional)
      .def("set_params", expt_set_params_v1_2e)
      .def("compute", expt_compute_v1_2e)
      .def("compute", expt_compute_v2_2e)

  ;




  void (Interaction_3_Body::*expt_set_coords_v1_3b)(VECTOR& r1_, VECTOR& r2_, VECTOR& r3_) = &Interaction_3_Body::set_coords;
  void (Interaction_3_Body::*expt_set_transl_v1_3b)(VECTOR& t1_, VECTOR& t2_, VECTOR& t3_) = &Interaction_3_Body::set_transl;
  void (Interaction_3_Body::*expt_set_forces_v1_3b)(VECTOR& f1_, VECTOR& f2_, VECTOR& f3_) = &Interaction_3_Body::set_forces;
  void (Interaction_3_Body::*expt_set_charges_v1_3b)(double& q1_, double& q2_, double& q3_) = &Interaction_3_Body::set_charges;

  void (Angle_Interaction::*expt_compute_v1_3b)() = &Angle_Interaction::compute;

  class_<Interaction_3_Body, bases<Interaction_N_Body> >("Interaction_3_Body",init<>())
      .def_readwrite("energy", &Interaction_N_Body::energy)
      .def_readwrite("stress_at", &Interaction_N_Body::stress_at)

      .def("set_coords", expt_set_coords_v1)
      .def("set_coords", expt_set_coords_v2)
      .def("set_transl", expt_set_transl_v1)
      .def("set_transl", expt_set_transl_v2)
      .def("set_forces", expt_set_forces_v1)
      .def("set_forces", expt_set_forces_v2)
      .def("set_charges", expt_set_charges_v1)
      .def("set_charges", expt_set_charges_v2)
      .def("set_hessian", expt_set_hessian_v1)
      .def("set_functional", expt_set_functional_v1)

      .def("set_coords", expt_set_coords_v1_3b)
      .def("set_transl", expt_set_transl_v1_3b)
      .def("set_forces", expt_set_forces_v1_3b)
      .def("set_charges", expt_set_charges_v1_3b)
  ;

  class_<Angle_Interaction, bases<Interaction_3_Body> >("Angle_Interaction",init<>())
      .def_readwrite("energy", &Interaction_N_Body::energy)
      .def_readwrite("stress_at", &Interaction_N_Body::stress_at)

      .def_readwrite("k_theta", &Angle_Interaction::k_theta)
      .def_readwrite("theta_0", &Angle_Interaction::theta_0)
      .def_readwrite("cos_theta_0", &Angle_Interaction::cos_theta_0)
      .def_readwrite("C0", &Angle_Interaction::C0)
      .def_readwrite("C1", &Angle_Interaction::C1)
      .def_readwrite("C2", &Angle_Interaction::C2)
      .def_readwrite("coordination", &Angle_Interaction::coordination)

      .def("set_functional", &Angle_Interaction::set_functional)
      .def("set_params", &Angle_Interaction::set_params)

      .def("compute", expt_compute_v1_3b)
  ;





  void (Interaction_4_Body::*expt_set_coords_v1_4b)(VECTOR& r1_, VECTOR& r2_, VECTOR& r3_, VECTOR& r4_) = &Interaction_4_Body::set_coords;
  void (Interaction_4_Body::*expt_set_transl_v1_4b)(VECTOR& t1_, VECTOR& t2_, VECTOR& t3_, VECTOR& t4_) = &Interaction_4_Body::set_transl;
  void (Interaction_4_Body::*expt_set_forces_v1_4b)(VECTOR& f1_, VECTOR& f2_, VECTOR& f3_, VECTOR& f4_) = &Interaction_4_Body::set_forces;
  void (Interaction_4_Body::*expt_set_charges_v1_4b)(double& q1_, double& q2_, double& q3_, double& q4_) = &Interaction_4_Body::set_charges;

  void (Dihedral_Interaction::*expt_compute_v1_4d)() = &Dihedral_Interaction::compute;
  void (OOP_Interaction::*expt_compute_v1_4o)() = &OOP_Interaction::compute;

  class_<Interaction_4_Body, bases<Interaction_N_Body> >("Interaction_4_Body",init<>())
      .def_readwrite("energy", &Interaction_N_Body::energy)
      .def_readwrite("stress_at", &Interaction_N_Body::stress_at)

      .def("set_coords", expt_set_coords_v1)
      .def("set_coords", expt_set_coords_v2)
      .def("set_transl", expt_set_transl_v1)
      .def("set_transl", expt_set_transl_v2)
      .def("set_forces", expt_set_forces_v1)
      .def("set_forces", expt_set_forces_v2)
      .def("set_charges", expt_set_charges_v1)
      .def("set_charges", expt_set_charges_v2)
      .def("set_hessian", expt_set_hessian_v1)
      .def("set_functional", expt_set_functional_v1)

      .def("set_coords", expt_set_coords_v1_4b)
      .def("set_transl", expt_set_transl_v1_4b)
      .def("set_forces", expt_set_forces_v1_4b)
      .def("set_charges", expt_set_charges_v1_4b)
  ;

  class_<Dihedral_Interaction, bases<Interaction_4_Body> >("Dihedral_Interaction",init<>())
      .def_readwrite("energy", &Interaction_N_Body::energy)
      .def_readwrite("stress_at", &Interaction_N_Body::stress_at)

      .def_readwrite("Vphi", &Dihedral_Interaction::Vphi)
      .def_readwrite("phi0", &Dihedral_Interaction::phi0)
      .def_readwrite("Vphi1", &Dihedral_Interaction::Vphi1)
      .def_readwrite("Vphi2", &Dihedral_Interaction::Vphi2)
      .def_readwrite("Vphi3", &Dihedral_Interaction::Vphi3)
      .def_readwrite("n", &Dihedral_Interaction::n)
      .def_readwrite("opt", &Dihedral_Interaction::opt)

      .def("set_functional", &Dihedral_Interaction::set_functional)
      .def("set_params", &Dihedral_Interaction::set_params)

      .def("compute", expt_compute_v1_4d)
  ;



  class_<OOP_Interaction, bases<Interaction_4_Body> >("OOP_Interaction",init<>())
      .def_readwrite("energy", &Interaction_N_Body::energy)
      .def_readwrite("stress_at", &Interaction_N_Body::stress_at)

      .def_readwrite("K", &OOP_Interaction::K)
      .def_readwrite("C0", &OOP_Interaction::C0)
      .def_readwrite("C1", &OOP_Interaction::C1)
      .def_readwrite("C2", &OOP_Interaction::C2)
      .def_readwrite("xi_0", &OOP_Interaction::xi_0)
      .def_readwrite("opt", &OOP_Interaction::opt)

      .def("set_functional", &OOP_Interaction::set_functional)
      .def("set_params", &OOP_Interaction::set_params)

      .def("compute", expt_compute_v1_4o)
  ;




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
      .def("show_interactions", &listHamiltonian_MM::show_interactions)

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

  export_hamiltonian_mm_objects();

}




}// namespace libhamiltonian_mm
}// namespace libatomistic
}// liblibra

