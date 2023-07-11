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
  \file libchemsys.cpp
  \brief The file implements Python export function
    
*/

#if defined(USING_PCH)
#include "../../pch.h"
#else
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#endif
#include "libchemsys.h"

/// liblibra namespace
namespace liblibra{

using namespace boost::python;

/// libchemobjects namespace
namespace libchemobjects{

/// libchemsys namespace
namespace libchemsys{


void export_Chemsys_objects(){
/** 
  \brief Exporter of libchemsys classes and functions

*/


//void (System::*CREATE_ATOM1)()      = &System::CREATE_ATOM;
void (System::*CREATE_ATOM2)(Atom)  = &System::CREATE_ATOM;
void (System::*LINK_ATOMS1)(Atom&,Atom&) = &System::LINK_ATOMS;
void (System::*LINK_ATOMS2)(int,int) = &System::LINK_ATOMS;

void (System::*init_box1)()         = &System::init_box;
void (System::*init_box2)(double,double,double) = &System::init_box;
void (System::*init_box3)(VECTOR,VECTOR,VECTOR) = &System::init_box;

void (System::*print_ent1)(std::string) = &System::print_ent;
void (System::*print_ent2)(std::string,int,std::string) = &System::print_ent;
void (System::*print_ent3)(std::string,boost::python::list) = &System::print_ent;
void (System::*print_ent4)(std::string,boost::python::list,int,std::string) = &System::print_ent;

void (System::*print_xyz1)(std::string,int) = &System::print_xyz;
void (System::*print_xyz2)(std::string,int,std::string,int) = &System::print_xyz;

std::string (System::*expt_get_xyz_v1)(int fold,std::string pbc_type,int frame) = &System::get_xyz;

int (System::*expt_Find_Angle_v1)(int,int) = &System::Find_Angle;
int (System::*expt_Find_Angle_v2)(int,int,int) = &System::Find_Angle;


void (System::*expt_init_fragment_velocities_v1)(double Temp, Random& rnd) = &System::init_fragment_velocities;
void (System::*expt_init_fragment_velocities_v2)(double Temp,VECTOR TOT_P,VECTOR TOT_L, Random& rnd) = &System::init_fragment_velocities;

void (System::*expt_init_atom_velocities_v1)(double Temp, Random& rnd) = &System::init_atom_velocities;
void (System::*expt_init_atom_velocities_v2)(double Temp,VECTOR TOT_P, Random& rnd) = &System::init_atom_velocities;

void (System::*expt_ROTATE_FRAGMENT_v1)(double phi, const VECTOR& dir, int fr_id, const VECTOR& pivot) = &System::ROTATE_FRAGMENT;
void (System::*expt_ROTATE_FRAGMENT_v2)(double phi, const VECTOR& dir, int fr_id, int center_indx) = &System::ROTATE_FRAGMENT;
void (System::*expt_ROTATE_FRAGMENT_v3)(double phi, const VECTOR& dir, int fr_id) = &System::ROTATE_FRAGMENT;
void (System::*expt_ROTATE_FRAGMENT_v4)(const MATRIX3x3& rot, int fr_id, const VECTOR& pivot) = &System::ROTATE_FRAGMENT;
void (System::*expt_ROTATE_FRAGMENT_v5)(const MATRIX3x3& rot, int fr_id, int center_indx) = &System::ROTATE_FRAGMENT;
void (System::*expt_ROTATE_FRAGMENT_v6)(const MATRIX3x3& rot, int fr_id) = &System::ROTATE_FRAGMENT;


  class_<System>("System",init<>())
      .def("__copy__", &generic__copy__<System>)
      .def("__deepcopy__", &generic__deepcopy__<System>)
      .def(init<const System&>())

      .def_readwrite("name", &System::name)
      .def_readwrite("id",&System::id)
      .def_readwrite("mass",&System::mass)
      .def_readwrite("Nf_t",&System::Nf_t)
      .def_readwrite("Nf_r",&System::Nf_r)

//      .def("set",&System::set)
      .def("show_info",&System::show_info)

      .def_readwrite("Number_of_atoms",&System::Number_of_atoms)
      .def_readwrite("Number_of_bonds",&System::Number_of_bonds)
      .def_readwrite("Number_of_angles",&System::Number_of_angles)
      .def_readwrite("Number_of_dihedrals",&System::Number_of_dihedrals)
      .def_readwrite("Number_of_impropers",&System::Number_of_impropers)
      .def_readwrite("Number_of_pairs",&System::Number_of_pairs)
      .def_readwrite("Number_of_fragments",&System::Number_of_fragments)
      .def_readwrite("Number_of_rings",&System::Number_of_rings)
      .def_readwrite("Number_of_molecules",&System::Number_of_molecules)

      .def_readwrite("Atoms",&System::Atoms)
      .def_readwrite("Bonds",&System::Bonds)
      .def_readwrite("Angles",&System::Angles)
      .def_readwrite("Dihedrals",&System::Dihedrals)
      .def_readwrite("Impropers",&System::Impropers)
      .def_readwrite("Pairs",&System::Pairs)
      .def_readwrite("Fragments",&System::Fragments)
      .def_readwrite("Rings",&System::Rings)
      .def_readwrite("Molecules",&System::Molecules)

      .def_readwrite("Frag_bonds",&System::Frag_bonds)
      .def_readwrite("Frag_angles",&System::Frag_angles)
      .def_readwrite("Frag_dihedrals",&System::Frag_dihedrals)
      .def_readwrite("Frag_impropers",&System::Frag_impropers)
      .def_readwrite("Frag_pairs",&System::Frag_pairs)
      .def_readwrite("Surface_atoms",&System::Surface_atoms)

      .def_readwrite("Number_of_frag_bonds",&System::Number_of_frag_bonds)
      .def_readwrite("Number_of_frag_angles",&System::Number_of_frag_angles)
      .def_readwrite("Number_of_frag_dihedrals",&System::Number_of_frag_dihedrals)
      .def_readwrite("Number_of_frag_impropers",&System::Number_of_frag_impropers)
      .def_readwrite("Number_of_frag_pairs",&System::Number_of_frag_pairs)
      .def_readwrite("Number_of_surface_atoms",&System::Number_of_surface_atoms)

      .def_readwrite("Box",&System::Box)
      .def_readwrite("Box_origin",&System::Box_origin)
 


  //---------- Defined in System_methods.cpp  ------------------

      .def("show_atoms",&System::show_atoms)
      .def("show_bonds",&System::show_bonds)
      .def("show_angles",&System::show_angles)
      .def("show_dihedrals",&System::show_dihedrals)
      .def("show_impropers",&System::show_impropers)
      .def("show_pairs",&System::show_pairs)
      .def("show_frag_bonds",&System::show_frag_bonds)
      .def("show_frag_angles",&System::show_frag_angles)
      .def("show_frag_dihedrals",&System::show_frag_dihedrals)
      .def("show_frag_impropers",&System::show_frag_impropers)
      .def("show_frag_pairs",&System::show_frag_pairs)
      .def("show_fragments",&System::show_fragments)
      .def("show_rings",&System::show_rings)
      .def("show_molecules",&System::show_molecules)


      .def("get_atom_index_by_atom_id", &System::get_atom_index_by_atom_id)
      .def("get_fragment_index_by_fragment_id", &System::get_fragment_index_by_fragment_id)
      .def("get_molecule_index_by_molecule_id", &System::get_molecule_index_by_molecule_id)
      .def("Find_Bond", &System::Find_Bond)
      .def("Find_Frag_Pair", &System::Find_Frag_Pair)
      .def("Find_Angle", expt_Find_Angle_v1)
      .def("Find_Angle", expt_Find_Angle_v2)
      .def("Find_Dihedral", &System::Find_Dihedral)
      .def("Find_Improper", &System::Find_Improper)
      .def("is_12pair", &System::is_12pair)
      .def("is_13pair", &System::is_13pair)
      .def("is_14pair", &System::is_14pair)
      .def("is_group_pair", &System::is_group_pair)

      

  //----------- Defined in System_methods1.cpp -----------

      .def("Generate_Connectivity_Matrix", &System::Generate_Connectivity_Matrix)
      .def("Assign_Rings", &System::Assign_Rings)
      .def("DIVIDE_GRAPH", &System::DIVIDE_GRAPH)


  //----------- Defined in System_methods2.cpp ------------------

      .def("update_max_id",&System::update_max_id)
      .def("CREATE_ATOM",CREATE_ATOM2)
      .def("LINK_ATOMS",LINK_ATOMS2)
      .def("GROUP_ATOMS",&System::GROUP_ATOMS)
      .def("UPDATE_FRAG_TOPOLOGY", &System::UPDATE_FRAG_TOPOLOGY)
      .def("ADD_ATOM_TO_FRAGMENT", &System::ADD_ATOM_TO_FRAGMENT)
      .def("CREATE_BONDS", &System::CREATE_BONDS)
      .def("CLONE_MOLECULE", &System::CLONE_MOLECULE)



  //----------- Defined in System_methods3.cpp ------------------

      .def("move_atom_by_index",&System::move_atom_by_index)
      .def("move_fragment_by_index",&System::move_fragment_by_index)
      .def("move_molecule_by_index",&System::move_molecule_by_index)

      .def("update_atoms_for_fragment", &System::update_atoms_for_fragment)
      .def("update_fragments_for_molecule", &System::update_fragments_for_molecule)
      .def("update_atoms_for_molecule", &System::update_atoms_for_molecule)

      .def("rotate_atoms_of_fragment", &System::rotate_atoms_of_fragment)
      .def("rotate_fragments_of_molecule", &System::rotate_fragments_of_molecule)
      .def("rotate_atoms_of_molecule", &System::rotate_atoms_of_molecule)

      .def("TRANSLATE_ATOM", &System::TRANSLATE_ATOM)
      .def("TRANSLATE_FRAGMENT", &System::TRANSLATE_FRAGMENT)
      .def("TRANSLATE_MOLECULE", &System::TRANSLATE_MOLECULE)
      .def("ROTATE_FRAGMENT", expt_ROTATE_FRAGMENT_v1)
      .def("ROTATE_FRAGMENT", expt_ROTATE_FRAGMENT_v2)
      .def("ROTATE_FRAGMENT", expt_ROTATE_FRAGMENT_v3)
      .def("ROTATE_FRAGMENT", expt_ROTATE_FRAGMENT_v4)
      .def("ROTATE_FRAGMENT", expt_ROTATE_FRAGMENT_v5)
      .def("ROTATE_FRAGMENT", expt_ROTATE_FRAGMENT_v6)
      .def("ROTATE_MOLECULE", &System::ROTATE_MOLECULE)



  //----------- Defined in System_methods4.cpp ------------------

      .def("determine_functional_groups",&System::determine_functional_groups)


  //----------- Defined in System_methods5.cpp (extractors/converters) -------

      .def("extract_atomic_q", &System::extract_atomic_q)
      .def("set_atomic_q", &System::set_atomic_q)

      .def("extract_atomic_p", &System::extract_atomic_p)
      .def("set_atomic_p", &System::set_atomic_p)

      .def("extract_atomic_v", &System::extract_atomic_v)
      .def("set_atomic_v", &System::set_atomic_v)

      .def("extract_atomic_f", &System::extract_atomic_f)
      .def("set_atomic_f", &System::set_atomic_f)

      .def("extract_atomic_mass", &System::extract_atomic_mass)
      .def("set_atomic_mass", &System::set_atomic_mass)


      .def("extract_fragment_q", &System::extract_fragment_q)
      .def("set_fragment_q", &System::set_fragment_q)

      .def("extract_fragment_p", &System::extract_fragment_p)
      .def("set_fragment_p", &System::set_fragment_p)

      .def("extract_fragment_v", &System::extract_fragment_v)
      .def("set_fragment_v", &System::set_fragment_v)

      .def("extract_fragment_f", &System::extract_fragment_f)
      .def("set_fragment_f", &System::set_fragment_f)

      .def("extract_fragment_mass", &System::extract_fragment_mass)
      .def("set_fragment_mass", &System::set_fragment_mass)


  //------------- Defined in System_methods6.cpp ------------------
      .def("cool_atoms", &System::cool_atoms)
      .def("cool_fragments", &System::cool_fragments)
      .def("cool", &System::cool)

      .def("zero_atom_forces",&System::zero_atom_forces)
      .def("zero_fragment_forces",&System::zero_fragment_forces)
      .def("zero_fragment_torques",&System::zero_fragment_torques)
      .def("zero_forces",&System::zero_forces)
      .def("zero_forces_and_torques",&System::zero_forces_and_torques)
      .def("update_fragment_forces",&System::update_fragment_forces)
      .def("update_fragment_torques",&System::update_fragment_torques)
      .def("update_fragment_forces_and_torques",&System::update_fragment_forces_and_torques)

      .def("save_forces", &System::save_forces)
      .def("save_torques", &System::save_torques)
      .def("load_forces", &System::load_forces)
      .def("load_torques", &System::load_torques)
    
      .def("save_stress", &System::save_stress)
      .def("increment_stress", &System::increment_stress)
      .def("save_respa_state", &System::save_respa_state)
      .def("load_respa_state", &System::load_respa_state)

      .def("init_fragments",&System::init_fragments)
      .def("init_molecules",&System::init_molecules)

      .def("init_box_origin", &System::init_box_origin)
      .def("init_box",init_box1)
      .def("init_box",init_box2)
      .def("init_box",init_box3)
      .def("apply_frag_pbc",&System::apply_frag_pbc)

      .def("fix_fragment_translation", &System::fix_fragment_translation)
      .def("fix_fragment_rotation", &System::fix_fragment_rotation)
      .def("fix_fragment", &System::fix_fragment)

      .def("ekin_tr",&System::ekin_tr)
      .def("ekin_tr_atom",&System::ekin_tr_atom)
      .def("ekin_tr_int",&System::ekin_tr_int)
      .def("ekin_rot",&System::ekin_rot)
      .def("volume",&System::volume)
      .def("pressure_tensor", &System::pressure_tensor)

      .def("init_fragment_velocities", expt_init_fragment_velocities_v1)
      .def("init_fragment_velocities", expt_init_fragment_velocities_v2)

      .def("init_atom_velocities", expt_init_atom_velocities_v1)
      .def("init_atom_velocities", expt_init_atom_velocities_v2)



  //---------------- Defined in System_methods7.cpp -----------------
      .def("print_ent",print_ent1)
      .def("print_ent",print_ent2)
      .def("print_ent",print_ent3)
      .def("print_ent",print_ent4)

      .def("print_xyz",print_xyz1)
      .def("print_xyz",print_xyz2)

      .def("get_xyz", expt_get_xyz_v1);

  ;



}// export_Chemsys_objects




#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygchemsys){
#else
BOOST_PYTHON_MODULE(libchemsys){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

//  export_Mathematics_objects();
  export_Chemsys_objects();

}


}// namespace libchemsys
}// namespace libchemobjects

}// liblibra

