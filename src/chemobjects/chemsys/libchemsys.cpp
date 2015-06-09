#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "libchemsys.h"

using namespace boost::python;


namespace libchemobjects{
namespace libchemsys{


void export_Chemsys_objects(){

//void (System::*CREATE_ATOM1)()      = &System::CREATE_ATOM;
void (System::*CREATE_ATOM2)(Atom)  = &System::CREATE_ATOM;
void (System::*LINK_ATOMS1)(Atom&,Atom&) = &System::LINK_ATOMS;
void (System::*LINK_ATOMS2)(int,int) = &System::LINK_ATOMS;

//double (System::*energy1)()         = &System::energy;
//double (System::*energy2)(std::string) = &System::energy;

//void (System::*init_box1)()         = &System::init_box;
//void (System::*init_box2)(double,double,double) = &System::init_box;
//void (System::*init_box3)(VECTOR,VECTOR,VECTOR) = &System::init_box;

void (System::*show_interactions1)()            = &System::show_interactions;
void (System::*show_interactions2)(std::string) = &System::show_interactions;


void (System::*print_ent1)(std::string) = &System::print_ent;
void (System::*print_ent2)(std::string,int,std::string) = &System::print_ent;
void (System::*print_ent3)(std::string,boost::python::list) = &System::print_ent;
void (System::*print_ent4)(std::string,boost::python::list,int,std::string) = &System::print_ent;

void (System::*print_xyz1)(std::string,int) = &System::print_xyz;
void (System::*print_xyz2)(std::string,int,std::string,int) = &System::print_xyz;



  class_<System>("System",init<>())
      .def("__copy__", &generic__copy__<System>)
      .def("__deepcopy__", &generic__deepcopy__<System>)
      .def(init<const System&>())

      .def_readwrite("id",&System::id)

//      .def("set",&System::set)
      .def_readwrite("Number_of_atoms",&System::Number_of_atoms)
      .def_readwrite("Number_of_bonds",&System::Number_of_bonds)
      .def_readwrite("Number_of_angles",&System::Number_of_angles)
      .def_readwrite("Number_of_dihedrals",&System::Number_of_dihedrals)
      .def_readwrite("Number_of_impropers",&System::Number_of_impropers)
      .def_readwrite("Number_of_pairs",&System::Number_of_pairs)
      .def_readwrite("Number_of_fragments",&System::Number_of_fragments)
      .def_readwrite("Number_of_rings",&System::Number_of_rings)
      .def_readwrite("Number_of_molecules",&System::Number_of_molecules)

      .def_readwrite("Box",&System::Box)
      .def_readwrite("Box_origin",&System::Box_origin)
 
      .def_readwrite("Atoms",&System::Atoms)

      .def("show_info",&System::show_info)
      .def("move_atom_by_index",&System::move_atom_by_index)
      .def("move_fragment_by_index",&System::move_fragment_by_index)
      .def("move_molecule_by_index",&System::move_molecule_by_index)

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
      .def("show_interactions",show_interactions1)
      .def("show_interactions",show_interactions2)

//      .def("CREATE_ATOM",CREATE_ATOM1)
      .def("CREATE_ATOM",CREATE_ATOM2)
      .def("LINK_ATOMS",LINK_ATOMS2)
      .def("GROUP_ATOMS",&System::GROUP_ATOMS)

      .def("determine_functional_groups",&System::determine_functional_groups)
//      .def("set_interactions_for_atoms",&System::set_interactions_for_atoms)
//      .def("set_respa_types",&System::set_respa_types)
//      .def("show_interactions_statistics",&System::show_interactions_statistics)
/*
      .def("zero_atom_forces",&System::zero_atom_forces)
      .def("zero_fragment_forces",&System::zero_fragment_forces)
      .def("zero_fragment_torques",&System::zero_fragment_torques)
      .def("zero_forces",&System::zero_forces)
      .def("zero_forces_and_torques",&System::zero_forces_and_torques)
      .def("init_fragments",&System::init_fragments)
      .def("init_molecules",&System::init_molecules)
      .def("init_box",init_box1)
      .def("init_box",init_box2)
      .def("init_box",init_box3)
      .def("apply_frag_pbc",&System::apply_frag_pbc)

      .def("fix_fragment_translation", &System::fix_fragment_translation)
      .def("fix_fragment_rotation", &System::fix_fragment_rotation)
      .def("fix_fragment", &System::fix_fragment)


      .def("energy",energy1)
      .def("energy",energy2)
      .def("ekin_tr",&System::ekin_tr)
      .def("ekin_rot",&System::ekin_rot)
      .def("volume",&System::volume)
*/

      .def("print_ent",print_ent1)
      .def("print_ent",print_ent2)
      .def("print_ent",print_ent3)
      .def("print_ent",print_ent4)

      .def("print_xyz",print_xyz1)
      .def("print_xyz",print_xyz2)

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



