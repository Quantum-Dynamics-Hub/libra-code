#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "libmol.h"

using namespace boost::python;


namespace libchemobjects{
namespace libmol{


void export_Mol_objects(){


  void (Atom::*expt_save_v1)(boost::property_tree::ptree& pt,std::string path) = &Atom::save;
  void (Atom::*expt_save_v2)(std::string path) = &Atom::save;


  class_<Atom>("Atom",init<Universe&>())
      .def(init<Universe&, boost::python::dict>())
      .def("__copy__", &generic__copy__<Atom>)
      .def("__deepcopy__", &generic__deepcopy__<Atom>)


      // Topology
      .def_readwrite("Atom_id",&Atom::Atom_id)
      .def_readwrite("globAtom_Index",&Atom::globAtom_Index)
      .def_readwrite("locAtom_Index",&Atom::locAtom_Index)
      .def_readwrite("globGroup_Index",&Atom::globGroup_Index)
      .def_readwrite("globMolecule_Index",&Atom::globMolecule_Index)
      .def_readwrite("globAtom_Adjacent_Atoms",&Atom::globAtom_Adjacent_Atoms)

      // Geometry
      .def_readwrite("Atom_RB",&Atom::Atom_RB)
      .def_readwrite("Atom_RB_old",&Atom::Atom_RB_old)
      .def_readwrite("Atom_displ2",&Atom::Atom_displ2)

      // Properties
      .def_readwrite("Atom_element",&Atom::Atom_element)
      .def_readwrite("Atom_atomic_radius",&Atom::Atom_atomic_radius)
      .def_readwrite("Atom_charge",&Atom::Atom_charge)
      .def_readwrite("Atom_electronegativity",&Atom::Atom_electronegativity)

      .def_readwrite("Atom_formal_charge",&Atom::Atom_formal_charge)
      .def_readwrite("Atom_coordination",&Atom::Atom_coordination)
      .def_readwrite("Atom_functional_group",&Atom::Atom_functional_group)
      .def_readwrite("Atom_min_ring_size",&Atom::Atom_min_ring_size)

      // FF types
      .def_readwrite("Atom_ff_type",&Atom::Atom_ff_type)
//      .def_readwrite("Atom_ff_int_type",&Atom::Atom_ff_int_type)
//      .def_readwrite("Atom_is_surface_atom",&Atom::Atom_is_surface_atom)
//      .def_readwrite("Atom_surface_index",&Atom::Atom_surface_index)
//      .def_readwrite("Atom_is_basis_atom",&Atom::Atom_is_basis_atom)
//      .def_readwrite("Atom_is_C60_CT",&Atom::Atom_is_C60_CT)

      .def("set",&Atom::set)
      .def("show_info",&Atom::show_info)
      .def("save",expt_save_v2)
      .def("load",&Atom::load)
  ;

  class_<std::vector<Atom> >("AtomList")
      .def(vector_indexing_suite<std::vector<Atom> >())
  ;


//  void (Atom::*expt_save_v1)(boost::property_tree::ptree& pt,std::string path) = &Atom::save;
//  void (Atom::*expt_save_v2)(std::string path) = &Atom::save;


  class_<Group>("Group",init<>())
      .def("__copy__", &generic__copy__<Group>)
      .def("__deepcopy__", &generic__deepcopy__<Group>)

      // Topology
      .def_readwrite("globGroup_Size",&Group::globGroup_Size)
      .def_readwrite("locGroup_Size",&Group::locGroup_Size)
      .def_readwrite("Group_Size",&Group::Group_Size)
      .def_readwrite("globAtom_Index",&Group::globAtom_Index)
      .def_readwrite("locAtom_Index",&Group::locAtom_Index)
      .def_readwrite("globGroup_Index",&Group::globGroup_Index)
      .def_readwrite("locGroup_Index",&Group::locGroup_Index)
      .def_readwrite("globMolecule_Index",&Group::globMolecule_Index)

      .def_readwrite("Group_name",&Group::Group_name)
      .def_readwrite("Group_id",&Group::Group_id)
      .def_readwrite("Group_radius",&Group::Group_radius)
      .def_readwrite("Group_RB",&Group::Group_RB)
      .def_readwrite("Group_ff_type",&Group::Group_ff_type)

      .def_readwrite("Group_bond_order",&Group::Group_bond_order)
      .def_readwrite("Group_bond_alpha",&Group::Group_bond_alpha)

      .def("set",&Group::set)
      .def("show_inf",&Group::show_info)

  ;

  class_<std::vector<Group> >("GroupList")
      .def(vector_indexing_suite<std::vector<Group> >())
  ;


  class_<Molecule>("Molecule",init<>())
      .def("__copy__", &generic__copy__<Molecule>)
      .def("__deepcopy__", &generic__deepcopy__<Molecule>)

      // Topology
      .def_readwrite("globMolecule_Size",&Molecule::globMolecule_Size)
      .def_readwrite("locMolecule_Size",&Molecule::locMolecule_Size)
      .def_readwrite("Molecule_Size",&Molecule::Molecule_Size)
      .def_readwrite("globAtom_Index",&Molecule::globAtom_Index)
      .def_readwrite("locAtom_Index",&Molecule::locAtom_Index)
      .def_readwrite("globMolecule_Index",&Molecule::globMolecule_Index)
      .def_readwrite("locMolecule_Index",&Molecule::locMolecule_Index)

      .def_readwrite("Molecule_Number_of_bonds",&Molecule::Molecule_Number_of_bonds)
      .def_readwrite("Molecule_Number_of_angles",&Molecule::Molecule_Number_of_angles)
      .def_readwrite("Molecule_Number_of_dihedrals",&Molecule::Molecule_Number_of_dihedrals)
      .def_readwrite("Molecule_Number_of_impropers",&Molecule::Molecule_Number_of_impropers)

      .def_readwrite("Molecule_name",&Molecule::Molecule_name)
      .def_readwrite("Molecule_id",&Molecule::Molecule_id)
      .def_readwrite("Molecule_RB",&Molecule::Molecule_RB)

  ;

  class_<std::vector<Molecule> >("MoleculeList")
      .def(vector_indexing_suite<std::vector<Molecule> >())
  ;



}// export_Mol_objects




#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygmol){
#else
BOOST_PYTHON_MODULE(libmol){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

//  export_Mathematics_objects();
  export_Mol_objects();

}


}// namespace libmol
}// namespace libchemobjects



