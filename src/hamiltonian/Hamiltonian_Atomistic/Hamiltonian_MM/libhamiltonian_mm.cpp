#include <memory> // for std::auto_ptr<>
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "libhamiltonian_mm.h"

using namespace boost::python;
//using namespace libhamiltonian;

namespace libhamiltonian{
namespace libhamiltonian_atomistic{
namespace libhamiltonian_mm{

using namespace libforcefield;


void export_Hamiltonian_MM_objects(){

  export_forcefield_objects();

  // Also export self classes:
  //  void (Hamiltonian::*set_params)(boost::python::list) = &Hamiltonian::set_params;

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


