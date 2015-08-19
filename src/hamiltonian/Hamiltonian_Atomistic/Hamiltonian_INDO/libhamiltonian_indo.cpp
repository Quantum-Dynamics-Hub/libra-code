#include <memory> // for std::auto_ptr<>
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "libhamiltonian_indo.h"

using namespace boost::python;

namespace libhamiltonian{
namespace libhamiltonian_atomistic{
namespace libhamiltonian_indo{


void export_Hamiltonian_INDO_objects(){

//  export_forcefield_objects();

  // Also export self classes:
//  bool (listHamiltonian_INDO::*expt_is_active_v1)(Atom&,Atom&) = &listHamiltonian_MM::is_active;

  class_<Hamiltonian_INDO>("Hamiltonian_INDO",init<>())
      .def("__copy__", &generic__copy__<Hamiltonian_INDO>)
      .def("__deepcopy__", &generic__deepcopy__<Hamiltonian_INDO>)

//      .def("activate", &Hamiltonian_MM::activate)

  ;

  class_<listHamiltonian_INDO>("listHamiltonian_INDO",init<>())
      .def("__copy__", &generic__copy__<listHamiltonian_INDO>)
      .def("__deepcopy__", &generic__deepcopy__<listHamiltonian_INDO>)
//      .def("is_active", expt_is_active_v3)
  ;



}


#ifdef CYGWIN
BOOST_PYTHON_MODULE(cyghamiltonian_indo){
#else
BOOST_PYTHON_MODULE(libhamiltonian_indo){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

//  export_Mathematics_objects();
  export_Hamiltonian_INDO_objects();

}




}// namespace libhamiltonian_indo
}// namespace libhamiltonian_atomistic
}// namespace libhamiltonian


