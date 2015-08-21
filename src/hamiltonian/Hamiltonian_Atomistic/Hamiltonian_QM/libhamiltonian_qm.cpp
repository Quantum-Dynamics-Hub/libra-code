#include <memory> // for std::auto_ptr<>
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "libhamiltonian_qm.h"

using namespace boost::python;

namespace libhamiltonian{
namespace libhamiltonian_atomistic{
namespace libhamiltonian_qm{


void export_Hamiltonian_QM_objects(){

//  export_model_parameters_objects();




}


#ifdef CYGWIN
BOOST_PYTHON_MODULE(cyghamiltonian_qm){
#else
BOOST_PYTHON_MODULE(libhamiltonian_qm){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

  export_Hamiltonian_QM_objects();

}




}// namespace libhamiltonian_qm
}// namespace libhamiltonian_atomistic
}// namespace libhamiltonian


