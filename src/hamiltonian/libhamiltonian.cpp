#include <memory> // for std::auto_ptr<>
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "libhamiltonian.h"

using namespace boost::python;
//using namespace libhamiltonian;

namespace libhamiltonian{

using namespace libhamiltonian_generic;
using namespace libhamiltonian_model;

void export_Hamiltonian_objects(){


  export_hamiltonian_generic_objects();
  export_hamiltonian_model_objects();

}


#ifdef CYGWIN
BOOST_PYTHON_MODULE(cyghamiltonian){
#else
BOOST_PYTHON_MODULE(libhamiltonian){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

//  export_Mathematics_objects();
  export_Hamiltonian_objects();

}


}// namespace libhamiltonian


