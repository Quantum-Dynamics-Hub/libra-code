#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "libdyn.h"
using namespace boost::python;


namespace libdyn{

using namespace libnuclear;
using namespace libelectronic;
using namespace librigidbody;

void export_Dyn_objects(){


}// export_Dyn_objects()




#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygdyn){
#else
BOOST_PYTHON_MODULE(libdyn){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();
//  export_Mathematics_objects();
  export_Nuclear_objects();
  export_RigidBody_objects();
  export_Electronic_objects();
  export_Dyn_objects();

}


}// libmmath

