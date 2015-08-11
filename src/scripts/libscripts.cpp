#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "libscripts.h"
using namespace boost::python;


namespace libscripts{

using namespace libstate;


void export_scripts_objects(){

  export_state_objects();

}// export_scripts_objects()




#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygscripts){
#else
BOOST_PYTHON_MODULE(libscripts){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();
//  export_Mathematics_objects();
  export_scripts_objects();

}


}// libscripts


