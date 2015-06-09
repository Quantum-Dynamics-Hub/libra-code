#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "libchemobjects.h"
using namespace boost::python;


namespace libchemobjects{

using namespace libuniverse;
using namespace libmol;
using namespace libchemsys;


void export_chemobjects_objects(){

  export_Universe_objects();
  export_Mol_objects();
  export_Chemsys_objects();

}// export_chemobjects_objects()




#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygchemobjects){
#else
BOOST_PYTHON_MODULE(libchemobjects){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();
//  export_Mathematics_objects();
  export_chemobjects_objects();

}


}// libchemobjects


