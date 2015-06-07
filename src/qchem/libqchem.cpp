#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "libqchem.h"
using namespace boost::python;


namespace libqchem{

using namespace libmolint;
using namespace libqobjects;

void export_Qchem_objects(){

  export_molint_objects();
  export_qobjects_objects();

}// export_Qchem_objects()




#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygqchem){
#else
BOOST_PYTHON_MODULE(libqchem){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

  export_Qchem_objects();

}


}// libqchem

