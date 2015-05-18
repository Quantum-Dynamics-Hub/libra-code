#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "libdyn.h"
using namespace boost::python;


namespace libmol{

using namespace libatom;

void export_Mol_objects(){

  export_Atom_objects();

}// export_Mol_objects()




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


}// libmol

