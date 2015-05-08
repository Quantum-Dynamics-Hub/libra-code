#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "libpot.h"

using namespace libpot;
using namespace boost::python;


void export_Pot_objects(){


  def("max_vector", max_vector);
  def("apply_pbc", apply_pbc);


} // export_Pot_objects()



#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygpot){
#else
BOOST_PYTHON_MODULE(libpot){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

  export_Pot_objects();

}


