#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "libgraph.h"

using namespace boost::python;


namespace libmmath{
namespace libgraph{

void export_GRAPH_objects(){
}



#ifdef CYGWIN
BOOST_PYTHON_MODULE(cyggraph){
#else
BOOST_PYTHON_MODULE(libgraph){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();
  export_GRAPH_objects();

}



}// libgraph
}// libmmath

