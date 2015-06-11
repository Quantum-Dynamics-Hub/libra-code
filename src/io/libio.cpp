#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "libio.h"


namespace libio{


void export_io_objects(){

// Nothing to export to Python


}// export_io_objects()



#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygio){
#else
BOOST_PYTHON_MODULE(libio){
#endif

  export_io_objects();

}



}// namespace libio




