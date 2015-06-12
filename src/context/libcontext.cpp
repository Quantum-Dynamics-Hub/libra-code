#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "libcontext.h"

namespace libcontext{


void export_context_objects(){

  export_Context_objects();
  export_ctx_Control_Parameters_objects();

}// export_context_objects()



#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygcontext){
#else
BOOST_PYTHON_MODULE(libcontext){
#endif

  export_context_objects();

}



}// namespace libcontext




