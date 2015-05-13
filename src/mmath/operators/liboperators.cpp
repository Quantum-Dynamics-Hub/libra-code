#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "liboperators.h"

using namespace boost::python;


namespace libmmath{
namespace liboperators{


void export_Operators_objects(){


  boost::python::list (*expt_rotate1)(double,double,double) = &expt_rotate;
  double (*expt_shift1)(double, double) = &expt_shift;
  double (*expt_scale1)(double, double) = &expt_scale;


  def("rotate", expt_rotate1);
  def("shift", expt_shift1);
  def("scale", expt_scale1);


}


#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygoperators){
#else
BOOST_PYTHON_MODULE(liboperators){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

//  export_Mathematics_objects();
  export_Operators_objects();

}

}// liboperators
}// libmmath


