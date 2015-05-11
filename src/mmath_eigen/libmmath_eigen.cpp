#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "libmmath_eigen.h"

using namespace boost::python;
using namespace libmmath;



void export_mmath_eigen_objects(){


//boost::python::list (*expt_rotate1)(double,double,double) = &expt_rotate;
//double (*expt_shift1)(double, double) = &expt_shift;
//double (*expt_scale1)(double, double) = &expt_scale;
//  def("rotate", expt_rotate1);
//  def("shift", expt_shift1);
//  def("scale", expt_scale1);


}


#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygmmath_eigen){
#else
BOOST_PYTHON_MODULE(libmmath_eigen){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

//  export_Mathematics_objects();
  export_mmath_eigen_objects();

}


