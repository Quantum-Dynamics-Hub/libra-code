#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "libhamiltonian.h"

using namespace boost::python;
using namespace libhamiltonian;



void export_Hamiltonian_objects(){

  boost::python::list (*expt_SAC_Ham1)(double,boost::python::list) = &SAC_Ham;
  boost::python::list (*expt_DAC_Ham1)(double,boost::python::list) = &DAC_Ham;
  boost::python::list (*expt_ECWR_Ham1)(double,boost::python::list) = &ECWR_Ham;

//boost::python::list (*expt_rotate1)(double,double,double) = &expt_rotate;
//double (*expt_shift1)(double, double) = &expt_shift;
//double (*expt_scale1)(double, double) = &expt_scale;
//  def("rotate", expt_rotate1);
//  def("shift", expt_shift1);
  def("SAC_Ham", expt_SAC_Ham1);
  def("DAC_Ham", expt_DAC_Ham1);
  def("ECWR_Ham", expt_ECWR_Ham1);


}


#ifdef CYGWIN
BOOST_PYTHON_MODULE(cyghamiltonian){
#else
BOOST_PYTHON_MODULE(libhamiltonian){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

//  export_Mathematics_objects();
  export_Hamiltonian_objects();

}


