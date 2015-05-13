#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "libelectronic.h"

using namespace boost::python;
using namespace libdyn::libelectronic;



void export_Electronic_objects(){

//  def("SAC_Ham", expt_SAC_Ham1);
//  void (Hamiltonian_Model::*set_params)(boost::python::list) = &Hamiltonian_Model::set_params;


  class_<Electronic>("Electronic",init<>())
      .def(init<int>())
      .def(init<int,int>())
      .def("__copy__", &generic__copy__<Electronic>)
      .def("__deepcopy__", &generic__deepcopy__<Electronic>)

//      .def("nac", &Hamiltonian_Model::nac)
 
  ;



}


#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygelectronic){
#else
BOOST_PYTHON_MODULE(libelectronic){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

//  export_Mathematics_objects();
  export_Electronic_objects();

}


