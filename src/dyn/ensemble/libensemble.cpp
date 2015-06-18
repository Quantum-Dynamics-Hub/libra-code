#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "libensemble.h"

using namespace boost::python;


namespace libdyn{
namespace libensemble{


void export_Ensemble_objects(){

//  def("SAC_Ham", expt_SAC_Ham1);
//  void (Electronic::*expt_propagate_electronic)(double,Hamiltonian&) = &Electronic::propagate_electronic;


  class_<Ensemble>("Ensemble",init<>())
      .def(init<int,int,int>())
      .def("__copy__", &generic__copy__<Ensemble>)
      .def("__deepcopy__", &generic__deepcopy__<Ensemble>)

      .def_readwrite("ntraj", &Ensemble::ntraj)
      .def_readwrite("nnucl", &Ensemble::nnucl)
      .def_readwrite("nelec", &Ensemble::nelec)
      .def_readwrite("is_active", &Ensemble::is_active)


 
  ;



}


#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygensemble){
#else
BOOST_PYTHON_MODULE(libensemble){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

//  export_Mathematics_objects();
  export_Ensemble_objects();

}


}// namespace libensemble
}// namespace libdyn

