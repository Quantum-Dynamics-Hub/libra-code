#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "libwfcgrid.h"

using namespace boost::python;
//using namespace libdyn::libelectronic;


namespace libdyn{
namespace libwfcgrid{


void export_Wfcgrid_objects(){

//  def("SAC_Ham", expt_SAC_Ham1);
//  void (Electronic::*expt_propagate_electronic)(double,Hamiltonian&) = &Electronic::propagate_electronic;


  class_<Wfcgrid>("Wfcgrid",init<double, double, double, double, double, double, int>())
      .def("__copy__", &generic__copy__<Wfcgrid>)
      .def("__deepcopy__", &generic__deepcopy__<Wfcgrid>)

      .def("init_wfc", &Wfcgrid::init_wfc)
      .def("print_map", &Wfcgrid::print_map)

      .def("update_potential", &Wfcgrid::update_potential)
      .def("update_propagator", &Wfcgrid::update_propagator)
      .def("update_propagator_K", &Wfcgrid::update_propagator_K)

      .def("propagate_exact_2D", &Wfcgrid::propagate_exact_2D)
 
  ;



}


#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygwfcgrid){
#else
BOOST_PYTHON_MODULE(libwfcgrid){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

//  export_Mathematics_objects();
  export_Wfcgrid_objects();

}


}// namespace libwfcgrid
}// namespace libdyn

