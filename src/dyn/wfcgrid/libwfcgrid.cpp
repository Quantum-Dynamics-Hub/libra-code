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
  boost::python::list (Wfcgrid::*expt_absorb_1D)(double dL) = &Wfcgrid::absorb_1D;


  class_<Wfcgrid>("Wfcgrid",init<double,double,double,int>())
      .def(init<double, double, double, double, double, double, int>())
      .def("__copy__", &generic__copy__<Wfcgrid>)
      .def("__deepcopy__", &generic__deepcopy__<Wfcgrid>)

      .def("init_wfc_1D", &Wfcgrid::init_wfc_1D)
      .def("init_wfc_2D", &Wfcgrid::init_wfc_2D)

      .def("print_wfc_1D", &Wfcgrid::print_wfc_1D)
      .def("print_wfc_2D", &Wfcgrid::print_wfc_2D)
      .def("print_reci_wfc_1D", &Wfcgrid::print_reci_wfc_1D)
      .def("print_ham_1D", &Wfcgrid::print_ham_1D)
      .def("print_expH_1D", &Wfcgrid::print_expH_1D)
      .def("print_expK_1D", &Wfcgrid::print_expK_1D)

      .def("print_populations_1D", &Wfcgrid::print_populations_1D)

      .def("update_potential_1D", &Wfcgrid::update_potential_1D)
      .def("update_potential_2D", &Wfcgrid::update_potential_2D)

      .def("update_propagator_1D", &Wfcgrid::update_propagator_1D)
      .def("update_propagator_2D", &Wfcgrid::update_propagator_2D)

      .def("update_propagator_K_1D", &Wfcgrid::update_propagator_K_1D)
      .def("update_propagator_K_2D", &Wfcgrid::update_propagator_K_2D)

      .def("propagate_exact_1D", &Wfcgrid::propagate_exact_1D)
      .def("propagate_exact_2D", &Wfcgrid::propagate_exact_2D)

      .def("absorb",expt_absorb_1D);
 
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

