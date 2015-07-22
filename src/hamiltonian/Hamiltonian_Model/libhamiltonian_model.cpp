#include <memory> // for std::auto_ptr<>
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "libhamiltonian_model.h"

using namespace boost::python;
//using namespace libhamiltonian;

namespace libhamiltonian{
namespace libhamiltonian_model{


void export_hamiltonian_model_objects(){

  boost::python::list (*expt_SAC_Ham1)(double,boost::python::list) = &SAC_Ham;
  boost::python::list (*expt_DAC_Ham1)(double,boost::python::list) = &DAC_Ham;
  boost::python::list (*expt_ECWR_Ham1)(double,boost::python::list) = &ECWR_Ham;
  boost::python::list (*expt_Marcus_Ham1)(double,boost::python::list) = &Marcus_Ham;
  boost::python::list (*expt_SEXCH_Ham1)(double,boost::python::list) = &SEXCH_Ham;
  boost::python::list (*expt_Rabi2_Ham1)(double,boost::python::list) = &Rabi2_Ham;

  def("SAC_Ham", expt_SAC_Ham1);
  def("DAC_Ham", expt_DAC_Ham1);
  def("ECWR_Ham", expt_ECWR_Ham1);
  def("Marcus_Ham", expt_Marcus_Ham1);
  def("SEXCH_Ham", expt_SEXCH_Ham1);
  def("Rabi2_Ham", expt_Rabi2_Ham1);


//  void (Hamiltonian_Model::*expt_set_params_v1)(boost::python::list) = &Hamiltonian_Model::set_params;
//  void (Hamiltonian_Model::*set_q)(boost::python::list) = &Hamiltonian_Model::set_q;
//  void (Hamiltonian_Model::*set_v)(boost::python::list) = &Hamiltonian_Model::set_v;



  class_<Hamiltonian_Model, bases<Hamiltonian> >("Hamiltonian_Model",init<int>())
      .def("__copy__", &generic__copy__<Hamiltonian_Model>)
      .def("__deepcopy__", &generic__deepcopy__<Hamiltonian_Model>)

//      .def("set_params", expt_set_params_v1)
//      .def("set_rep", &Hamiltonian_Model::set_rep)
//      .def("set_q", set_q)
//      .def("set_v", set_v)

//      .def("compute",          &Hamiltonian_Model::compute)
      .def("compute_diabatic", &Hamiltonian_Model::compute_diabatic)
      .def("compute_adiabatic",&Hamiltonian_Model::compute_adiabatic)

//      .def("H", &Hamiltonian_Model::H)
//      .def("dHdq", &Hamiltonian_Model::dHdq)
//      .def("Hvib", &Hamiltonian_Model::Hvib)
//      .def("D", &Hamiltonian_Model::D)
//      .def("nac", &Hamiltonian_Model::nac)
 
  ;




}


#ifdef CYGWIN
BOOST_PYTHON_MODULE(cyghamiltonian_model){
#else
BOOST_PYTHON_MODULE(libhamiltonian_model){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

//  export_Mathematics_objects();
  export_hamiltonian_model_objects();

}


}// namespace libhamiltonian_model
}// namespace libhamiltonian


