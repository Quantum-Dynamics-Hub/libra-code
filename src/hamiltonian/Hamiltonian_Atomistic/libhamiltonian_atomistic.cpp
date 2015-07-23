#include <memory> // for std::auto_ptr<>
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "libhamiltonian_atomistic.h"

using namespace boost::python;


namespace libhamiltonian{
namespace libhamiltonian_atomistic{


void export_hamiltonian_atomistic_objects(){


//  void (Hamiltonian_Atomistic::*set_params)(boost::python::list) = &Hamiltonian_Atomistic::set_params;
//  void (Hamiltonian_Atomistic::*set_q)(boost::python::list) = &Hamiltonian_Atomistic::set_q;
//  void (Hamiltonian_Atomistic::*set_v)(boost::python::list) = &Hamiltonian_Atomistic::set_v;



  class_<Hamiltonian_Atomistic, bases<Hamiltonian> >("Hamiltonian_Atomistic",init<int,int>())
      .def("__copy__", &generic__copy__<Hamiltonian_Atomistic>)
      .def("__deepcopy__", &generic__deepcopy__<Hamiltonian_Atomistic>)

/*
      .def("set_params", set_params)
      .def("set_rep", &Hamiltonian_Atomistic::set_rep)
      .def("set_q", set_q)
      .def("set_v", set_v)

      .def("compute",          &Hamiltonian_Atomistic::compute)
*/
      .def("compute_diabatic", &Hamiltonian_Atomistic::compute_diabatic)
      .def("compute_adiabatic",&Hamiltonian_Atomistic::compute_adiabatic)

/*
      .def("H", &Hamiltonian_Atomistic::H)
      .def("dHdq", &Hamiltonian_Atomistic::dHdq)
      .def("Hvib", &Hamiltonian_Atomistic::Hvib)
      .def("D", &Hamiltonian_Atomistic::D)
      .def("nac", &Hamiltonian_Atomistic::nac)
*/
 
  ;


}


#ifdef CYGWIN
BOOST_PYTHON_MODULE(cyghamiltonian_atomistic){
#else
BOOST_PYTHON_MODULE(libhamiltonian_atomistic){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

  export_hamiltonian_atomistic_objects();

}


}// namespace libhamiltonian_atomistic
}// namespace libhamiltonian


