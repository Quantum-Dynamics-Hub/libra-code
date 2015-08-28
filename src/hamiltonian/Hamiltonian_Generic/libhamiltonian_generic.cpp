/*********************************************************************************
* Copyright (C) 2015 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/

#include <memory> // for std::auto_ptr<>
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "libhamiltonian_generic.h"

using namespace boost::python;
//using namespace libhamiltonian;

namespace libhamiltonian{
namespace libhamiltonian_generic{


void export_hamiltonian_generic_objects(){



  void (Hamiltonian::*set_params)(boost::python::list) = &Hamiltonian::set_params;
  void (Hamiltonian::*set_q)(boost::python::list) = &Hamiltonian::set_q;
  void (Hamiltonian::*set_v)(boost::python::list) = &Hamiltonian::set_v;


  class_<Hamiltonian>("Hamiltonian",init<>())
      .def("__copy__", &generic__copy__<Hamiltonian>)
      .def("__deepcopy__", &generic__deepcopy__<Hamiltonian>)

      .def("set_rep", &Hamiltonian::set_rep)
      .def("set_params", set_params)
      .def("set_q", set_q)
      .def("set_v", set_v)

      .def("compute",          &Hamiltonian::compute)
//      .def("compute_diabatic", &Hamiltonian_Model::compute_diabatic)
//      .def("compute_adiabatic",&Hamiltonian_Model::compute_adiabatic)

      .def("H", &Hamiltonian::H)
      .def("dHdq", &Hamiltonian::dHdq)
      .def("Hvib", &Hamiltonian::Hvib)
      .def("D", &Hamiltonian::D)
      .def("nac", &Hamiltonian::nac)
  ;

  class_< HamiltonianList >("HamiltonianList")
      .def(vector_indexing_suite< HamiltonianList >())
  ;


}


#ifdef CYGWIN
BOOST_PYTHON_MODULE(cyghamiltonian_generic){
#else
BOOST_PYTHON_MODULE(libhamiltonian_generic){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

//  export_Mathematics_objects();
  export_hamiltonian_generic_objects();

}



}// namespace libhamiltonian_generic
}// namespace libhamiltonian


