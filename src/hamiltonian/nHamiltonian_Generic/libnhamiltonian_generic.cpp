/*********************************************************************************
* Copyright (C) 2015-2017 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file libhamiltonian_generic.cpp
  \brief The file implements Python export function
    
*/

#include <memory> // for std::auto_ptr<>
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "libnhamiltonian_generic.h"

/// liblibra namespace
namespace liblibra{


using namespace boost::python;

/// libhamiltonian namespace
namespace libhamiltonian{

/// libhamiltonian_generic namespace
namespace libhamiltonian_generic{


void export_nhamiltonian_generic_objects(){
/** 
  \brief Exporter of the libnhamiltonian_generic classes and functions

*/


//  void (Hamiltonian::*set_params)(boost::python::list) = &Hamiltonian::set_params;

  CMATRIX (nHamiltonian::*expt_get_d2ham_dia_v1)(int i) = &nHamiltonian::get_d2ham_dia;
  CMATRIX (nHamiltonian::*expt_get_d2ham_dia_v2)(int i, int j) = &nHamiltonian::get_d2ham_dia;

  CMATRIX (nHamiltonian::*expt_get_d2ham_adi_v1)(int i) = &nHamiltonian::get_d2ham_adi;
  CMATRIX (nHamiltonian::*expt_get_d2ham_adi_v2)(int i, int j) = &nHamiltonian::get_d2ham_adi;


  // for internal model types
  void (nHamiltonian::*expt_compute_diabatic_v1)(int model, vector<double>& q, vector<double>& params)
  = &nHamiltonian::compute_diabatic; 

  // for models defined in Python
  void (nHamiltonian::*expt_compute_diabatic_v2)(bp::object py_funct, bp::object q, bp::object params)
  = &nHamiltonian::compute_diabatic;




  class_<nHamiltonian>("nHamiltonian",init<int,int,int>())
//      .def("__copy__", &generic__copy__<Hamiltonian>)
//      .def("__deepcopy__", &generic__deepcopy__<Hamiltonian>)

      .def("set_ovlp_dia_by_ref", &nHamiltonian::set_ovlp_dia_by_ref)
      .def("set_ovlp_dia_by_val", &nHamiltonian::set_ovlp_dia_by_val)

//      .def("set_ampl_dia_by_ref", &nHamiltonian::set_ampl_dia_by_ref)
//      .def("set_ampl_dia_by_val", &nHamiltonian::set_ampl_dia_by_val)

      .def("set_dc1_dia_by_ref", &nHamiltonian::set_dc1_dia_by_ref)
      .def("set_dc1_dia_by_val", &nHamiltonian::set_dc1_dia_by_val)

      .def("set_ham_dia_by_ref", &nHamiltonian::set_ham_dia_by_ref)
      .def("set_ham_dia_by_val", &nHamiltonian::set_ham_dia_by_val)

      .def("set_nac_dia_by_ref", &nHamiltonian::set_nac_dia_by_ref)
      .def("set_nac_dia_by_val", &nHamiltonian::set_nac_dia_by_val)

      .def("set_hvib_dia_by_ref", &nHamiltonian::set_hvib_dia_by_ref)
      .def("set_hvib_dia_by_val", &nHamiltonian::set_hvib_dia_by_val)

      .def("set_d1ham_dia_by_ref", &nHamiltonian::set_d1ham_dia_by_ref)
      .def("set_d1ham_dia_by_val", &nHamiltonian::set_d1ham_dia_by_val)

      .def("set_d2ham_dia_by_ref", &nHamiltonian::set_d2ham_dia_by_ref)
      .def("set_d2ham_dia_by_val", &nHamiltonian::set_d2ham_dia_by_val)

//      .def("set_den_mat_dia_by_ref", &nHamiltonian::set_den_mat_dia_by_ref)
//      .def("set_den_mat_dia_by_val", &nHamiltonian::set_den_mat_dia_by_val)



//      .def("set_ampl_adi_by_ref", &nHamiltonian::set_ampl_adi_by_ref)
//      .def("set_ampl_adi_by_val", &nHamiltonian::set_ampl_adi_by_val)

      .def("set_dc1_adi_by_ref", &nHamiltonian::set_dc1_adi_by_ref)
      .def("set_dc1_adi_by_val", &nHamiltonian::set_dc1_adi_by_val)

      .def("set_ham_adi_by_ref", &nHamiltonian::set_ham_adi_by_ref)
      .def("set_ham_adi_by_val", &nHamiltonian::set_ham_adi_by_val)

      .def("set_nac_adi_by_ref", &nHamiltonian::set_nac_adi_by_ref)
      .def("set_nac_adi_by_val", &nHamiltonian::set_nac_adi_by_val)

      .def("set_hvib_adi_by_ref", &nHamiltonian::set_hvib_adi_by_ref)
      .def("set_hvib_adi_by_val", &nHamiltonian::set_hvib_adi_by_val)

      .def("set_d1ham_adi_by_ref", &nHamiltonian::set_d1ham_adi_by_ref)
      .def("set_d1ham_adi_by_val", &nHamiltonian::set_d1ham_adi_by_val)

      .def("set_d2ham_adi_by_ref", &nHamiltonian::set_d2ham_adi_by_ref)
      .def("set_d2ham_adi_by_val", &nHamiltonian::set_d2ham_adi_by_val)

//      .def("set_den_mat_adi_by_ref", &nHamiltonian::set_den_mat_adi_by_ref)
//      .def("set_den_mat_adi_by_val", &nHamiltonian::set_den_mat_adi_by_val)

      .def("set_basis_transform_by_ref", &nHamiltonian::set_basis_transform_by_ref)
      .def("set_basis_transform_by_val", &nHamiltonian::set_basis_transform_by_val)


      .def("get_ovlp_dia", &nHamiltonian::get_ovlp_dia)
//      .def("get_ampl_dia", &nHamiltonian::get_ampl_dia)
      .def("get_dc1_dia", &nHamiltonian::get_dc1_dia)
      .def("get_ham_dia", &nHamiltonian::get_ham_dia)
      .def("get_nac_dia", &nHamiltonian::get_nac_dia)
      .def("get_hvib_dia", &nHamiltonian::get_hvib_dia)
      .def("get_d1ham_dia", &nHamiltonian::get_d1ham_dia)
      .def("get_d2ham_dia", expt_get_d2ham_dia_v1)
      .def("get_d2ham_dia", expt_get_d2ham_dia_v2)
//      .def("get_den_mat_dia", &nHamiltonian::get_den_mat_dia)

//      .def("get_ampl_adi", &nHamiltonian::get_ampl_adi)
      .def("get_dc1_adi", &nHamiltonian::get_dc1_adi)
      .def("get_ham_adi", &nHamiltonian::get_ham_adi)
      .def("get_nac_adi", &nHamiltonian::get_nac_adi)
      .def("get_hvib_adi", &nHamiltonian::get_hvib_adi)
      .def("get_d1ham_adi", &nHamiltonian::get_d1ham_dia)
      .def("get_d2ham_adi", expt_get_d2ham_adi_v1)
      .def("get_d2ham_adi", expt_get_d2ham_adi_v2)
//      .def("get_den_mat_adi", &nHamiltonian::get_den_mat_adi)

      .def("get_basis_transform", &nHamiltonian::get_basis_transform)


      .def("compute_diabatic", expt_compute_diabatic_v1)
      .def("compute_diabatic", expt_compute_diabatic_v2)

      .def("compute_adiabatic", &nHamiltonian::compute_adiabatic)
      .def("ampl_adi2dia", &nHamiltonian::ampl_adi2dia)
      .def("ampl_dia2adi", &nHamiltonian::ampl_dia2adi)


      .def("Ehrenfest_energy_adi", &nHamiltonian::Ehrenfest_energy_adi)
      .def("Ehrenfest_energy_dia", &nHamiltonian::Ehrenfest_energy_dia)

      .def("forces_tens_adi", &nHamiltonian::forces_tens_adi)
      .def("forces_tens_dia", &nHamiltonian::forces_tens_dia)

      .def("forces_adi", &nHamiltonian::forces_adi)
      .def("forces_dia", &nHamiltonian::forces_dia)

      .def("Ehrenfest_forces_tens_adi", &nHamiltonian::Ehrenfest_forces_tens_adi)
      .def("Ehrenfest_forces_tens_dia", &nHamiltonian::Ehrenfest_forces_tens_dia)

      .def("Ehrenfest_forces_adi", &nHamiltonian::Ehrenfest_forces_adi)
      .def("Ehrenfest_forces_dia", &nHamiltonian::Ehrenfest_forces_dia)

      .def("compute_nac_dia", &nHamiltonian::compute_nac_dia)
      .def("compute_nac_adi", &nHamiltonian::compute_nac_adi)

      .def("compute_hvib_dia", &nHamiltonian::compute_hvib_dia)
      .def("compute_hvib_adi", &nHamiltonian::compute_hvib_adi)




  ;

  class_< nHamiltonianList >("nHamiltonianList")
      .def(vector_indexing_suite< nHamiltonianList >())
  ;


}


#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygnhamiltonian_generic){
#else
BOOST_PYTHON_MODULE(libnhamiltonian_generic){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

//  export_Mathematics_objects();
  export_nhamiltonian_generic_objects();

}



}// namespace libhamiltonian_generic
}// namespace libhamiltonian
}// liblibra

