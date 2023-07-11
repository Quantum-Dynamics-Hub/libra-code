/*********************************************************************************
* Copyright (C) 2019 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file libwfcgrid2.cpp
  \brief The file implements Python export function
    
*/

#if defined(USING_PCH)
#include "../../pch.h"
#else
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#endif

#include "libwfcgrid2.h"

/// liblibra namespace
namespace liblibra{

using namespace boost::python;

/// libdyn namespace
namespace libdyn{

/// libwfcgrid namespace
namespace libwfcgrid2{


void export_Wfcgrid2_objects(){
/** 
  \brief Exporter of libwfcgrid2 classes and functions

*/

  // Gaussian wavepackets
  CMATRIX (*expt_Gaussian_v1)
  (vector<double>& x, vector<double>& x0, vector<double>& px0, vector<double>& dx, int init_state, int nstates, complex<double> scl) = &Gaussian;

  // Harmonic oscillator
  CMATRIX (*expt_HO_v1)
  (vector<double>& x, vector<double>& x0, vector<double>& px0, 
  vector<double>& alpha, int init_state, int nstates, vector<int>& nu, complex<double> scl) = &HO;

  def("Gaussian", expt_Gaussian_v1);
  def("HO", expt_HO_v1);


  /*
  boost::python::list (Wfcgrid::*expt_absorb_1D)(double dL) = &Wfcgrid::absorb_1D;

  void (Wfcgrid::*expt_update_potential_1D_v1)(Hamiltonian& ham) = &Wfcgrid::update_potential_1D;
  void (Wfcgrid::*expt_update_potential_1D_v2)(bp::object py_funct, bp::object params) = &Wfcgrid::update_potential_1D;
  void (Wfcgrid::*expt_update_potential_2D_v1)(Hamiltonian& ham) = &Wfcgrid::update_potential_2D;
  void (Wfcgrid::*expt_update_potential_2D_v2)(bp::object py_funct, bp::object params) = &Wfcgrid::update_potential_2D;

  void (Wfcgrid::*expt_update_propagator_K_2D_v1)(double dt, double m1, double m2) = &Wfcgrid::update_propagator_K_2D;
  void (Wfcgrid::*expt_update_propagator_K_2D_v2)(double dt, double m0) = &Wfcgrid::update_propagator_K_2D;


  // Arbitrary wavefunction
  void (Wfcgrid::*expt_init_wfc_1D_ARB_v1)
  (bp::object py_funct, bp::object params) = &Wfcgrid::init_wfc_1D_ARB;
  */

  // Gaussian wavepacket
  void (Wfcgrid2::*expt_add_wfc_Gau_v1)
  (vector<double>& x0, vector<double>& px0, vector<double>& dx0, int init_state, complex<double> weight, int rep) = &Wfcgrid2::add_wfc_Gau;

  // Harmonic oscillator wavepacket
  void (Wfcgrid2::*expt_add_wfc_HO_v1)
  (vector<double>& x0, vector<double>& px0, vector<double>& alpha, int init_state, vector<int>& nu, complex<double> weight, int rep) = &Wfcgrid2::add_wfc_HO;

  // Arbitrary wavefunction
  void (Wfcgrid2::*expt_add_wfc_ARB_v1)
  (bp::object py_funct, bp::object params, int rep) = &Wfcgrid2::add_wfc_ARB;


  MATRIX (Wfcgrid2::*expt_get_pops_v1)
  (int rep) = &Wfcgrid2::get_pops;

  MATRIX (Wfcgrid2::*expt_get_pops_v2)
  (int rep, vector<double>& bmin, vector<double>& bmax) = &Wfcgrid2::get_pops;



  // Auxiliary functions
  int (*expt_points_on_same_line_v1)
  (int idof, vector<int>& pt1, vector<int>& pt2) = &points_on_same_line;

  def("points_on_same_line", expt_points_on_same_line_v1);


  void (*expt_dvr0_v1)
  (MATRIX& T, int npts, double rmin, double rmax, double mass) = &dvr0;
  def("dvr0", expt_dvr0_v1);

  void (*expt_dvr1_v1)
  (MATRIX& T, int npts, double dr, double mass) = &dvr1;
  def("dvr1", expt_dvr1_v1);

  void (*expt_dvr2_v1)
  (MATRIX& T, int npts, double dr, double mass) = &dvr2;
  def("dvr2", expt_dvr2_v1);

/*
  vector<CMATRIX> (*expt_expV_v1)
  (vector<CMATRIX>& V, complex<double> dt) = &expV;
  def("expV", expt_expV_v1);

  vector<CMATRIX> (*expt_expV_diag_v1)
  (vector<CMATRIX>& V, complex<double> dt) = &expV_diag;
  def("expV_diag", expt_expV_diag_v1);
*/


  class_<Wfcgrid2>("Wfcgrid2",init<vector<double>&, vector<double>&, vector<double>&, int>())
      .def(init<const Wfcgrid2&>())
      .def("__copy__", &generic__copy__<Wfcgrid2>)
      .def("__deepcopy__", &generic__deepcopy__<Wfcgrid2>)

      .def_readwrite("nstates", &Wfcgrid2::nstates)
      .def_readwrite("ndof", &Wfcgrid2::ndof)
      .def_readwrite("Npts", &Wfcgrid2::Npts)
      .def_readwrite("npts", &Wfcgrid2::npts)
      .def_readwrite("rmin", &Wfcgrid2::rmin)
      .def_readwrite("rmax", &Wfcgrid2::rmax)
      .def_readwrite("dr",   &Wfcgrid2::dr)
      .def_readwrite("kmin", &Wfcgrid2::kmin)
      .def_readwrite("dk",   &Wfcgrid2::dk)

      .def_readwrite("gmap",   &Wfcgrid2::gmap)

      .def_readwrite("PSI_dia", &Wfcgrid2::PSI_dia)
      .def_readwrite("reciPSI_dia", &Wfcgrid2::reciPSI_dia)
      .def_readwrite("PSI_adi", &Wfcgrid2::PSI_adi)
      .def_readwrite("reciPSI_adi", &Wfcgrid2::reciPSI_adi)


      /*
      .def_readwrite("DtreciPSI", &Wfcgrid::DtreciPSI)
      .def_readwrite("DxPSI", &Wfcgrid::DxPSI)
      .def_readwrite("DyPSI", &Wfcgrid::DyPSI)
      .def_readwrite("KxreciPSI", &Wfcgrid::KxreciPSI)
      .def_readwrite("KyreciPSI", &Wfcgrid::KyreciPSI)
      */

      .def_readwrite("Hdia", &Wfcgrid2::Hdia)
      .def_readwrite("Hadi", &Wfcgrid2::Hadi)
      .def_readwrite("Vcomplex", &Wfcgrid2::Vcomplex)
      .def_readwrite("NAC1", &Wfcgrid2::NAC1)
      .def_readwrite("NAC2", &Wfcgrid2::NAC2)
      .def_readwrite("U",    &Wfcgrid2::U)
      .def_readwrite("expH", &Wfcgrid2::expH)
      .def_readwrite("expK", &Wfcgrid2::expK)


      .def("imap", &Wfcgrid2::imap)

      /**  Wfcgrid2_ColbertMiller    */
      .def("T_PSI", &Wfcgrid2::T_PSI)
      .def("T_PSI_adi", &Wfcgrid2::T_PSI_adi)
      .def("T_PSI_dia", &Wfcgrid2::T_PSI_dia)
      .def("operator_T", &Wfcgrid2::operator_T)
      .def("expT_PSI", &Wfcgrid2::expT_PSI)
      .def("Colbert_Miller_SOFT", &Wfcgrid2::Colbert_Miller_SOFT)

      .def("nubla_PSI", &Wfcgrid2::nubla_PSI)
      .def("nubla_PSI_adi", &Wfcgrid2::nubla_PSI_adi)
      .def("nubla_PSI_dia", &Wfcgrid2::nubla_PSI_dia)

      .def("V_PSI", &Wfcgrid2::V_PSI)
      .def("V_PSI_adi", &Wfcgrid2::V_PSI_adi)
      .def("V_PSI_dia", &Wfcgrid2::V_PSI_dia)

      .def("H_PSI_adi", &Wfcgrid2::H_PSI_adi)
      .def("H_PSI_dia", &Wfcgrid2::H_PSI_dia)

      .def("Colbert_Miller_propagate_adi1", &Wfcgrid2::Colbert_Miller_propagate_adi1)
      .def("Colbert_Miller_propagate_adi2", &Wfcgrid2::Colbert_Miller_propagate_adi2)
      .def("Colbert_Miller_propagate_dia1", &Wfcgrid2::Colbert_Miller_propagate_dia1)
      .def("Colbert_Miller_propagate_dia2", &Wfcgrid2::Colbert_Miller_propagate_dia2)



      /**  Wfcgrid2_direct    */
      .def("direct_allocate_tmp_vars", &Wfcgrid2::direct_allocate_tmp_vars)
      .def("direct_propagate_adi2", &Wfcgrid2::direct_propagate_adi2)
      .def("direct_propagate_adi1", &Wfcgrid2::direct_propagate_adi1)
      .def("direct_propagate_dia2", &Wfcgrid2::direct_propagate_dia2)
      .def("direct_propagate_dia1", &Wfcgrid2::direct_propagate_dia1)


      /**  Wfcgrid2_initialize    */
      .def("add_wfc_Gau", expt_add_wfc_Gau_v1)
      .def("add_wfc_HO", expt_add_wfc_HO_v1)
      .def("add_wfc_ARB", expt_add_wfc_ARB_v1)

      /**  Wfcgrid2_properties    */
      .def("norm", &Wfcgrid2::norm)
      .def("e_kin", &Wfcgrid2::e_kin)
      .def("e_pot", &Wfcgrid2::e_pot)
      .def("e_tot", &Wfcgrid2::e_tot)
      .def("get_pow_q", &Wfcgrid2::get_pow_q)
      .def("get_pow_p", &Wfcgrid2::get_pow_p)
      .def("get_den_mat", &Wfcgrid2::get_den_mat)
      .def("get_pops", expt_get_pops_v1)
      .def("get_pops", expt_get_pops_v2)

      /**  Wfcgrid2_SOFT    */
      .def("update_propagator_H", &Wfcgrid2::update_propagator_H)
      .def("update_propagator_K", &Wfcgrid2::update_propagator_K)
      .def("SOFT_propagate", &Wfcgrid2::SOFT_propagate)

      /**  Wfcgrid2_transforms    */
      .def("update_reciprocal", &Wfcgrid2::update_reciprocal)
      .def("update_real", &Wfcgrid2::update_real)
      .def("normalize", &Wfcgrid2::normalize)
      .def("reshape_wfc_1D", &Wfcgrid2::reshape_wfc_1D)
      .def("reshape_wfc_2D", &Wfcgrid2::reshape_wfc_2D)

      /**  Wfcgrid2_updates    */
      .def("update_Hamiltonian", &Wfcgrid2::update_Hamiltonian)
      .def("update_adiabatic", &Wfcgrid2::update_adiabatic)
      .def("update_diabatic", &Wfcgrid2::update_diabatic)

      /**  Wfcgrid2_io    */
      .def("print_wfc_1D", &Wfcgrid2::print_wfc_1D)
      .def("print_reci_wfc_1D", &Wfcgrid2::print_reci_wfc_1D)

      .def("print_wfc_2D", &Wfcgrid2::print_wfc_2D)
      .def("print_reci_wfc_2D", &Wfcgrid2::print_reci_wfc_2D)


       /*
      .def("print_wfc_1D", &Wfcgrid::print_wfc_1D)
      .def("print_wfc_2D", &Wfcgrid::print_wfc_2D)
      .def("print_reci_wfc_1D", &Wfcgrid::print_reci_wfc_1D)
      .def("print_reci_wfc_2D", &Wfcgrid::print_reci_wfc_2D)
      .def("print_ham_1D", &Wfcgrid::print_ham_1D)
      .def("print_expH_1D", &Wfcgrid::print_expH_1D)
      .def("print_expK_1D", &Wfcgrid::print_expK_1D)

      .def("print_populations_1D", &Wfcgrid::print_populations_1D)
      .def("print_populations_2D", &Wfcgrid::print_populations_2D)

      .def("absorb_1D",expt_absorb_1D)

      .def("flux_1D", &Wfcgrid::flux_1D)

      .def("get_x_1D", &Wfcgrid::get_x_1D)
      .def("get_pow_x_1D", &Wfcgrid::get_pow_x_1D)
      .def("get_px_1D", &Wfcgrid::get_px_1D)
      .def("get_pow_px_1D", &Wfcgrid::get_pow_px_1D)

      .def("get_x_2D", &Wfcgrid::get_x_2D)
      .def("get_y_2D", &Wfcgrid::get_y_2D)
      .def("get_px_2D", &Wfcgrid::get_px_2D)
      .def("get_py_2D", &Wfcgrid::get_py_2D)

      */
 
  ;



}


#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygwfcgrid2){
#else
BOOST_PYTHON_MODULE(libwfcgrid2){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

//  export_Mathematics_objects();
  export_Wfcgrid2_objects();

}


}// namespace libwfcgrid2
}// namespace libdyn
}// liblibra
