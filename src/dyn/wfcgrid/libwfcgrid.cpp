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
  \file libwfcgrid.cpp
  \brief The file implements Python export function
    
*/

#if defined(USING_PCH)
#include "../../pch.h"
#else
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#endif 

#include "libwfcgrid.h"

/// liblibra namespace
namespace liblibra{

using namespace boost::python;

/// libdyn namespace
namespace libdyn{

/// libwfcgrid namespace
namespace libwfcgrid{


void export_Wfcgrid_objects(){
/** 
  \brief Exporter of libwfcgrid classes and functions

*/

  int (*expt_compute_imapping_v1)(vector<int>& inp, vector<int>& npts) = compute_imapping;
  def("compute_imapping", expt_compute_imapping_v1);
 
  vector<int> (*expt_compute_mapping_v1)(int indx, vector<int>& npts) = &compute_mapping;
  def("compute_mapping", expt_compute_mapping_v1);

  vector<vector<int> > (*expt_compute_mapping_v2)(vector<vector<int> >& inp, vector<int>& npts) = &compute_mapping;
  def("compute_mapping", expt_compute_mapping_v2);

  vector<int> (*expt_compute_hyperplane_v1)(vector<int>& npts, int idim_const, int ipt_const) = &compute_hyperplane;
  def("compute_hyperplane", expt_compute_hyperplane_v1);


  boost::python::list (Wfcgrid::*expt_absorb_1D)(double dL) = &Wfcgrid::absorb_1D;

  void (Wfcgrid::*expt_update_potential_1D_v1)(Hamiltonian& ham) = &Wfcgrid::update_potential_1D;
  void (Wfcgrid::*expt_update_potential_1D_v2)(bp::object py_funct, bp::object params) = &Wfcgrid::update_potential_1D;
  void (Wfcgrid::*expt_update_potential_2D_v1)(Hamiltonian& ham) = &Wfcgrid::update_potential_2D;
  void (Wfcgrid::*expt_update_potential_2D_v2)(bp::object py_funct, bp::object params) = &Wfcgrid::update_potential_2D;


  void (Wfcgrid::*expt_update_propagator_K_2D_v1)(double dt, double m1, double m2) = &Wfcgrid::update_propagator_K_2D;
  void (Wfcgrid::*expt_update_propagator_K_2D_v2)(double dt, double m0) = &Wfcgrid::update_propagator_K_2D;


  // Gaussian wavepackets
  void (Wfcgrid::*expt_init_wfc_1D_v1)(double x0, double px0, double dx, int init_state) = &Wfcgrid::init_wfc_1D;
  void (Wfcgrid::*expt_init_wfc_1D_v2)
  (vector<double>& x0, vector<double>& px0, vector<double>& dx, vector<int>& init_state, vector<complex<double> >& weights) = &Wfcgrid::init_wfc_1D;

  // Harmonic oscillator
  void (Wfcgrid::*expt_init_wfc_1D_HO_v1)
  (vector<int>& init_state, vector<int>& nu, vector<complex<double> >& weights,
   vector<double>& x0, vector<double>& px0, vector<double>& alpha) = &Wfcgrid::init_wfc_1D_HO;

  // Arbitrary wavefunction
  void (Wfcgrid::*expt_init_wfc_1D_ARB_v1)
  (bp::object py_funct, bp::object params) = &Wfcgrid::init_wfc_1D_ARB;


  class_<Wfcgrid>("Wfcgrid",init<double,double,double,int>())
      .def(init<double, double, double, double, double, double, int>())
      .def(init<const Wfcgrid&>())
      .def("__copy__", &generic__copy__<Wfcgrid>)
      .def("__deepcopy__", &generic__deepcopy__<Wfcgrid>)

      .def_readwrite("nstates", &Wfcgrid::nstates)
      .def_readwrite("Nx", &Wfcgrid::Nx)
      .def_readwrite("Ny", &Wfcgrid::Ny)
      .def_readwrite("xmin", &Wfcgrid::xmin)
      .def_readwrite("ymin", &Wfcgrid::ymin)
      .def_readwrite("xmax", &Wfcgrid::xmax)
      .def_readwrite("ymax", &Wfcgrid::ymax)
      .def_readwrite("dx", &Wfcgrid::dx)
      .def_readwrite("dy", &Wfcgrid::dy)
      .def_readwrite("kxmin", &Wfcgrid::kxmin)
      .def_readwrite("kymin", &Wfcgrid::kymin)

      .def_readwrite("PSI", &Wfcgrid::PSI)
      .def_readwrite("reciPSI", &Wfcgrid::reciPSI)
      .def_readwrite("DtreciPSI", &Wfcgrid::DtreciPSI)
      .def_readwrite("DxPSI", &Wfcgrid::DxPSI)
      .def_readwrite("DyPSI", &Wfcgrid::DyPSI)
      .def_readwrite("KxreciPSI", &Wfcgrid::KxreciPSI)
      .def_readwrite("KyreciPSI", &Wfcgrid::KyreciPSI)

      .def_readwrite("H", &Wfcgrid::H)
      .def_readwrite("Dx", &Wfcgrid::Dx)
      .def_readwrite("Dy", &Wfcgrid::Dy)
      .def_readwrite("expH", &Wfcgrid::expH)
      .def_readwrite("expK", &Wfcgrid::expK)


      .def("init_wfc_1D", expt_init_wfc_1D_v1)
      .def("init_wfc_1D", expt_init_wfc_1D_v2)
      .def("init_wfc_1D_HO", expt_init_wfc_1D_HO_v1)
      .def("init_wfc_1D_ARB", expt_init_wfc_1D_ARB_v1)

      .def("init_wfc_2D", &Wfcgrid::init_wfc_2D)

      .def("print_wfc_1D", &Wfcgrid::print_wfc_1D)
      .def("print_wfc_2D", &Wfcgrid::print_wfc_2D)
      .def("print_reci_wfc_1D", &Wfcgrid::print_reci_wfc_1D)
      .def("print_reci_wfc_2D", &Wfcgrid::print_reci_wfc_2D)
      .def("print_ham_1D", &Wfcgrid::print_ham_1D)
      .def("print_expH_1D", &Wfcgrid::print_expH_1D)
      .def("print_expK_1D", &Wfcgrid::print_expK_1D)

      .def("print_populations_1D", &Wfcgrid::print_populations_1D)
      .def("print_populations_2D", &Wfcgrid::print_populations_2D)

      .def("update_potential_1D", expt_update_potential_1D_v1)
      .def("update_potential_1D", expt_update_potential_1D_v2)
      .def("update_potential_2D", expt_update_potential_2D_v1)
      .def("update_potential_2D", expt_update_potential_2D_v2)

      .def("update_propagator_1D", &Wfcgrid::update_propagator_1D)
      .def("update_propagator_2D", &Wfcgrid::update_propagator_2D)

      .def("update_propagator_K_1D", &Wfcgrid::update_propagator_K_1D)
      .def("update_propagator_K_2D", expt_update_propagator_K_2D_v1)
      .def("update_propagator_K_2D", expt_update_propagator_K_2D_v2)

      .def("propagate_exact_1D", &Wfcgrid::propagate_exact_1D)
      .def("propagate_exact_2D", &Wfcgrid::propagate_exact_2D)

      .def("normalize_wfc_1D", &Wfcgrid::normalize_wfc_1D)

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

      .def("norm_1D", &Wfcgrid::norm_1D)
      .def("e_kin_1D", &Wfcgrid::e_kin_1D)
      .def("e_pot_1D", &Wfcgrid::e_pot_1D)
      .def("e_tot_1D", &Wfcgrid::e_tot_1D)

      .def("norm_2D", &Wfcgrid::norm_2D)
      .def("e_kin_2D", &Wfcgrid::e_kin_2D)
      .def("e_pot_2D", &Wfcgrid::e_pot_2D)
      .def("e_tot_2D", &Wfcgrid::e_tot_2D)

 
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
}// liblibra
