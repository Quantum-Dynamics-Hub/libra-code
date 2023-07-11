/*********************************************************************************
* Copyright (C) 2015-2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file libcontrol_parameters.cpp
  \brief The file implements Python export function
    
*/

#if defined(USING_PCH)
#include "../pch.h"
#else
#include <memory> // for std::auto_ptr<>
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#endif 

#include "libcontrol_parameters.h"
#include "../math_linalg/liblinalg.h"


/// liblibra namespace
namespace liblibra{

using namespace liblinalg;
using namespace boost::python;


/// libcontrol_parameters namespace
namespace libcontrol_parameters{



void export_Control_Parameters_objects(){
/** 
  \brief Exporter of the libcontrol_parameters classes and functions

*/



  class_<Control_Parameters>("Control_Parameters",init<>())
      .def("__copy__", &generic__copy__<Control_Parameters>)
      .def("__deepcopy__", &generic__deepcopy__<Control_Parameters>)

      .def_readwrite("runtype", &Control_Parameters::runtype)
      .def_readwrite("hamiltonian", &Control_Parameters::hamiltonian)
      .def_readwrite("spin_method", &Control_Parameters::spin_method)
      .def_readwrite("DF", &Control_Parameters::DF)

      .def_readwrite("guess_type", &Control_Parameters::guess_type)

      .def_readwrite("scf_algo", &Control_Parameters::scf_algo)
      .def_readwrite("use_disk", &Control_Parameters::use_disk)
      .def_readwrite("use_rosh", &Control_Parameters::use_rosh)
      .def_readwrite("do_annihilate", &Control_Parameters::do_annihilate)
      .def_readwrite("pop_opt", &Control_Parameters::pop_opt)
      .def_readwrite("use_diis", &Control_Parameters::use_diis)
      .def_readwrite("diis_max", &Control_Parameters::diis_max)
      .def_readwrite("diis_start_iter", &Control_Parameters::diis_start_iter)
      .def_readwrite("use_level_shift", &Control_Parameters::use_level_shift)
      .def_readwrite("shift_magnitude", &Control_Parameters::shift_magnitude)
      .def_readwrite("use_damping", &Control_Parameters::use_damping)
      .def_readwrite("damping_start", &Control_Parameters::damping_start)
      .def_readwrite("damping_const", &Control_Parameters::damping_const)
      .def_readwrite("etol", &Control_Parameters::etol)
      .def_readwrite("den_tol", &Control_Parameters::den_tol)
      .def_readwrite("Niter", &Control_Parameters::Niter)
      .def_readwrite("degen_tol", &Control_Parameters::degen_tol)

      .def_readwrite("parameters", &Control_Parameters::parameters)
      .def_readwrite("eht_params_format", &Control_Parameters::eht_params_format)
      .def_readwrite("eht_formula", &Control_Parameters::eht_formula)
      .def_readwrite("eht_sce_formula", &Control_Parameters::eht_sce_formula)
      .def_readwrite("eht_fock_opt", &Control_Parameters::eht_fock_opt)
      .def_readwrite("eht_electrostatics", &Control_Parameters::eht_electrostatics)


      .def_readwrite("compute_vertical_ip", &Control_Parameters::compute_vertical_ip)
      .def_readwrite("compute_vertical_ea", &Control_Parameters::compute_vertical_ea)

      .def_readwrite("md_dt", &Control_Parameters::md_dt)
      .def_readwrite("md_nsteps", &Control_Parameters::md_nsteps)

      .def_readwrite("opt_dt", &Control_Parameters::opt_dt)
      .def_readwrite("opt_nsteps", &Control_Parameters::opt_nsteps)

      .def_readwrite("compute_dipole", &Control_Parameters::compute_dipole)

      .def_readwrite("compute_dos", &Control_Parameters::compute_dos)
      .def_readwrite("dos_opt", &Control_Parameters::dos_opt)
      .def_readwrite("dos_prefix", &Control_Parameters::dos_prefix)

      .def_readwrite("compute_charge_density", &Control_Parameters::compute_charge_density)
      .def_readwrite("nx_grid", &Control_Parameters::nx_grid)
      .def_readwrite("ny_grid", &Control_Parameters::ny_grid)
      .def_readwrite("nz_grid", &Control_Parameters::nz_grid)
      .def_readwrite("charge_density_prefix", &Control_Parameters::charge_density_prefix)
      .def_readwrite("orbs", &Control_Parameters::orbs)

      .def_readwrite("nac_md_trajectory_filename", &Control_Parameters::nac_md_trajectory_filename)
      .def_readwrite("nac_prefix", &Control_Parameters::nac_prefix)
      .def_readwrite("nac_min_frame", &Control_Parameters::nac_min_frame)
      .def_readwrite("nac_max_frame", &Control_Parameters::nac_max_frame)
      .def_readwrite("nac_min_orbs", &Control_Parameters::nac_min_orbs)
      .def_readwrite("nac_max_orbs", &Control_Parameters::nac_max_orbs)
      .def_readwrite("nac_dt", &Control_Parameters::nac_dt)
      .def_readwrite("nac_opt", &Control_Parameters::nac_opt)

      .def_readwrite("scan_mov_at", &Control_Parameters::scan_mov_at)
      .def_readwrite("scan_ref_at", &Control_Parameters::scan_ref_at)
      .def_readwrite("scan_dir", &Control_Parameters::scan_dir)
      .def_readwrite("scan_dxmin", &Control_Parameters::scan_dxmin)
      .def_readwrite("scan_dxmax", &Control_Parameters::scan_dxmax)
      .def_readwrite("scan_dx", &Control_Parameters::scan_dx)


      .def_readwrite("compute_excitations", &Control_Parameters::compute_excitations)
      .def_readwrite("num_excitations", &Control_Parameters::num_excitations)
      .def_readwrite("excitations_opt", &Control_Parameters::excitations_opt)
      .def_readwrite("spectral_width", &Control_Parameters::spectral_width)
      .def_readwrite("excitations", &Control_Parameters::excitations)  

      .def_readwrite("t1", &Control_Parameters::t1)
      .def_readwrite("t2", &Control_Parameters::t2)
      .def_readwrite("t3", &Control_Parameters::t3)
      .def_readwrite("x_period", &Control_Parameters::x_period)
      .def_readwrite("y_period", &Control_Parameters::y_period)
      .def_readwrite("z_period", &Control_Parameters::z_period)

      .def_readwrite("Natoms", &Control_Parameters::Natoms)
      .def_readwrite("charge", &Control_Parameters::charge)
      .def_readwrite("spin", &Control_Parameters::spin)
      .def_readwrite("coordinates", &Control_Parameters::coordinates)

      .def_readwrite("fragments", &Control_Parameters::fragments)
      .def_readwrite("frag_size", &Control_Parameters::frag_size)
      .def_readwrite("frag_name", &Control_Parameters::frag_name)
      .def_readwrite("frag_charge", &Control_Parameters::frag_charge)

  ;

  void (*expt_get_parameters_from_file_v1)
  (std::string, Control_Parameters&) = &get_parameters_from_file;


  def("get_parameters_from_file",expt_get_parameters_from_file_v1);


}


#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygcontrol_parameters){
#else
BOOST_PYTHON_MODULE(libcontrol_parameters){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

  export_Control_Parameters_objects();

}



}// namespace libcontrol_parameters
}// namespace liblibra


