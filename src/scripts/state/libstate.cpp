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

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "libstate.h"

/// liblibra namespace
namespace liblibra{

using namespace boost::python;


namespace libscripts{
namespace libstate{


void export_state_objects(){

  class_<MD>("MD",init<>())
      .def(init<boost::python::dict>())
      .def_readwrite("ensemble",&MD::ensemble)
      .def_readwrite("integrator",&MD::integrator)
      .def_readwrite("dt",&MD::dt)
      .def_readwrite("n_medium",&MD::n_medium)
      .def_readwrite("n_fast",&MD::n_fast)
      .def_readwrite("n_outer",&MD::n_outer)
      .def_readwrite("max_time",&MD::max_time)
      .def_readwrite("max_step",&MD::max_step)
      .def_readwrite("curr_time",&MD::curr_time)
      .def_readwrite("curr_step",&MD::curr_step)
      .def_readwrite("terec_exp_size",&MD::terec_exp_size)
      .def_readwrite("use_vlist",&MD::use_vlist)
      .def_readwrite("vlist_upd_freq",&MD::vlist_upd_freq)
      .def_readwrite("vlist_time",&MD::vlist_time)

      .def("show_info",&MD::show_info)
  ;


  void (State::*expt_init_md_v1)(Nuclear& mol, Electronic& el, Hamiltonian& ham, Random& rnd) = &State::init_md;
  void (State::*expt_init_md_v2)(              Electronic& el, Hamiltonian& ham, Random& rnd) = &State::init_md;

  void (State::*expt_run_md_v1)(Nuclear& mol, Electronic& el, Hamiltonian& ham) = &State::run_md;
  void (State::*expt_run_md_v2)(Electronic& el, Hamiltonian& ham) = &State::run_md;


  class_<State>("State",init<>())
      .def_readwrite("E_kin",&State::E_kin)
      .def_readwrite("E_kin_tr",&State::E_kin_tr)
      .def_readwrite("E_kin_rot",&State::E_kin_rot)
      .def_readwrite("E_pot",&State::E_pot)
      .def_readwrite("E_tot",&State::E_tot)
      .def_readwrite("H",&State::H)
      .def_readwrite("curr_T",&State::curr_T)
      .def_readwrite("curr_P",&State::curr_P)
      .def_readwrite("curr_V",&State::curr_V)
      .def_readwrite("H_NP",&State::H_NP)
      .def_readwrite("L_tot",&State::L_tot)
      .def_readwrite("P_tot",&State::P_tot)

      .def("set_system",&State::set_system)
      .def("set_thermostat",&State::set_thermostat)
      .def("set_barostat",&State::set_barostat)
      .def("set_md",&State::set_md)
      .def("show_info",&State::show_info)

       // Defined in State_methods.cpp
      .def("update",&State::update)
      .def("cool",&State::cool)

      .def("init_md",expt_init_md_v1)
      .def("init_md",expt_init_md_v2)

      .def("run_md",expt_run_md_v1)
      .def("run_md",expt_run_md_v2)
//      .def("run_md_respa",&State::run_md_respa)
  ;



}


#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygstate){
#else
BOOST_PYTHON_MODULE(libstate){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

  export_state_objects();

}



}// namespace libstate
}// namespace libscripts
}// liblibra
