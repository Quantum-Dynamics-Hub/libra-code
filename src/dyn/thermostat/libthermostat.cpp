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

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "libthermostat.h"

using namespace boost::python;


namespace libdyn{
namespace libthermostat{



void export_Thermostat_objects(){

  void (Thermostat::*expt_set_Nf_t_v1)(int nf_t) = &Thermostat::set_Nf_t;
  void (Thermostat::*expt_set_Nf_t_v2)(double nf_t) = &Thermostat::set_Nf_t;
  void (Thermostat::*expt_set_Nf_r_v1)(int nf_t) = &Thermostat::set_Nf_r;
  void (Thermostat::*expt_set_Nf_r_v2)(double nf_t) = &Thermostat::set_Nf_r;
  void (Thermostat::*expt_set_Nf_b_v1)(int nf_t) = &Thermostat::set_Nf_b;
  void (Thermostat::*expt_set_Nf_b_v2)(double nf_t) = &Thermostat::set_Nf_b;

  void (Thermostat::*expt_update_thermostat_forces_v1)(double, double, double) = &Thermostat::update_thermostat_forces;
  void (Thermostat::*expt_update_thermostat_forces_v2)(double, double, double, int) = &Thermostat::update_thermostat_forces;



  class_<Thermostat>("Thermostat",init<>())
      .def(init<boost::python::dict>())
      .def_readwrite("s_var",&Thermostat::s_var)
      .def_readwrite("Ps",&Thermostat::Ps)
      .def_readwrite("Q",&Thermostat::Q)
      .def_readwrite("NHC_size",&Thermostat::NHC_size)
      .def_readwrite("nu_therm",&Thermostat::nu_therm)
      .def_readwrite("Temperature",&Thermostat::Temperature)
      .def_readwrite("thermostat_type",&Thermostat::thermostat_type)

      .def("show_info",&Thermostat::show_info)

      .def("energy", &Thermostat::energy)
      .def("set_Nf_t", expt_set_Nf_t_v1)
      .def("set_Nf_t", expt_set_Nf_t_v2)
      .def("set_Nf_r", expt_set_Nf_r_v1)
      .def("set_Nf_r", expt_set_Nf_r_v2)
      .def("set_Nf_b", expt_set_Nf_b_v1)
      .def("set_Nf_b", expt_set_Nf_b_v2)
      .def("get_Nf_t", &Thermostat::get_Nf_t)
      .def("get_Nf_r", &Thermostat::get_Nf_r)
      .def("get_Nf_b", &Thermostat::get_Nf_b)

      .def("get_s_var", &Thermostat::get_s_var)
      .def("get_ksi_t", &Thermostat::get_ksi_t)
      .def("get_ksi_r", &Thermostat::get_ksi_r)
      .def("get_ksi_b", &Thermostat::get_ksi_b)

      .def("propagate_sPs", &Thermostat::propagate_sPs)
      .def("propagate_Ps",&Thermostat::propagate_Ps)
      .def("vel_scale",&Thermostat::vel_scale)
      .def("ang_vel_scale",&Thermostat::ang_vel_scale)

      .def("update_thermostat_forces", expt_update_thermostat_forces_v1)
      .def("update_thermostat_forces", expt_update_thermostat_forces_v2)

      .def("init_nhc", &Thermostat::init_nhc)
      .def("propagate_nhc", &Thermostat::propagate_nhc)
      .def("cool", &Thermostat::cool)

  ;

}// export_Thermostat_objects



#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygthermostat){
#else
BOOST_PYTHON_MODULE(libthermostat){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

  export_Thermostat_objects();

}


}// namespace libthermostat
}// namespace libdyn

