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
#include "libdyn.h"
using namespace boost::python;


namespace libdyn{

using namespace libnuclear;
using namespace libelectronic;
using namespace librigidbody;
using namespace libthermostat;
using namespace libbarostat;
using namespace libwfcgrid;
using namespace libensemble;


void export_Dyn_objects(){

  export_Nuclear_objects();
  export_RigidBody_objects();
  export_Electronic_objects();
  export_Thermostat_objects();
  export_Barostat_objects();
  export_Wfcgrid_objects();
  export_Ensemble_objects();

  double (*expt_compute_kinetic_energy_v1)(Nuclear& mol) = &compute_kinetic_energy;
  double (*expt_compute_kinetic_energy_v2)(Ensemble& ens) = &compute_kinetic_energy;
  double (*expt_compute_potential_energy_v1)(Nuclear& mol, Electronic& el, Hamiltonian& ham, int opt) = &compute_potential_energy;
  double (*expt_compute_potential_energy_v2)(Ensemble& ens, int opt) = &compute_potential_energy;

  void (*expt_compute_forces_v1)(Nuclear& mol, Electronic& el, Hamiltonian& ham, int opt) = &compute_forces;
  void (*expt_compute_forces_v2)(Ensemble& ens, int opt) = &compute_forces;


  void (*expt_compute_hopping_probabilities_fssh_v1)
  (Nuclear& mol, Electronic& el, Hamiltonian& ham, MATRIX& g,
   double dt, int use_boltz_factor,double T) = &compute_hopping_probabilities_fssh;
  void (*expt_compute_hopping_probabilities_fssh_v2)
  (Ensemble& ens, int i, MATRIX& g,
   double dt, int use_boltz_factor,double T) = &compute_hopping_probabilities_fssh;


  void (*expt_compute_hopping_probabilities_gfsh_v1)
  (Nuclear& mol, Electronic& el, Hamiltonian& ham, MATRIX& g,
   double dt, int use_boltz_factor,double T) = &compute_hopping_probabilities_gfsh;
  void (*expt_compute_hopping_probabilities_gfsh_v2)
  (Ensemble& ens, int i, MATRIX& g,
   double dt, int use_boltz_factor,double T) = &compute_hopping_probabilities_gfsh;


  void (*expt_compute_hopping_probabilities_mssh_v1)
  (Nuclear& mol, Electronic& el, Hamiltonian& ham, MATRIX& g,
   double dt, int use_boltz_factor,double T) = &compute_hopping_probabilities_mssh;
  void (*expt_compute_hopping_probabilities_mssh_v2)
  (Ensemble& ens, int i, MATRIX& g,
   double dt, int use_boltz_factor,double T) = &compute_hopping_probabilities_mssh;


  void (*expt_compute_hopping_probabilities_esh_v1)
  (Ensemble& ens, MATRIX* g, double dt, int use_boltz_factor,double T) = compute_hopping_probabilities_esh;



  int (*expt_hop_v1)
  (int initstate, Nuclear& mol, Hamiltonian& ham, double ksi, MATRIX& g, int do_rescaling, int rep, int do_reverse) = &hop;
  int (*expt_hop_v2)
  (int initstate, Ensemble& ens, int i, double ksi, MATRIX& g, int do_rescaling, int rep, int do_reverse) = &hop;


  int (*expt_rescale_velocities_adiabatic_v1)
  (Nuclear& mol, Hamiltonian& ham, int old_st, int do_reverse) = &rescale_velocities_adiabatic;

  int (*expt_rescale_velocities_diabatic_v1)
  (Nuclear& mol, Hamiltonian& ham, int old_st) = &rescale_velocities_diabatic;




  def("compute_kinetic_energy",expt_compute_kinetic_energy_v1);
  def("compute_kinetic_energy",expt_compute_kinetic_energy_v2);
  def("compute_potential_energy",expt_compute_potential_energy_v1);
  def("compute_potential_energy",expt_compute_potential_energy_v2);
  def("compute_forces",expt_compute_forces_v1);
  def("compute_forces",expt_compute_forces_v2);

  def("compute_hopping_probabilities_fssh",expt_compute_hopping_probabilities_fssh_v1);
  def("compute_hopping_probabilities_fssh",expt_compute_hopping_probabilities_fssh_v2);
  def("compute_hopping_probabilities_gfsh",expt_compute_hopping_probabilities_gfsh_v1);
  def("compute_hopping_probabilities_gfsh",expt_compute_hopping_probabilities_gfsh_v2);
  def("compute_hopping_probabilities_mssh",expt_compute_hopping_probabilities_mssh_v1);
  def("compute_hopping_probabilities_mssh",expt_compute_hopping_probabilities_mssh_v2);
  def("compute_hopping_probabilities_esh", expt_compute_hopping_probabilities_esh_v1);

  def("hop", expt_hop_v1);
  def("hop", expt_hop_v2);

  def("rescale_velocities_adiabatic", expt_rescale_velocities_adiabatic_v1);
  def("rescale_velocities_diabatic", expt_rescale_velocities_diabatic_v1);



}// export_Dyn_objects()




#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygdyn){
#else
BOOST_PYTHON_MODULE(libdyn){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

  export_Dyn_objects();

}


}// libdyn

