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
  \file libdyn.cpp
  \brief The file implements Python export function
    
*/

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "libdyn.h"


/// liblibra namespace
namespace liblibra{

using namespace boost::python;

/// libdyn namespace
namespace libdyn{

using namespace libnuclear;
using namespace libelectronic;
//using namespace librigidbody;
using namespace libthermostat;
using namespace libbarostat;
using namespace libwfcgrid;
using namespace libensemble;
using namespace libgwp;


void export_Dyn_objects(){
/** 
  \brief Exporter of libdyn classes and functions

*/


  export_Nuclear_objects();
//  export_RigidBody_objects();
  export_Electronic_objects();
  export_Thermostat_objects();
  export_Barostat_objects();
  export_Wfcgrid_objects();
  export_Ensemble_objects();
  export_gwp_objects();


  double (*expt_Ehrenfest_dia_v1)(CMATRIX& C, CMATRIX& H, vector<CMATRIX>& dHdR, vector<double>& f, int opt) = &Ehrenfest_dia;
  double (*expt_Ehrenfest_adi_v1)(CMATRIX& C, CMATRIX& E, vector<CMATRIX>& dEdR, vector<CMATRIX>& D, vector<double>& f, int opt) = &Ehrenfest_adi;
  def("Ehrenfest_dia",expt_Ehrenfest_dia_v1);
  def("Ehrenfest_adi",expt_Ehrenfest_adi_v1);



  double (*expt_compute_kinetic_energy_v1)(Nuclear& mol) = &compute_kinetic_energy;
  double (*expt_compute_kinetic_energy_v2)(Ensemble& ens) = &compute_kinetic_energy;
  double (*expt_compute_potential_energy_v1)(Nuclear& mol, Electronic& el, Hamiltonian& ham, int opt) = &compute_potential_energy;
  double (*expt_compute_potential_energy_v2)(Ensemble& ens, int opt) = &compute_potential_energy;

  double (*expt_compute_forces_v1)(Nuclear& mol, Electronic& el, Hamiltonian& ham, int opt) = &compute_forces;
  double (*expt_compute_forces_v2)(Ensemble& ens, int opt) = &compute_forces;


  def("compute_kinetic_energy",expt_compute_kinetic_energy_v1);
  def("compute_kinetic_energy",expt_compute_kinetic_energy_v2);
  def("compute_potential_energy",expt_compute_potential_energy_v1);
  def("compute_potential_energy",expt_compute_potential_energy_v2);
  def("compute_forces",expt_compute_forces_v1);
  def("compute_forces",expt_compute_forces_v2);


  void (*expt_propagate_ensemble_v1)(double dt,Ensemble& ens,int opt) = &propagate_ensemble;
  def("propagate_ensemble", expt_propagate_ensemble_v1);



  void (*expt_compute_hopping_probabilities_fssh_v1)
  (Nuclear& mol, Electronic& el, Hamiltonian& ham, MATRIX& g,
   double dt, int use_boltz_factor,double T) = &compute_hopping_probabilities_fssh;
  void (*expt_compute_hopping_probabilities_fssh_v2)
  (Ensemble& ens, int i, MATRIX& g,
   double dt, int use_boltz_factor,double T) = &compute_hopping_probabilities_fssh;
  MATRIX (*expt_compute_hopping_probabilities_fssh_v3)
  (CMATRIX& Coeff, CMATRIX& Hvib, double dt, int use_boltz_factor,double T) = &compute_hopping_probabilities_fssh;


  void (*expt_compute_hopping_probabilities_gfsh_v1)
  (Nuclear& mol, Electronic& el, Hamiltonian& ham, MATRIX& g,
   double dt, int use_boltz_factor,double T) = &compute_hopping_probabilities_gfsh;

  void (*expt_compute_hopping_probabilities_gfsh_v2)
  (Ensemble& ens, int i, MATRIX& g,
   double dt, int use_boltz_factor,double T) = &compute_hopping_probabilities_gfsh;

  MATRIX (*expt_compute_hopping_probabilities_gfsh_v3)
  (CMATRIX& Coeff, CMATRIX& Hvib, double dt, int use_boltz_factor,double T) = &compute_hopping_probabilities_gfsh;


  void (*expt_compute_hopping_probabilities_mssh_v1)
  (Nuclear& mol, Electronic& el, Hamiltonian& ham, MATRIX& g,
   double dt, int use_boltz_factor,double T) = &compute_hopping_probabilities_mssh;
  void (*expt_compute_hopping_probabilities_mssh_v2)
  (Ensemble& ens, int i, MATRIX& g,
   double dt, int use_boltz_factor,double T) = &compute_hopping_probabilities_mssh;
  MATRIX (*expt_compute_hopping_probabilities_mssh_v3)(CMATRIX& Coeff, CMATRIX& Hvib, int use_boltz_factor,double T)  = &compute_hopping_probabilities_mssh;

  void (*expt_compute_hopping_probabilities_esh_v1)
  (Ensemble& ens, MATRIX* g, double dt, int use_boltz_factor,double T) = compute_hopping_probabilities_esh;


  def("compute_hopping_probabilities_fssh",expt_compute_hopping_probabilities_fssh_v1);
  def("compute_hopping_probabilities_fssh",expt_compute_hopping_probabilities_fssh_v2);
  def("compute_hopping_probabilities_fssh",expt_compute_hopping_probabilities_fssh_v3);
  def("compute_hopping_probabilities_gfsh",expt_compute_hopping_probabilities_gfsh_v1);
  def("compute_hopping_probabilities_gfsh",expt_compute_hopping_probabilities_gfsh_v2);
  def("compute_hopping_probabilities_gfsh",expt_compute_hopping_probabilities_gfsh_v3);

  def("compute_hopping_probabilities_mssh",expt_compute_hopping_probabilities_mssh_v1);
  def("compute_hopping_probabilities_mssh",expt_compute_hopping_probabilities_mssh_v2);
  def("compute_hopping_probabilities_mssh",expt_compute_hopping_probabilities_mssh_v3);
  def("compute_hopping_probabilities_esh", expt_compute_hopping_probabilities_esh_v1);


  int (*expt_hop_v1)
  (int initstate, Nuclear& mol, Hamiltonian& ham, double ksi, MATRIX& g, int do_rescaling, int rep, int do_reverse) = &hop;
  int (*expt_hop_v2)
  (int initstate, Ensemble& ens, int i, double ksi, MATRIX& g, int do_rescaling, int rep, int do_reverse) = &hop;
  vector<int> (*expt_hop_v3)
  (int ntraj, vector<int> initstate, vector<Nuclear>& mol, vector<Hamiltonian>& ham, 
   vector<double> ksi, vector<MATRIX>& g, int do_rescaling, int rep, int do_reverse) = &hop;
  int (*expt_hop_v4)(int initstate, MATRIX& g, double ksi) = &hop;

  def("hop", expt_hop_v1);
  def("hop", expt_hop_v2);
  def("hop", expt_hop_v3);
  def("hop", expt_hop_v4);



  int (*expt_rescale_velocities_adiabatic_v1)
  (Nuclear& mol, Hamiltonian& ham, int old_st, int do_reverse) = &rescale_velocities_adiabatic;
  int (*expt_rescale_velocities_diabatic_v1)
  (Nuclear& mol, Hamiltonian& ham, int old_st) = &rescale_velocities_diabatic;

  def("rescale_velocities_adiabatic", expt_rescale_velocities_adiabatic_v1);
  def("rescale_velocities_diabatic", expt_rescale_velocities_diabatic_v1);




  int (*expt_ida_v1)(CMATRIX& Coeff, int old_st, int new_st, double E_old, double E_new, double T, double ksi) = &ida;

  def("ida", expt_ida_v1);


  MATRIX (*expt_coherence_intervals_v1)(const CMATRIX& Coeff, const MATRIX& rates) = &coherence_intervals;

  int (*expt_dish_v1)(Electronic& el, MATRIX& t_m, const MATRIX& tau_m, const CMATRIX& Hvib,
          int use_boltz_flag, double Ekin, double T, double ksi1, double ksi2) = &dish;

  int (*expt_dish_v2)(Electronic& el, Nuclear& mol, Hamiltonian& ham, 
          MATRIX& t_m, const MATRIX& tau_m, int use_boltz_flag, double T, double ksi1, double ksi2) = &dish;

  def("coherence_intervals", expt_coherence_intervals_v1);
  def("dish", expt_dish_v1);
  def("dish", expt_dish_v2);



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
}// liblibra

