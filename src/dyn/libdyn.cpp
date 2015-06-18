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


void export_Dyn_objects(){

  export_Nuclear_objects();
  export_RigidBody_objects();
  export_Electronic_objects();
  export_Thermostat_objects();
  export_Barostat_objects();
  export_Wfcgrid_objects();

  double (*expt_compute_kinetic_energy_v1)(Nuclear& mol) = &compute_kinetic_energy;
  double (*expt_compute_potential_energy_v1)(Nuclear& mol, Electronic& el, Hamiltonian& ham, int opt) = &compute_potential_energy;
  void (*expt_compute_forces_v1)(Nuclear& mol, Electronic& el, Hamiltonian& ham, int opt) = &compute_forces;

  void (*expt_compute_hopping_probabilities_fssh_v1)
  (Nuclear& mol, Electronic& el, Hamiltonian& ham, MATRIX& g,
   double dt, int use_boltz_factor,double T) = &compute_hopping_probabilities_fssh;

  void (*expt_compute_hopping_probabilities_gfsh_v1)
  (Nuclear& mol, Electronic& el, Hamiltonian& ham, MATRIX& g,
   double dt, int use_boltz_factor,double T) = &compute_hopping_probabilities_gfsh;

  void (*expt_compute_hopping_probabilities_mssh_v1)
  (Nuclear& mol, Electronic& el, Hamiltonian& ham, MATRIX& g,
   double dt, int use_boltz_factor,double T) = &compute_hopping_probabilities_mssh;


  int (*expt_hop_v1)
  (int initstate, Nuclear& mol, Hamiltonian& ham, double ksi, MATRIX& g, int do_rescaling, int rep) = &hop;



  def("compute_kinetic_energy",expt_compute_kinetic_energy_v1);
  def("compute_potential_energy",expt_compute_potential_energy_v1);
  def("compute_forces",expt_compute_forces_v1);

  def("compute_hopping_probabilities_fssh",expt_compute_hopping_probabilities_fssh_v1);
  def("compute_hopping_probabilities_gfsh",expt_compute_hopping_probabilities_gfsh_v1);
  def("compute_hopping_probabilities_mssh",expt_compute_hopping_probabilities_mssh_v1);

  def("hop", expt_hop_v1);

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

