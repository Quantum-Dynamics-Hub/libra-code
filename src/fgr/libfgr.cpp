/*********************************************************************************
* Copyright (C) 2019 Xiang Sun, Alexey V. Akimov
* Copyright (C) 2018 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/

#if defined(USING_PCH)
#include "../pch.h"
#else
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#endif

#include "libfgr.h"

//#define BOOST_PYTHON_MAX_ARITY 30


/// liblibra namespace
namespace liblibra{


using namespace boost::python;
using namespace libmeigen;
using namespace librandom;
using namespace liblinalg;

namespace libfgr{


void export_fgr_objects(){


/*
  class_<ivr_params>("ivr_params",init<int>())
    .def(init<const ivr_params&>()) 

    .def("set_qIn",&ivr_params::set_qIn)
    .def("get_invTuningP",&ivr_params::get_invTuningP)

    .def_readwrite("Ndof",&ivr_params::Ndof)
         
  ;


  class_<ivr_observable>("ivr_observables",init<>())
    .def(init<int, int>()) 

    .def_readwrite("observable_type",&ivr_observable::observable_type)
    .def_readwrite("observable_label",&ivr_observable::observable_label)

  ;



*/

  double (*expt_eq_shift_v1)(double Er, double Omega) = &eq_shift;
  double (*expt_reorganization_energy_v1)(double y0, double Omega) = &reorganization_energy;
  double (*expt_reorganization_energy_v2)(vector<double>& omega_nm, vector<double>& req_nm) = &reorganization_energy;
  double (*expt_diabat_crossing_v1)(double dE, double Er, double y0) = &diabat_crossing;
  double (*expt_coupling_Condon_v1)(double gamma, double dE, double Er, double y0) = &coupling_Condon;
  double (*expt_coupling_non_Condon_v1)(double y, double gamma, double dE, double Er, double y0) = &coupling_non_Condon;
  vector<double> (*expt_normal_modes_v1)(vector<double>& omega, vector<double>& coeff, MATRIX& T) = &normal_modes;
  vector<double> (*expt_compute_req_v1)(vector<double>& omega, vector<double>& coeff, double y0, MATRIX& T) = &compute_req;
  vector<double> (*expt_compute_TT_scaled_v1)(MATRIX& T, double scl) = &compute_TT_scaled;
  double (*expt_LVC2GOA_dE_v1)(double E0, double E1, vector<double>& omega_nm, vector<double>& d1, vector<double>& d2) = &LVC2GOA_dE;
  vector<double> (*expt_LVC2GOA_req_v1)(vector<double>& omega_nm, vector<double>& d1, vector<double>& d2) = &LVC2GOA_req;


  def("eq_shift", expt_eq_shift_v1);
  def("reorganization_energy", expt_reorganization_energy_v1);
  def("reorganization_energy", expt_reorganization_energy_v2);
  def("diabat_crossing", expt_diabat_crossing_v1);
  def("coupling_Condon", expt_coupling_Condon_v1);
  def("coupling_non_Condon", expt_coupling_non_Condon_v1);
  def("normal_modes", expt_normal_modes_v1);
  def("compute_req", expt_compute_req_v1);
  def("compute_TT_scaled", expt_compute_TT_scaled_v1);
  def("LVC2GOA_dE", expt_LVC2GOA_dE_v1);
  def("LVC2GOA_req", expt_LVC2GOA_req_v1);



  complex<double> (*expt_Integrand_NE_exact_v1)
  (double tp, double tau, double omega_DA, double omega, double req, double shift, double beta) = &Integrand_NE_exact;
  complex<double> (*expt_Linear_NE_exact_v1)
  (double tp, double tau, double gamma, double omega, double req, double shift, double beta) = &Linear_NE_exact;
  complex<double> (*expt_ACF_NE_exact_v1)
  (double tp, double tau, double omega_DA, double V, vector<double>& omega_nm, vector<double>& gamma_nm,
   vector<double>& req_nm, vector<double>& shift_NE, double beta, int type) = ACF_NE_exact;

  def("Integrand_NE_exact", expt_Integrand_NE_exact_v1);
  def("Linear_NE_exact", expt_Linear_NE_exact_v1);
  def("ACF_NE_exact", expt_ACF_NE_exact_v1);


  complex<double> (*expt_Integrand_NE_LSC_v1)
  (double tp, double tau, double omega_DA, double omega, double req, double shift, double beta) = &Integrand_NE_LSC;
  complex<double> (*expt_Linear_NE_LSC_v1)
  (double tp, double tau, double gamma, double omega, double req, double shift, double beta) = &Linear_NE_LSC;
  complex<double> (*expt_ACF_NE_LSC_v1)
  (double tp, double tau, double omega_DA, double V, vector<double>& omega_nm, vector<double>& gamma_nm,
   vector<double>& req_nm, vector<double>& shift_NE, double beta, int type) = ACF_NE_LSC;

  def("Integrand_NE_LSC", expt_Integrand_NE_LSC_v1);
  def("Linear_NE_LSC", expt_Linear_NE_LSC_v1);
  def("ACF_NE_LSC", expt_ACF_NE_LSC_v1);


  complex<double> (*expt_Integrand_NE_CAV_v1)
  (double tp, double tau, double omega_DA, double omega, double req, double shift, double beta) = &Integrand_NE_CAV;
  complex<double> (*expt_Linear_NE_CAV_v1)
  (double tp, double tau, double gamma, double omega, double req, double shift, double beta) = &Linear_NE_CAV;
  complex<double> (*expt_ACF_NE_CAV_v1)
  (double tp, double tau, double omega_DA, double V, vector<double>& omega_nm, vector<double>& gamma_nm,
   vector<double>& req_nm, vector<double>& shift_NE, double beta, int type) = ACF_NE_CAV;

  def("Integrand_NE_CAV", expt_Integrand_NE_CAV_v1);
  def("Linear_NE_CAV", expt_Linear_NE_CAV_v1);
  def("ACF_NE_CAV", expt_ACF_NE_CAV_v1);


  complex<double> (*expt_Integrand_NE_CD_v1)
  (double tp, double tau, double omega_DA, double omega, double req, double shift, double beta) = &Integrand_NE_CD;
  complex<double> (*expt_Linear_NE_CD_v1)
  (double tp, double tau, double gamma, double omega, double req, double shift, double beta) = &Linear_NE_CD;
  complex<double> (*expt_ACF_NE_CD_v1)
  (double tp, double tau, double omega_DA, double V, vector<double>& omega_nm, vector<double>& gamma_nm,
   vector<double>& req_nm, vector<double>& shift_NE, double beta, int type) = ACF_NE_CD;

  def("Integrand_NE_CD", expt_Integrand_NE_CD_v1);
  def("Linear_NE_CD", expt_Linear_NE_CD_v1);
  def("ACF_NE_CD", expt_ACF_NE_CD_v1);


  complex<double> (*expt_Integrand_NE_W0_v1)
  (double tp, double tau, double omega_DA, double omega, double req, double shift, double beta) = &Integrand_NE_W0;
  complex<double> (*expt_Linear_NE_W0_v1)
  (double tp, double tau, double gamma, double omega, double req, double shift, double beta) = &Linear_NE_W0;
  complex<double> (*expt_ACF_NE_W0_v1)
  (double tp, double tau, double omega_DA, double V, vector<double>& omega_nm, vector<double>& gamma_nm,
   vector<double>& req_nm, vector<double>& shift_NE, double beta, int type) = ACF_NE_W0;

  def("Integrand_NE_W0", expt_Integrand_NE_W0_v1);
  def("Linear_NE_W0", expt_Linear_NE_W0_v1);
  def("ACF_NE_W0", expt_ACF_NE_W0_v1);


  complex<double> (*expt_Integrand_NE_C0_v1)
  (double tp, double tau, double omega_DA, double omega, double req, double shift, double beta) = &Integrand_NE_C0;
  complex<double> (*expt_Linear_NE_C0_v1)
  (double tp, double tau, double gamma, double omega, double req, double shift, double beta) = &Linear_NE_C0;
  complex<double> (*expt_ACF_NE_C0_v1)
  (double tp, double tau, double omega_DA, double V, vector<double>& omega_nm, vector<double>& gamma_nm,
   vector<double>& req_nm, vector<double>& shift_NE, double beta, int type) = ACF_NE_C0;

  def("Integrand_NE_C0", expt_Integrand_NE_C0_v1);
  def("Linear_NE_C0", expt_Linear_NE_C0_v1);
  def("ACF_NE_C0", expt_ACF_NE_C0_v1);


  double (*expt_NEFGRL_rate_v1)
  (double t, double omega_DA, double V,
   vector<double>& omega_nm, vector<double>& gamma_nm,
   vector<double>& req_nm, vector<double>& shift_NE, int method, double beta, int type, double dtau) = &NEFGRL_rate;

  def("NEFGRL_rate", expt_NEFGRL_rate_v1);


  MATRIX (*expt_NEFGRL_population_v1)
  (double omega_DA, double V,
   vector<double>& omega_nm, vector<double>& gamma_nm,
   vector<double>& req_nm, vector<double>& shift_NE,
   int method, double beta, int type, double dtau, double tmax, double dt) = &NEFGRL_population;

  def("NEFGRL_population", expt_NEFGRL_population_v1);



} // export_fgr_objects()





#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygfgr){
#else
BOOST_PYTHON_MODULE(libfgr){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

  export_fgr_objects();

}

}// namespace libfgr
}// liblibra

