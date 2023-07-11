/*********************************************************************************
* Copyright (C) 2018-2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
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

#include "libivr.h"

//#define BOOST_PYTHON_MAX_ARITY 30


/// liblibra namespace
namespace liblibra{


using namespace boost::python;
using namespace libmeigen;
using namespace librandom;
using namespace liblinalg;

namespace libivr{


void export_ivr_objects(){


  class_<ivr_params>("ivr_params",init<int>())
    .def(init<const ivr_params&>()) 

    .def("set_qIn",&ivr_params::set_qIn)
    .def("set_pIn",&ivr_params::set_pIn)
    .def("set_Width0",&ivr_params::set_Width0)
    .def("set_WidthT",&ivr_params::set_WidthT)
    .def("set_TuningQ",&ivr_params::set_TuningQ)
    .def("set_TuningP",&ivr_params::set_TuningP)

    .def("get_qIn",&ivr_params::get_qIn)
    .def("get_pIn",&ivr_params::get_pIn)
    .def("get_Width0",&ivr_params::get_Width0)
    .def("get_WidthT",&ivr_params::get_WidthT)
    .def("get_invWidth0",&ivr_params::get_invWidth0)
    .def("get_invWidthT",&ivr_params::get_invWidthT)
    .def("get_TuningQ",&ivr_params::get_TuningQ)
    .def("get_TuningP",&ivr_params::get_TuningP)
    .def("get_invTuningQ",&ivr_params::get_invTuningQ)
    .def("get_invTuningP",&ivr_params::get_invTuningP)

    .def_readwrite("Ndof",&ivr_params::Ndof)
         
  ;


  class_<ivr_observable>("ivr_observables",init<>())
    .def(init<int, int>()) 

    .def_readwrite("observable_type",&ivr_observable::observable_type)
    .def_readwrite("observable_label",&ivr_observable::observable_label)

  ;



  ///============ Sampling Initial Distributions ===========================
  ///  In ivr_sampling.cpp

  MATRIX (*expt_ivr_Husimi_v1)(MATRIX& qIn, MATRIX& pIn, MATRIX& Width0, Random& rnd) = &ivr_Husimi;
  vector<MATRIX> (*expt_ivr_Husimi_v2)(MATRIX& qIn, MATRIX& pIn, MATRIX& Width0, Random& rnd, int sample_size) = &ivr_Husimi;

  def("ivr_Husimi", expt_ivr_Husimi_v1);
  def("ivr_Husimi", expt_ivr_Husimi_v2);



  MATRIX (*expt_ivr_LSC_v1)(MATRIX& qIn, MATRIX& pIn, MATRIX& Width0, Random& rnd) = &ivr_LSC;
  vector<MATRIX> (*expt_ivr_LSC_v2)(MATRIX& qIn, MATRIX& pIn, MATRIX& Width0, Random& rnd, int sample_size) = &ivr_LSC;

  def("ivr_LSC", expt_ivr_LSC_v1);
  def("ivr_LSC", expt_ivr_LSC_v2);



  MATRIX (*expt_ivr_DHK_v1)(MATRIX& qIn, MATRIX& pIn, MATRIX& Width0, Random& rnd) = &ivr_DHK;
  vector<MATRIX> (*expt_ivr_DHK_v2)(MATRIX& qIn, MATRIX& pIn, MATRIX& Width0, Random& rnd, int sample_size) = &ivr_DHK;

  def("ivr_DHK", expt_ivr_DHK_v1);
  def("ivr_DHK", expt_ivr_DHK_v2);


  MATRIX (*expt_ivr_FB_MQC_v1)(MATRIX& qIn, MATRIX& pIn, MATRIX& Width0, MATRIX& TuningQ, MATRIX& TuningP, int flag, Random& rnd) = &ivr_FB_MQC;
  vector<MATRIX> (*expt_ivr_FB_MQC_v2)(MATRIX& qIn, MATRIX& pIn, MATRIX& Width0, MATRIX& TuningQ, MATRIX& TuningP, int flag, Random& rnd, int sample_size) = &ivr_FB_MQC;

  def("ivr_FB_MQC", expt_ivr_FB_MQC_v1);
  def("ivr_FB_MQC", expt_ivr_FB_MQC_v2);


  MATRIX (*expt_ivr_FF_MQC_v1)(MATRIX& qIn, MATRIX& pIn, MATRIX& Width0, MATRIX& TuningQ, MATRIX& TuningP, Random& rnd) = &ivr_FF_MQC;
  vector<MATRIX> (*expt_ivr_FF_MQC_v2)(MATRIX& qIn, MATRIX& pIn, MATRIX& Width0, MATRIX& TuningQ, MATRIX& TuningP, Random& rnd, int sample_size) = &ivr_FF_MQC;

  def("ivr_FF_MQC", expt_ivr_FF_MQC_v1);
  def("ivr_FF_MQC", expt_ivr_FF_MQC_v2);


  ///============ Matrix Elements and Overlaps ===========================
  ///  In ivr_matrix_elements.cpp


  complex<double> (*expt_CS_overlap_v1)
  (MATRIX& q, MATRIX& p, MATRIX& qIn, MATRIX& pIn, 
  MATRIX& Width0, MATRIX& invWidth0) = &CS_overlap;

  def("CS_overlap", expt_CS_overlap_v1);


  complex<double> (*expt_mat_elt_FB_B_v1)(MATRIX& q, MATRIX& p, int opt, int lab) = &mat_elt_FB_B;
  def("mat_elt_FB_B", expt_mat_elt_FB_B_v1);


  complex<double> (*expt_mat_elt_FF_B_v1)
  (MATRIX& q, MATRIX& p, MATRIX& qp, MATRIX& pp, 
  MATRIX& WidthT, MATRIX& invWidthT, int opt, int lab) = &mat_elt_FF_B;

  def("mat_elt_FF_B", expt_mat_elt_FF_B_v1);


  complex<double> (*expt_mat_elt_HUS_B_v1)(MATRIX& q, MATRIX& p, int opt, int lab) = &mat_elt_HUS_B;

  def("mat_elt_HUS_B", expt_mat_elt_HUS_B_v1);


  complex<double> (*expt_mat_elt_LSC_B_v1)(MATRIX& q, MATRIX& p, int opt, int lab) = &mat_elt_LSC_B;

  def("mat_elt_LSC_B", expt_mat_elt_LSC_B_v1);


  ///============ SC Prefactors  ===========================
  ///  In ivr_prefactors.cpp
  complex<double> (*expt_MQC_prefactor_FB_G_v1)
  (vector<CMATRIX>& Mfwd, vector<CMATRIX>& Mbck, ivr_params& prms) = &MQC_prefactor_FB_G;

  def("MQC_prefactor_FB_G", expt_MQC_prefactor_FB_G_v1);


  complex<double> (*expt_MQC_prefactor_FF_G_v1)
  (vector<CMATRIX>& Mfwd, vector<CMATRIX>& Mbck, ivr_params& prms) = &MQC_prefactor_FF_G;

  def("MQC_prefactor_FF_G", expt_MQC_prefactor_FF_G_v1);


  vector<complex<double> > (*expt_DHK_prefactor_v1)
  (vector<CMATRIX>& Mfwd, vector<CMATRIX>& Mbck, ivr_params& prms) = &DHK_prefactor;

  def("DHK_prefactor", expt_DHK_prefactor_v1);


  ///============ Propagators  ===========================
  ///  In ivr_propagators.cpp
  void (*expt_Integrator_v1)(MATRIX& q, MATRIX& p, vector<MATRIX>& M, double& action, MATRIX& mass, double dt) = &Integrator;
  def("Integrator", expt_Integrator_v1);



  ///============ TCF calculators  ===========================
  ///  In ivr_timecorr.cpp
  void (*expt_compute_tcf_v1)(vector< complex<double> >& TCF, vector<int>& MCnum,
                   vector<MATRIX>& q, vector<MATRIX>& p, vector<int>& status,
                   int ivr_opt, int observable_type, int observable_label) = &compute_tcf;
/*
  void (*expt_compute_tcf_v2)(vector< complex<double> >& TCF, vector<int>& MCnum,
                 MATRIX& qIn, MATRIX& pIn,
                 MATRIX& Width0, MATRIX& invWidth0, MATRIX& WidthT, MATRIX& invWidthT,
                 MATRIX& TuningQ, MATRIX& invTuningQ, MATRIX& TuningP, MATRIX& invTuningP,
                 vector<MATRIX>& q,  vector<MATRIX>& p,  vector<int>& status, vector<double>& action,vector< vector<MATRIX> >& Mono,
                 vector<MATRIX>& qp, vector<MATRIX>& pp, vector<int>& statusp,vector<double>& actionp,vector< vector<MATRIX> >& Monop,
                 int ivr_opt, int observable_type, int observable_label) = &compute_tcf;
*/
  def("compute_tcf", expt_compute_tcf_v1);
//  def("compute_tcf", expt_compute_tcf_v2);

  void (*expt_normalize_tcf_v1)(vector< complex<double> >& TCF, vector<int>& MCnum) = &normalize_tcf;
  def("normalize_tcf", expt_normalize_tcf_v1);



} // export_ivr_objects()





#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygivr){
#else
BOOST_PYTHON_MODULE(libivr){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

  export_ivr_objects();

}

}// namespace libivr
}// liblibra

