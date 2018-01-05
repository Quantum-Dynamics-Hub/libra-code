/*********************************************************************************
* Copyright (C) 2018 Alexey V. Akimov
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

#include "libivr.h"

/// liblibra namespace
namespace liblibra{


using namespace boost::python;
using namespace libmeigen;
using namespace librandom;
using namespace liblinalg;

namespace libivr{


void export_ivr_objects(){

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


  complex<double> (*expt_CS_overlap_v1)(MATRIX& q, MATRIX& p, MATRIX& qIn, MATRIX& pIn, MATRIX& Width0, MATRIX& invWidth0) = &CS_overlap;
  def("CS_overlap", expt_CS_overlap_v1);


  complex<double> (*expt_mat_elt_FB_B_v1)(MATRIX& q, MATRIX& p, int opt, int lab) = &mat_elt_FB_B;
  def("mat_elt_FB_B", expt_mat_elt_FB_B_v1);


  complex<double> (*expt_mat_elt_FF_B_v1)(MATRIX& q, MATRIX& p, MATRIX& qp, MATRIX& pp, MATRIX& WidthT, MATRIX& invWidthT, int opt, int lab) = &mat_elt_FF_B;
  def("mat_elt_FF_B", expt_mat_elt_FF_B_v1);


  complex<double> (*expt_mat_elt_HUS_B_v1)(MATRIX& q, MATRIX& p, int opt, int lab) = &mat_elt_HUS_B;
  def("mat_elt_HUS_B", expt_mat_elt_HUS_B_v1);


  complex<double> (*expt_mat_elt_LSC_B_v1)(MATRIX& q, MATRIX& p, int opt, int lab) = &mat_elt_LSC_B;
  def("mat_elt_LSC_B", expt_mat_elt_LSC_B_v1);


  ///============ SC Prefactors  ===========================
  ///  In ivr_prefactors.cpp
  complex<double> (*expt_MQC_prefactor_FB_G_v1)
  (vector<CMATRIX>& Mfwd, vector<CMATRIX>& Mbck, 
   CMATRIX& Width0,  CMATRIX& invWidth0,  CMATRIX& WidthT,  CMATRIX& invWidthT,
   CMATRIX& TuningQ, CMATRIX& invTuningQ, CMATRIX& TuningP, CMATRIX& invTuningP
  ) = &MQC_prefactor_FB_G;

  def("MQC_prefactor_FB_G", expt_MQC_prefactor_FB_G_v1);


  complex<double> (*expt_MQC_prefactor_FF_G_v1)
  (vector<CMATRIX>& Mfwd, vector<CMATRIX>& Mbck, 
   CMATRIX& Width0,  CMATRIX& invWidth0,  CMATRIX& WidthT,  CMATRIX& invWidthT,
   CMATRIX& TuningQ, CMATRIX& invTuningQ, CMATRIX& TuningP, CMATRIX& invTuningP
  ) = &MQC_prefactor_FF_G;

  def("MQC_prefactor_FF_G", expt_MQC_prefactor_FF_G_v1);


  vector<complex<double> > (*expt_DHK_prefactor_v1)
  (vector<CMATRIX>& Mfwd, vector<CMATRIX>& Mbck, 
   CMATRIX& Width0,  CMATRIX& invWidth0,  CMATRIX& WidthT,  CMATRIX& invWidthT,
   CMATRIX& TuningQ, CMATRIX& invTuningQ, CMATRIX& TuningP, CMATRIX& invTuningP
  ) = &DHK_prefactor;

  def("DHK_prefactor", expt_DHK_prefactor_v1);


  ///============ Propagators  ===========================
  ///  In ivr_propagators.cpp
  void (*expt_Integrator_v1)(MATRIX& q, MATRIX& p, vector<MATRIX>& M, double& action, MATRIX& mass, double dt) = &Integrator;
  def("Integrator", expt_Integrator_v1);


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

