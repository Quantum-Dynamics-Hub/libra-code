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
  \file libhamiltonian.cpp
  \brief The file implements Python export function
    
*/

#if defined(USING_PCH)
#include "../pch.h"
#else
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#endif 

#include "libmodels.h"


/// liblibra namespace
namespace liblibra{

using namespace boost::python;

/// libhamiltonian namespace
namespace libmodels{


void export_models_objects(){
/** 
  \brief Exporter of the libhamiltonian_model classes and functions

*/

  boost::python::list (*expt_SAC_Ham1)(double,boost::python::list) = &SAC_Ham;
  def("SAC_Ham", expt_SAC_Ham1);

  boost::python::list (*expt_DAC_Ham1)(double,boost::python::list) = &DAC_Ham;
  def("DAC_Ham", expt_DAC_Ham1);

  boost::python::list (*expt_ECWR_Ham1)(double,boost::python::list) = &ECWR_Ham;
  def("ECWR_Ham", expt_ECWR_Ham1);

  boost::python::list (*expt_Marcus_Ham1)(double,boost::python::list) = &Marcus_Ham;
  def("Marcus_Ham", expt_Marcus_Ham1);

  boost::python::list (*expt_SEXCH_Ham1)(double,boost::python::list) = &SEXCH_Ham;
  def("SEXCH_Ham", expt_SEXCH_Ham1);

  boost::python::list (*expt_Rabi2_Ham1)(double,boost::python::list) = &Rabi2_Ham;
  def("Rabi2_Ham", expt_Rabi2_Ham1);

  boost::python::list (*expt_sin_Ham1)(double,boost::python::list) = &sin_Ham;
  def("sin_Ham", expt_sin_Ham1);

  boost::python::list (*expt_sin_2D_Ham1)(double,double,boost::python::list) = &sin_2D_Ham;
  def("sin_2D_Ham", expt_sin_2D_Ham1);

  boost::python::list (*expt_cubic_Ham1)(double,boost::python::list) = &cubic_Ham;
  def("cubic_Ham", expt_cubic_Ham1);

  boost::python::list (*expt_double_well_Ham1)(double,boost::python::list) = &double_well_Ham;
  def("double_well_Ham", expt_double_well_Ham1);


  void (*expt_model_SAC_v1)(CMATRIX& Hdia, CMATRIX& Sdia, vector<CMATRIX>& d1ham_dia, vector<CMATRIX>& dc1_dia,
               vector<double>& q, vector<double>& params) = &model_SAC;
  def("model_SAC", expt_model_SAC_v1);

  void (*expt_model_DAC_v1)(CMATRIX& Hdia, CMATRIX& Sdia, vector<CMATRIX>& d1ham_dia, vector<CMATRIX>& dc1_dia,
               vector<double>& q, vector<double>& params) = &model_DAC;
  def("model_DAC", expt_model_DAC_v1);

  void (*expt_model_ECWR_v1)(CMATRIX& Hdia, CMATRIX& Sdia, vector<CMATRIX>& d1ham_dia, vector<CMATRIX>& dc1_dia,
               vector<double>& q, vector<double>& params) = &model_ECWR;
  def("model_ECWR", expt_model_ECWR_v1);


  vector<double> (*expt_set_params_1S_1D_poly4_v1)(std::string model) = &set_params_1S_1D_poly4;

  void (*expt_model_1S_1D_poly2_v1)(CMATRIX& Hdia, CMATRIX& Sdia, vector<CMATRIX>& d1ham_dia, vector<CMATRIX>& dc1_dia,
                          vector<double>& q, vector<double>& params) = &model_1S_1D_poly2;

  void (*expt_model_1S_1D_poly4_v1)(CMATRIX& Hdia, CMATRIX& Sdia, vector<CMATRIX>& d1ham_dia, vector<CMATRIX>& dc1_dia,
                          vector<double>& q, vector<double>& params) = &model_1S_1D_poly4;

  def("set_params_1S_1D_poly4", expt_set_params_1S_1D_poly4_v1);
  def("model_1S_1D_poly2", expt_model_1S_1D_poly2_v1);
  def("model_1S_1D_poly4", expt_model_1S_1D_poly4_v1);



  vector<double> (*expt_set_params_2S_1D_sin_v1)(std::string model) = &set_params_2S_1D_sin;

  void (*expt_model_2S_1D_sin_v1)(CMATRIX& Hdia, CMATRIX& Sdia, vector<CMATRIX>& d1ham_dia, vector<CMATRIX>& dc1_dia,
                     vector<double>& q, vector<double>& params) = &model_2S_1D_sin;

  def("set_params_2S_1D_sin", expt_set_params_2S_1D_sin_v1);
  def("model_2S_1D_sin", expt_model_2S_1D_sin_v1);



  vector<double> (*expt_set_params_2S_2D_sin_v1)(std::string model) = &set_params_2S_2D_sin;

  void (*expt_model_2S_2D_sin_v1)(CMATRIX& Hdia, CMATRIX& Sdia, vector<CMATRIX>& d1ham_dia, vector<CMATRIX>& dc1_dia,
                     vector<double>& q, vector<double>& params) = &model_2S_2D_sin;

  def("set_params_2S_2D_sin", expt_set_params_2S_2D_sin_v1);
  def("model_2S_2D_sin", expt_model_2S_2D_sin_v1);


  vector<double> (*expt_set_params_2S_1D_tanh_v1)(std::string model) = &set_params_2S_1D_tanh;

  void (*expt_model_2S_1D_tanh_v1)(CMATRIX& Hdia, CMATRIX& Sdia, vector<CMATRIX>& d1ham_dia, vector<CMATRIX>& dc1_dia,
                     vector<double>& q, vector<double>& params) = &model_2S_1D_tanh;

  def("set_params_2S_1D_tanh", expt_set_params_2S_1D_tanh_v1);
  def("model_2S_1D_tanh", expt_model_2S_1D_tanh_v1);



}


#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygmodels){
#else
BOOST_PYTHON_MODULE(libmodels){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/

  export_models_objects();

}


}// namespace libmodels
}// liblibra

