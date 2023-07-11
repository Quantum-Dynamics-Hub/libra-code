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

#if defined(USING_PCH)
#include "../pch.h"
#else
#include <memory> // for std::auto_ptr<>
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>
#endif 

#include "libmodel_parameters.h"


/// liblibra namespace
namespace liblibra{


using namespace boost::python;

namespace libmodel_parameters{


void export_Model_Parameters_objects(){

  class_<HF_integrals>("HF_integrals",init<>())
      .def("set_JK_values", &HF_integrals::set_JK_values)
      .def("get_JK_values", &HF_integrals::get_JK_values)

  ;

  class_< HF_integralsList >("HF_integralsList")
      .def(vector_indexing_suite< HF_integralsList >())
  ;
  class_< HF_integralsMap >("HF_integralsMap")
      .def(vector_indexing_suite< HF_integralsMap >())
  ;


  class_<EHT_K>("EHT_K",init<>())
      .def("set_PSPS_value", &EHT_K::set_PSPS_value)

      .def("set_PPa_value", &EHT_K::set_PPa_value)
      .def("set_PP0_value", &EHT_K::set_PP0_value)
      .def("set_PP1_value", &EHT_K::set_PP1_value)
      .def("set_PP2_value", &EHT_K::set_PP2_value)
      .def("get_PPa_value", &EHT_K::get_PPa_value)
      .def("get_PP0_value", &EHT_K::get_PP0_value)
      .def("get_PP1_value", &EHT_K::get_PP1_value)
      .def("get_PP2_value", &EHT_K::get_PP2_value)

      .def("set_K_value",  &EHT_K::set_K_value)
      .def("get_K_value",  &EHT_K::get_K_value)
      .def("set_C_value",  &EHT_K::set_C_value)
      .def("get_C_value",  &EHT_K::get_C_value)

      .def("show",&EHT_K::show)  
      .def_readwrite("data", &EHT_K::data)
  ;


  class_< EHT_KList >("EHT_KList")
      .def(vector_indexing_suite< EHT_KList >())
  ;
  class_< EHT_KMap >("EHT_KMap")
      .def(vector_indexing_suite< EHT_KMap >())
  ;



  class_<mEHT_K>("mEHT_K",init<>())

      .def_readwrite("size",   &mEHT_K::size)
//      .def_readwrite("eht_K",  &mEHT_K::eht_K)
//      .def_readwrite("eht_C",  &mEHT_K::eht_C)
      .def_readwrite("eht_PPa", &mEHT_K::eht_PPa)
      .def_readwrite("eht_PP0", &mEHT_K::eht_PP0)
      .def_readwrite("eht_PP1", &mEHT_K::eht_PP1)
      .def_readwrite("eht_PP2", &mEHT_K::eht_PP2)


      .def("set_mapping", &mEHT_K::set_mapping)
      .def("set_mapping1", &mEHT_K::set_mapping1)

      .def("get_K_value",  &mEHT_K::get_K_value)
      .def("get_C_value",  &mEHT_K::get_C_value)
  ;


  class_< mEHT_KList >("mEHT_KList")
      .def(vector_indexing_suite< mEHT_KList >())
  ;
  class_< mEHT_KMap >("mEHT_KMap")
      .def(vector_indexing_suite< mEHT_KMap >())
  ;



  class_<pElement>("pElement",init<>())
//      .def(_s)
//  void _set(std::string en,int z){ elt_name = en; Z = z;  }
//  void _set(std::string en,int z,int nv){ elt_name = en; Z = z; Nval = nv; }
//  void _set(std::string en,int z,int nv,double zeff){ elt_name = en; Z = z; Nval = nv; Zeff = zeff; }
//  void _set(std::string en,int z,int nv,double zeff,map<std::string, double> ip){ elt_name = en; Z = z; Nval = nv; IP = ip; Zeff = zeff;}
//  void set_mass(double m_){ mass = m_; }


      .def_readwrite("elt_name",   &pElement::elt_name)
      .def_readwrite("Z",   &pElement::Z)
      .def_readwrite("PQN",   &pElement::PQN)
      .def_readwrite("Nval",   &pElement::Nval)
      .def_readwrite("Zeff",   &pElement::Zeff)
      .def_readwrite("mass",   &pElement::mass)

      .def_readwrite("IP",   &pElement::IP)
      .def_readwrite("EA",   &pElement::EA)

      .def_readwrite("Nquant",   &pElement::Nquant)
      .def_readwrite("Nzeta",   &pElement::Nzeta)
      .def_readwrite("zetas",   &pElement::zetas)
      .def_readwrite("coeffs",   &pElement::coeffs)

      .def_readwrite("J_param1",   &pElement::J_param1)
      .def_readwrite("J_param2",   &pElement::J_param2)
      .def_readwrite("J_param3",   &pElement::J_param3)
      .def_readwrite("J_param4",   &pElement::J_param4)

      .def_readwrite("G1",   &pElement::G1)
      .def_readwrite("F2",   &pElement::F2)
      .def_readwrite("beta0",  &pElement::beta0)
      .def_readwrite("Zeta",   &pElement::Zeta)

  ;


  class_< pElementList >("pElementList")
      .def(vector_indexing_suite< pElementList >())
  ;
  class_< pElementMap >("pElementMap")
      .def(vector_indexing_suite< pElementMap >())
  ;

  class_< map<std::string,pElement> >("StringpElementMap")
      .def(map_indexing_suite<map<std::string,pElement> >())
  ;


  class_<OrbParams>("OrbParams",init<>())

      .def_readwrite("IP",   &OrbParams::IP)
      .def_readwrite("EA",   &OrbParams::EA)
      .def_readwrite("Nquant",   &OrbParams::Nquant)
      .def_readwrite("Nzeta",   &OrbParams::Nzeta)
      .def_readwrite("coeffs",   &OrbParams::zetas)

      .def_readwrite("J_param1",   &OrbParams::J_param1)
      .def_readwrite("J_param2",   &OrbParams::J_param2)
      .def_readwrite("J_param3",   &OrbParams::J_param3)
      .def_readwrite("J_param4",   &OrbParams::J_param4)

      .def_readwrite("G1",   &OrbParams::G1)
      .def_readwrite("F2",   &OrbParams::F2)
      .def_readwrite("beta0",   &OrbParams::beta0)

  ;

  class_< OrbParamsList >("OrbParamsList")
      .def(vector_indexing_suite< OrbParamsList >())
  ;
  class_< OrbParamsMap >("OrbParamsMap")
      .def(vector_indexing_suite< OrbParamsMap >())
  ;




  class_<Model_Parameters>("Model_Parameters",init<>())
      .def("__copy__", &generic__copy__<Model_Parameters>)
      .def("__deepcopy__", &generic__deepcopy__<Model_Parameters>)

      .def_readwrite("PT", &Model_Parameters::PT)
      .def_readwrite("orb_params", &Model_Parameters::orb_params)
      .def_readwrite("eht_k", &Model_Parameters::eht_k)
      .def_readwrite("meht_k", &Model_Parameters::meht_k)
      .def_readwrite("hf_int", &Model_Parameters::hf_int)

      .def("set_PT_mapping", &Model_Parameters::set_PT_mapping)

  ;



//  void (*expt_set_parameters_eht_v1)
//  (Control_Parameters&, Model_Parameters&) = &set_parameters_eht;

  void (*expt_set_parameters_hf_v1)
  (Control_Parameters&, Model_Parameters&, vector<AO>&) = &set_parameters_hf;

  void (*expt_set_parameters_indo_v1)
  (Control_Parameters&, Model_Parameters&) = &set_parameters_indo;

  void (*expt_set_parameters_eht_v1)
  (Control_Parameters& prms, Model_Parameters& modprms) = &set_parameters_eht;

  void (*expt_set_parameters_eht_mapping_v1)
  (Model_Parameters& modprms, const vector<AO>& basis_ao) = &set_parameters_eht_mapping;

  void (*expt_set_parameters_eht_mapping1_v1)
  (Model_Parameters& modprms, int nat, vector<std::string>& mol_at_types) = &set_parameters_eht_mapping1;


  def("set_parameters_hf", expt_set_parameters_hf_v1);
  def("set_parameters_indo", expt_set_parameters_indo_v1);
  def("set_parameters_eht", expt_set_parameters_eht_v1);

  def("set_parameters_eht_mapping", expt_set_parameters_eht_mapping_v1);
  def("set_parameters_eht_mapping1", expt_set_parameters_eht_mapping1_v1);



}


#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygmodel_parameters){
#else
BOOST_PYTHON_MODULE(libmodel_parameters){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

  export_Model_Parameters_objects();

}



}// namespace libmodel_parameters
}// namespace liblibra


