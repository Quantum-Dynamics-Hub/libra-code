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

#include <memory> // for std::auto_ptr<>
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "libhamiltonian_qm.h"

using namespace boost::python;

namespace libhamiltonian{
namespace libhamiltonian_atomistic{
namespace libhamiltonian_qm{


void export_Hamiltonian_QM_objects(){

  //----------- Electronic_Structure ------------
  class_<Electronic_Structure>("Electronic_Structure",init<>())
      .def(init<int>())
      .def("__copy__", &generic__copy__<Electronic_Structure>)
      .def("__deepcopy__", &generic__deepcopy__<Electronic_Structure>)

      .def_readwrite("Norb", &Electronic_Structure::Norb)
      .def_readwrite("Nocc_alp", &Electronic_Structure::Nocc_alp)
      .def_readwrite("Nocc_bet", &Electronic_Structure::Nocc_bet)
      .def_readwrite("Nelec", &Electronic_Structure::Nelec)

      .def_readwrite("bands_alp", &Electronic_Structure::bands_alp)
      .def_readwrite("bands_bet", &Electronic_Structure::bands_bet)
      .def_readwrite("occ_alp", &Electronic_Structure::occ_alp)
      .def_readwrite("occ_bet", &Electronic_Structure::occ_bet)

      .def_readwrite("Mull_orb_pop_net", &Electronic_Structure::Mull_orb_pop_net)
      .def_readwrite("Mull_orb_pop_gross", &Electronic_Structure::Mull_orb_pop_gross)

      .def("set_P_alp", &Electronic_Structure::set_P_alp)
      .def("set_P_bet", &Electronic_Structure::set_P_bet)
      .def("set_P", &Electronic_Structure::set_P)
      .def("get_P_alp", &Electronic_Structure::get_P_alp)
      .def("get_P_bet", &Electronic_Structure::get_P_bet)
      .def("get_P", &Electronic_Structure::get_P)

      .def("set_C_alp", &Electronic_Structure::set_C_alp)
      .def("set_C_bet", &Electronic_Structure::set_C_bet)
      .def("get_C_alp", &Electronic_Structure::get_C_alp)
      .def("get_C_bet", &Electronic_Structure::get_C_bet)


      .def("set_Sao", &Electronic_Structure::set_Sao)
      .def("set_Hao", &Electronic_Structure::set_Hao)
      .def("get_Sao", &Electronic_Structure::get_Sao)
      .def("get_Hao", &Electronic_Structure::get_Hao)


      .def("set_Fao_alp", &Electronic_Structure::set_Fao_alp)
      .def("set_Fao_bet", &Electronic_Structure::set_Fao_bet)
      .def("get_Fao_alp", &Electronic_Structure::get_Fao_alp)
      .def("get_Fao_bet", &Electronic_Structure::get_Fao_bet)


      .def("set_dFao_alp_dP_alp", &Electronic_Structure::set_dFao_alp_dP_alp)
      .def("set_dFao_alp_dP_bet", &Electronic_Structure::set_dFao_alp_dP_bet)
      .def("set_dFao_bet_dP_alp", &Electronic_Structure::set_dFao_bet_dP_alp)
      .def("set_dFao_bet_dP_bet", &Electronic_Structure::set_dFao_bet_dP_bet)
      .def("get_dFao_alp_dP_alp", &Electronic_Structure::get_dFao_alp_dP_alp)
      .def("get_dFao_alp_dP_bet", &Electronic_Structure::get_dFao_alp_dP_bet)
      .def("get_dFao_bet_dP_alp", &Electronic_Structure::get_dFao_bet_dP_alp)
      .def("get_dFao_bet_dP_bet", &Electronic_Structure::get_dFao_bet_dP_bet)


      .def("set_E_alp", &Electronic_Structure::set_E_alp)
      .def("set_E_bet", &Electronic_Structure::set_E_bet)
      .def("get_E_alp", &Electronic_Structure::get_E_alp)
      .def("get_E_bet", &Electronic_Structure::get_E_bet)



  ;


  //----------- INDO -----------------
  void (*expt_indo_core_parameters_v1)
  ( System& syst, vector<AO>& basis_ao, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  int opt, int DF) = &indo_core_parameters;

  void (*expt_Hamiltonian_core_indo_v1)
  ( System& syst, vector<AO>& basis_ao, 
  Control_Parameters& prms, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  MATRIX& Hao, MATRIX& Sao, int DF) = &Hamiltonian_core_indo;

  void (*expt_get_integrals_v1)
  (int i,int j,vector<AO>& basis_ao, double eri_aa, double G1, double F2, double& ii_jj,double& ij_ij) = &get_integrals;

  void (*expt_Hamiltonian_Fock_indo_v1)
  (Electronic_Structure& el, System& syst, vector<AO>& basis_ao,
   Control_Parameters& prms,Model_Parameters& modprms,
   vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map
  ) = &Hamiltonian_Fock_indo;



  def("indo_core_parameters", expt_indo_core_parameters_v1);
  def("Hamiltonian_core_indo", expt_Hamiltonian_core_indo_v1);
  def("get_integrals",expt_get_integrals_v1);
  def("Hamiltonian_Fock_indo",expt_Hamiltonian_Fock_indo_v1);



  //----------- HF -----------------
  void (*expt_Hamiltonian_core_hf_v1)
  ( System& syst, vector<AO>& basis_ao, 
  Control_Parameters& prms, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  MATRIX& Hao, MATRIX& Sao, int DF) = &Hamiltonian_core_hf;


  void (*expt_Hamiltonian_Fock_hf_v1)
  (Electronic_Structure& el, System& syst, vector<AO>& basis_ao,
   Control_Parameters& prms,Model_Parameters& modprms,
   vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map
  ) = &Hamiltonian_Fock_hf;




  def("Hamiltonian_core_hf", expt_Hamiltonian_core_hf_v1);
  def("Hamiltonian_Fock_hf", expt_Hamiltonian_Fock_hf_v1);



  //----------- Hamiltonian_QM -----------------
  void (*expt_Hamiltonian_core_v1)
  ( System& syst, vector<AO>& basis_ao, 
  Control_Parameters& prms, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  MATRIX& Hao, MATRIX& Sao, int DF) = &Hamiltonian_core;

  void (*expt_Hamiltonian_Fock_v1)
  (Electronic_Structure& el, System& syst, vector<AO>& basis_ao,
   Control_Parameters& prms,Model_Parameters& modprms,
   vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map
  ) = &Hamiltonian_Fock;


  def("Hamiltonian_core", expt_Hamiltonian_core_v1);
  def("Hamiltonian_Fock", expt_Hamiltonian_Fock_v1);



  //--------------- SCF.*** ------------------------------------
  double (*expt_scf_v1)
  (Electronic_Structure& el, System& syst, vector<AO>& basis_ao,
   Control_Parameters& prms,Model_Parameters& modprms,
   vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map, int BM) = &scf;

  double (*expt_scf_none_v1)
  (Electronic_Structure& el, System& syst, vector<AO>& basis_ao,
   Control_Parameters& prms,Model_Parameters& modprms,
   vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map, int BM) = &scf_none;

  double (*expt_scf_oda_v1)
  (Electronic_Structure& el, System& syst, vector<AO>& basis_ao,
   Control_Parameters& prms,Model_Parameters& modprms,
   vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map, int BM) = &scf_oda;

  double (*expt_scf_oda_disk_v1)(Electronic_Structure& el, System& syst, vector<AO>& basis_ao,
           Control_Parameters& prms,Model_Parameters& modprms,
           vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map, int BM) = &scf_oda_disk;

  def("scf", expt_scf_v1);
  def("scf_none", expt_scf_none_v1);
  def("scf_oda", expt_scf_oda_v1);
  def("scf_oda_disk", expt_scf_oda_disk_v1);


}


#ifdef CYGWIN
BOOST_PYTHON_MODULE(cyghamiltonian_qm){
#else
BOOST_PYTHON_MODULE(libhamiltonian_qm){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

  export_Hamiltonian_QM_objects();

}




}// namespace libhamiltonian_qm
}// namespace libhamiltonian_atomistic
}// namespace libhamiltonian


