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
  \file libhamiltonian_qm.cpp
  \brief The file implements Python export function
    
*/

#define BOOST_PYTHON_MAX_ARITY 30

//#if defined(USING_PCH)
//#include "../../../pch.h"
//#else
#include <memory> // for std::auto_ptr<>
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
//#endif 

#include "libhamiltonian_qm.h"


/// liblibra namespace
namespace liblibra{


using namespace boost::python;

namespace libatomistic{

/// libhamiltonian_qm namespace
namespace libhamiltonian_qm{


using namespace libbasis_setups;
using namespace libcontrol_parameters;
using namespace libmodel_parameters;




void export_hamiltonian_qm_objects(){
/** 
  \brief Exporter of the libhamiltonian_qm classes and functions

*/


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

      .def("get_bands_alp", &Electronic_Structure::get_bands_alp)
      .def("get_bands_bet", &Electronic_Structure::get_bands_bet)
      .def("get_occ_alp", &Electronic_Structure::get_occ_alp)
      .def("get_occ_bet", &Electronic_Structure::get_occ_bet)


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
  vector<int> (*expt_compute_sorb_indices_v1)
  ( int sz, vector<AO>& basis_ao, vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map
  ) = &compute_sorb_indices;

  void (*expt_compute_indo_core_parameters_v1)
  ( System& syst, vector<AO>& basis_ao, Model_Parameters& modprms,
    vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
    vector<int>& sorb_indx,
    int opt, int a, int b, double& eri, double& V_AB
  ) = &compute_indo_core_parameters;

  void (*expt_compute_indo_core_parameters_derivs_v1)
  ( System& syst, vector<AO>& basis_ao, Model_Parameters& modprms,
    vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
    vector<int>& sorb_indx,
    int opt, int a, int b, int c, VECTOR& deri, VECTOR& dV_AB
  ) = &compute_indo_core_parameters_derivs;

  void (*expt_indo_core_parameters_v1)
  ( System& syst, vector<AO>& basis_ao, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  int opt, int DF) = &indo_core_parameters;



  void (*expt_Hamiltonian_core_indo_v1)
  ( System& syst, vector<AO>& basis_ao, 
  Control_Parameters& prms, Model_Parameters& modprms,
  vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
  MATRIX& Hao, MATRIX& Sao, int DF) = &Hamiltonian_core_indo;

  void (*expt_Hamiltonian_core_deriv_indo_v1)
  ( System& syst, vector<AO>& basis_ao, 
    Control_Parameters& prms, Model_Parameters& modprms,
    vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
    MATRIX& Hao, MATRIX& Sao, int DF,
    int c,
    MATRIX& dHao_dx, MATRIX& dHao_dy, MATRIX& dHao_dz, 
    MATRIX& dSao_dx, MATRIX& dSao_dy, MATRIX& dSao_dz
  ) = &Hamiltonian_core_deriv_indo;



  void (*expt_get_integrals_v1)
  (int i,int j,vector<AO>& basis_ao, double eri_aa, double G1, double F2, double& ii_jj,double& ij_ij) = &get_integrals;

  void (*expt_Hamiltonian_Fock_indo_v1)
  (Electronic_Structure& el, System& syst, vector<AO>& basis_ao,
   Control_Parameters& prms,Model_Parameters& modprms,
   vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map
  ) = &Hamiltonian_Fock_indo;

  void (*expt_Hamiltonian_Fock_derivs_indo_v1)
  ( Electronic_Structure& el, System& syst, vector<AO>& basis_ao,
    Control_Parameters& prms, Model_Parameters& modprms,
    vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
    int c, 
    MATRIX& dHao_dx,     MATRIX& dHao_dy,     MATRIX& dHao_dz,
    MATRIX& dFao_alp_dx, MATRIX& dFao_alp_dy, MATRIX& dFao_alp_dz,
    MATRIX& dFao_bet_dx, MATRIX& dFao_bet_dy, MATRIX& dFao_bet_dz
  ) = &Hamiltonian_Fock_derivs_indo;



  def("compute_sorb_indices", expt_compute_sorb_indices_v1);
  def("compute_indo_core_parameters", expt_compute_indo_core_parameters_v1);
  def("compute_indo_core_parameters_derivs", expt_compute_indo_core_parameters_derivs_v1);
  def("indo_core_parameters", expt_indo_core_parameters_v1);

  def("Hamiltonian_core_indo", expt_Hamiltonian_core_indo_v1);
  def("Hamiltonian_core_deriv_indo", expt_Hamiltonian_core_deriv_indo_v1);

  def("get_integrals",expt_get_integrals_v1);
  def("Hamiltonian_Fock_indo",expt_Hamiltonian_Fock_indo_v1);
  def("Hamiltonian_Fock_derivs_indo",expt_Hamiltonian_Fock_derivs_indo_v1);


  //----------- EHT -----------------

  void (*expt_Hamiltonian_core_eht_v1)
  ( System& syst, vector<AO>& basis_ao, 
    Control_Parameters& prms, Model_Parameters& modprms,
    vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
    MATRIX& Hao, MATRIX& Sao, int DF
  ) = &Hamiltonian_core_eht;

  void (*expt_Hamiltonian_core_deriv_eht_v1)
  ( System& syst, vector<AO>& basis_ao, 
    Control_Parameters& prms, Model_Parameters& modprms,
    vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
    MATRIX& Hao, MATRIX& Sao, int DF,
    int c,
    MATRIX& dHao_dx, MATRIX& dHao_dy, MATRIX& dHao_dz, 
    MATRIX& dSao_dx, MATRIX& dSao_dy, MATRIX& dSao_dz
  ) = &Hamiltonian_core_deriv_eht;

  void (*expt_Hamiltonian_Fock_eht_v1)
  ( Electronic_Structure& el, System& syst, vector<AO>& basis_ao,
    Control_Parameters& prms, Model_Parameters& modprms,
    vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map
  ) = &Hamiltonian_Fock_eht;

  def("Hamiltonian_core_eht", expt_Hamiltonian_core_eht_v1);
  def("Hamiltonian_core_deriv_eht", expt_Hamiltonian_core_deriv_eht_v1);
  def("Hamiltonian_Fock_eht", expt_Hamiltonian_Fock_eht_v1);




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


  double (*expt_energy_and_forces_v1)
  ( Electronic_Structure& el, System& syst, vector<AO>& basis_ao,
    Control_Parameters& prms,Model_Parameters& modprms,
    vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map
  ) = &energy_and_forces;


  void (*expt_derivative_couplings_v1)
  ( Electronic_Structure& el, System& syst, vector<AO>& basis_ao,
    Control_Parameters& prms,Model_Parameters& modprms,
    vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
    MATRIX& Hao, MATRIX& Sao,   int Norb, int at_indx, 
    int x_period, int y_period, int z_period, VECTOR& t1, VECTOR& t2, VECTOR& t3,
    MATRIX& Dmo_a_x, MATRIX& Dmo_a_y, MATRIX& Dmo_a_z,
    MATRIX& Dmo_b_x, MATRIX& Dmo_b_y, MATRIX& Dmo_b_z
  ) = &derivative_couplings;

  void (*expt_derivative_couplings_v2)
  ( Electronic_Structure& el, System& syst, vector<AO>& basis_ao,
    Control_Parameters& prms,Model_Parameters& modprms,
    vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
    MATRIX& Hao, MATRIX& Sao,   int Norb, int at_indx, 
    int x_period, int y_period, int z_period, VECTOR& t1, VECTOR& t2, VECTOR& t3,
    MATRIX& Dmo_a_x, MATRIX& Dmo_a_y, MATRIX& Dmo_a_z,
    MATRIX& Dmo_b_x, MATRIX& Dmo_b_y, MATRIX& Dmo_b_z,
    MATRIX& dEa_dx,  MATRIX& dEa_dy,  MATRIX& dEa_dz,
    MATRIX& dEb_dx,  MATRIX& dEb_dy,  MATRIX& dEb_dz
  ) = &derivative_couplings1;

  void (*expt_derivative_couplings_v3)
  ( Electronic_Structure& el, System& syst, vector<AO>& basis_ao,
    Control_Parameters& prms,Model_Parameters& modprms,
    vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
    MATRIX& Hao, MATRIX& Sao,   int Norb, int at_indx, 
    int x_period, int y_period, int z_period, VECTOR& t1, VECTOR& t2, VECTOR& t3,
    MATRIX& Dmo_a_x, MATRIX& Dmo_a_y, MATRIX& Dmo_a_z,
    MATRIX& Dmo_b_x, MATRIX& Dmo_b_y, MATRIX& Dmo_b_z
  ) = &derivative_couplings1;




  VECTOR (*expt_force_v1)
  ( Electronic_Structure& el, System& syst, vector<AO>& basis_ao,
    Control_Parameters& prms,Model_Parameters& modprms,
    vector< vector<int> >& atom_to_ao_map, vector<int>& ao_to_atom_map,
    MATRIX& Hao, MATRIX& Sao, int Norb, int at_indx, 
    int x_period, int y_period, int z_period, VECTOR& t1, VECTOR& t2, VECTOR& t3
  ) = &force;





  def("Hamiltonian_core", expt_Hamiltonian_core_v1);
  def("Hamiltonian_Fock", expt_Hamiltonian_Fock_v1);
  def("derivative_couplings", expt_derivative_couplings_v1);
  def("derivative_couplings", expt_derivative_couplings_v2);
  def("derivative_couplings", expt_derivative_couplings_v3);
  def("force", expt_force_v1);
  def("energy_and_forces", expt_energy_and_forces_v1);



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



  class_<listHamiltonian_QM>("listHamiltonian_QM",init<>())
      .def(init<std::string, System&>())
      .def(init<const listHamiltonian_QM&>())
      .def("__copy__", &generic__copy__<listHamiltonian_QM>)
      .def("__deepcopy__", &generic__deepcopy__<listHamiltonian_QM>)

      .def_readwrite("Norb", &listHamiltonian_QM::Norb)
      .def_readwrite("Nelec", &listHamiltonian_QM::Nelec)
      .def_readwrite("prms", &listHamiltonian_QM::prms)
      .def_readwrite("modprms", &listHamiltonian_QM::modprms)
      .def_readwrite("basis_ao", &listHamiltonian_QM::basis_ao)
      .def_readwrite("atom_to_ao_map", &listHamiltonian_QM::atom_to_ao_map)
      .def_readwrite("ao_to_atom_map", &listHamiltonian_QM::ao_to_atom_map)


      .def("init", &listHamiltonian_QM::init)
      .def("compute_scf", &listHamiltonian_QM::compute_scf)
      .def("get_parameters_from_file", &listHamiltonian_QM::get_parameters_from_file)
      .def("get_electronic_structure", &listHamiltonian_QM::get_electronic_structure)
      .def("set_electronic_structure", &listHamiltonian_QM::set_electronic_structure)
      .def("compute_overlap", &listHamiltonian_QM::compute_overlap)
      .def("compute_core_Hamiltonian", &listHamiltonian_QM::compute_core_Hamiltonian)
      .def("energy_and_forces", &listHamiltonian_QM::energy_and_forces)

      .def("excite_alp", &listHamiltonian_QM::excite_alp)
      .def("excite_bet", &listHamiltonian_QM::excite_bet)

  ;

//  class_< listHamiltonian_QM >("listHamiltonian_QM")
//      .def(vector_indexing_suite< listHamiltonian_QM >())
//  ;

//  class_< listHamiltonian_QM >("listHamiltonian_QM")
//      .def(vector_indexing_suite< listHamiltonian_QM >())
//  ;



}


#ifdef CYGWIN
BOOST_PYTHON_MODULE(cyghamiltonian_qm){
#else
BOOST_PYTHON_MODULE(libhamiltonian_qm){
#endif

  // Register converters:
  // See here: https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
  //to_python_converter<std::vector<DATA>, VecToList<DATA> >();

  export_hamiltonian_qm_objects();

}




}// namespace libhamiltonian_qm
}// namespace libatomistic
}// liblibra

